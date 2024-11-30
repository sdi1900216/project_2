#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Delaunay_mesher_2.h>
#include "triangulation.h"
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/centroid.h>
#include <CGAL/convex_hull_2.h>
#include <unordered_map>
#include <cmath>

//////////////////////////////////////////////////////////

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef K::Line_2 Line;
typedef CDT::Point Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CDT::Face_handle Face_handle;
using namespace std;
///////////////////////////////////////////////////

int is_obtuse_angle(Point &A, Point &B, Point &C)
{
    return angle(A, B, C) == CGAL::OBTUSE;
}

int count_Obtuse_Angles(CDT &cdt)
{
    int count = 0;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point a = fit->vertex(0)->point();
        Point b = fit->vertex(1)->point();
        Point c = fit->vertex(2)->point();
        // cout << "Checking triangle with points: (" << a.x() << ", " << a.y() << "), (" << b.x() << ", " << b.y() << "), (" << c.x() << ", " << c.y() << ")" << endl;

        if (is_obtuse_angle(a, b, c) || is_obtuse_angle(b, c, a) || is_obtuse_angle(c, a, b))
        {
            count++;
            // cout << "Obtuse angle found at vertex " << is_Obtuse << endl;
        }
    }
    return count;
}

Point project_point(Point &A, Point &B, Point &P)
{ // επιστρέφει την προβολή από ένα σημείο P στην πλευρά που σχηματίζουν τα Α-Β
    Line line(A, B);
    return line.projection(P);
}
///////////////////////////////////////////////////////////////

Polygon_2 find_convex_polygon_around_obtuse_triangle(CDT &cdt, Face_handle face)
{
    set<Face_handle> visited_faces;
    queue<Face_handle> face_queue;
    Polygon_2 convex_polygon;
    // Αρχικοποίηση με το αρχικό αμβλυγώνιο τρίγωνο
    face_queue.push(face);
    visited_faces.insert(face);

    while (!face_queue.empty())
    {
        Face_handle current_face = face_queue.front();
        face_queue.pop();

        // Προσθήκη σημείων του current_face στο polygon
        for (int i = 0; i < 3; i++)
        {
            convex_polygon.push_back(current_face->vertex(i)->point());
        }
        // Έλεγχος για τους γειτονικούς κόμβους
        for (int i = 0; i < 3; i++)
        {
            Face_handle neighbor = current_face->neighbor(i);
            if (!cdt.is_infinite(neighbor) && visited_faces.find(neighbor) == visited_faces.end())
            {
                Point a = neighbor->vertex(0)->point();
                Point b = neighbor->vertex(1)->point();
                Point c = neighbor->vertex(2)->point();

                if (is_obtuse_angle(a, b, c) || is_obtuse_angle(b, c, a) || is_obtuse_angle(c, a, b)) // Έλεγχος αν το γειτονικό τρίγωνο είναι αμβλυγώνιο
                {
                    face_queue.push(neighbor);
                    visited_faces.insert(neighbor);
                }
            }
        }
    }
    // Δημιουργία του κυρτού περιβλήματος των σημείων του πολυγώνου
    Polygon_2 convex_hull;
    convex_hull_2(convex_polygon.vertices_begin(), convex_polygon.vertices_end(), back_inserter(convex_hull));
    return convex_hull;
}

// Μέθοδος εισαγωγής Steiner points στο εσωτερικό κυρτών πολυγώνων που σχηματίζονται από αμβλυγώνια τρίγωνα
Point insert_Steiner_point_in_convex_polygons(CDT &cdt, Polygon_2 &region_boundary)
{
    int steiner_counter = 0;

    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point a = fit->vertex(0)->point();
        Point b = fit->vertex(1)->point();
        Point c = fit->vertex(2)->point();
        // Ελέγχουμε αν το τρίγωνο έχει αμβλεία γωνία
        if (is_obtuse_angle(a, b, c) || is_obtuse_angle(b, c, a) || is_obtuse_angle(c, a, b))
        {
            // Δημιουργία κυρτού πολυγώνου γύρω από το αμβλυγώνιο τρίγωνο και τους γείτονές του
            Polygon_2 convex_polygon = find_convex_polygon_around_obtuse_triangle(cdt, fit);

            // Υπολογισμός του centroid του κυρτού πολυγώνου
            Point centroid = CGAL::centroid(convex_polygon.vertices_begin(), convex_polygon.vertices_end());
            return centroid;
        }
    }
}

///////////////////////////////////////////////////////////

// Συνάρτηση που επιστρέφει σημείο Steiner για μία από τις 5 στρατηγικές
Point select_steiner_point(Point &a, Point &b, Point &c, int strategy, CDT &cdt, Polygon_2 convex_hull)
{
    switch (strategy)
    {
    case 0:
    { // Περίκεντρο check
        return circumcenter(a, b, c);
    }
    case 1:
    { // Κέντρο βάρους ενός τριγώνου
        double cx = (a.x() + b.x() + c.x()) / 3.0;
        double cy = (a.y() + b.y() + c.y()) / 3.0;
        return Point(cx, cy);
    }
    case 2:
    { // Μέσο της μεγαλύτερης ακμής check
        double d_ab = squared_distance(a, b);
        double d_bc = squared_distance(b, c);
        double d_ca = squared_distance(c, a);
        if (d_ab >= d_bc && d_ab >= d_ca)
        {
            return midpoint(a, b);
        }
        else if (d_bc >= d_ab && d_bc >= d_ca)
        {
            return midpoint(b, c);
        }
        else
        {
            return midpoint(c, a);
        }
    }
    case 3:
    { // Προβολή της κορυφής της αμβλείας γωνίας στην απέναντι πλευρά
        if (is_obtuse_angle(b, a, c))
        {
            return project_point(b, c, a); // προβολή του A στην πλευρά B-C
        }
        else if (is_obtuse_angle(a, b, c))
        {
            return project_point(a, c, b); // του B στην πλευρά A-C
        }
        else if (is_obtuse_angle(a, c, b))
        {
            return project_point(a, b, c); // και του C στην πλευρά AB
        }
    }
    case 4:
    {
        return insert_Steiner_point_in_convex_polygons(cdt, convex_hull);
    }
    default:
        throw invalid_argument("Invalid strategy selected.");
    }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////

// State definition
struct State
{
    CDT cdt;
    int obtuse_count;
    int steiner_points;
    std::vector<Point> steiner_locations;
    std::vector<int> strategies;

    // Εδώ χρησιμοποιούμε την compareStates για τη σύγκριση
    bool operator==(const State &other) const
    {
        return obtuse_count == other.obtuse_count &&
               steiner_points == other.steiner_points &&
               steiner_locations == other.steiner_locations &&
               strategies == other.strategies;
    }
};

// Συνάρτηση για τη σύγκριση δύο καταστάσεων State
bool compareStates(const State &a, const State &b)
{
    return a.obtuse_count == b.obtuse_count &&
           a.steiner_points == b.steiner_points &&
           a.steiner_locations == b.steiner_locations &&
           a.strategies == b.strategies;
}

struct StateHash
{
    std::size_t operator()(const State &s) const
    {
        size_t hash_val = 0;
        // Hash για obtuse_count, steiner_points, και άλλα χαρακτηριστικά
        hash_val ^= std::hash<int>{}(s.obtuse_count) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        hash_val ^= std::hash<int>{}(s.steiner_points) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        for (const auto &loc : s.steiner_locations)
        {
            hash_val ^= std::hash<double>{}(loc.x()) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            hash_val ^= std::hash<double>{}(loc.y()) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        }
        return hash_val;
    }
};

State bfs_triangulation(CDT &initial_cdt, Polygon_2 &convex_hull, int &best_obtuse, CDT &best_cdt, int max_depth)
{
    std::queue<State> queue;
    std::unordered_set<State, StateHash> visited; // Χρησιμοποιούμε custom hash για State

    // Αρχικοποίηση με την αρχική κατάσταση
    int initial_obtuse = count_Obtuse_Angles(initial_cdt);
    State initial_state = {initial_cdt, initial_obtuse, 0, {}, {}};
    queue.push(initial_state);
    visited.insert(initial_state);
    cout << "HERE" << endl;
    best_obtuse = initial_obtuse;
    best_cdt = initial_cdt;

    // Εξερεύνηση μέσω BFS
    while (!queue.empty())
    {
        State current_state = queue.front();
        queue.pop();

        // Αν η τρέχουσα κατάσταση είναι βέλτιστη, ενημερώνουμε τη βέλτιστη λύση
        if (current_state.obtuse_count < best_obtuse)
        {
            best_obtuse = current_state.obtuse_count;
            best_cdt = current_state.cdt;
            // cout << "HERE" << endl;
        }

        // Αν φτάσουμε στο μέγιστο βάθος ή δεν έχουμε άλλες αμβλείες γωνίες, σταματάμε
        if (current_state.steiner_points >= max_depth || current_state.obtuse_count == 0)
            return current_state;

        // Εξερεύνηση όλων των τριγώνων με αμβλείες γωνίες
        for (auto fit = current_state.cdt.finite_faces_begin(); fit != current_state.cdt.finite_faces_end(); ++fit)
        {
            Point a = fit->vertex(0)->point();
            Point b = fit->vertex(1)->point();
            Point c = fit->vertex(2)->point();

            if (CGAL::angle(a, b, c) == CGAL::OBTUSE || CGAL::angle(b, c, a) == CGAL::OBTUSE || CGAL::angle(c, a, b))
            {
                // Δοκιμή όλων των στρατηγικών
                for (int strategy = 0; strategy < 5; ++strategy)
                {
                    Point steiner = select_steiner_point(a, b, c, strategy, current_state.cdt, convex_hull);

                    // Έλεγχος αν το σημείο είναι μέσα στο κυρτό περίβλημα
                    if (convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDED_SIDE || convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDARY)
                    {
                        CDT temp_cdt = current_state.cdt;
                        temp_cdt.insert(steiner);
                        // cout << "HERE4" << endl;
                        int new_obtuse = count_Obtuse_Angles(temp_cdt);

                        // Δημιουργούμε μια νέα κατάσταση και ελέγχουμε αν υπάρχει ήδη
                        State new_state = {temp_cdt, new_obtuse, current_state.steiner_points + 1, {}, {}};

                        // Αν η νέα κατάσταση δεν έχει επισκεφθεί ξανά, την προσθέτουμε
                        if (visited.find(new_state) == visited.end())
                        {
                            queue.push(new_state);
                            visited.insert(new_state);
                            // cout << "HERE--5" << endl;
                        }
                    }
                }
            }
        }
    }
}

// Κύρια συνάρτηση
void triangulate(const vector<int> &points_x, const vector<int> &points_y, const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints)
{
    CDT cdt;
    vector<Point> points;

    for (size_t i = 0; i < points_x.size(); ++i)
    {
        points.push_back(Point(points_x[i], points_y[i]));
    }

    Polygon_2 convex_hull;
    for (size_t i : region_boundary)
    {
        convex_hull.push_back(points[i]);
    }

    // Προσθήκη constraints
    for (size_t i = 0; i < region_boundary.size(); ++i)
    {
        int next = (i + 1) % region_boundary.size();
        cdt.insert_constraint(points[region_boundary[i]], points[region_boundary[next]]);
    }
    for (const auto &constraint : additional_constraints)
    {
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }

    int best_obtuse = count_Obtuse_Angles(cdt);
    cout << "Initial obtuse angles: " << best_obtuse << endl;
    CDT best_cdt;
    int max_depth = 1000;

    State best = bfs_triangulation(cdt, convex_hull, best_obtuse, best_cdt, max_depth);
    cout << "Final obtuse angles: " << best.obtuse_count << endl;
    CGAL::draw(best.cdt);
}

//////////////////////////////////////////////////////////////////////////
/*
//////////////////////////////// τριγωνοποίηση
void triangulate(const vector<int> &points_x, const vector<int> &points_y, const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints)
{
   CDT cdt;
   // προσθήκη σημείων και ακμών του γράφου
   vector<Point> points;
   for (size_t i = 0; i < points_x.size(); ++i)
   {
       points.push_back(Point(points_x[i], points_y[i]));
   }
   // Δημιουργία του πολυγώνου του κυρτού περιβλήματος
   Polygon_2 convex_hull;
   for (size_t i : region_boundary)
   {
       convex_hull.push_back(points[i]);
   }

   for (int i = 0; i < points_x.size(); i++)
       cout << points_x[i] << " ";
   for (int i = 0; i < points_y.size(); i++)
       cout << points_y[i] << " ";

   // προσθήκη των ακμών του κυρτού περιβλήματος για όριο
   for (size_t i = 0; i < region_boundary.size(); ++i)
   {
       int next = (i + 1) % region_boundary.size();
       cdt.insert_constraint(points[region_boundary[i]], points[region_boundary[next]]);
   }

   // προσθήκη πρόσθετων constraints
   for (const auto &constraint : additional_constraints)
   {
       cdt.insert_constraint(points[constraint.first], points[constraint.second]);
   }

   int steiner_counter = 0;

   int best_obtuse_count = count_Obtuse_Angles(cdt); // μετρητής για το πόσες αμβλείες γωνίες υπάρχουν
   int original_graph_count = best_obtuse_count;
   CDT best_cdt = cdt;
   cout << "Number of obtuse angles in the triangulation (before using any methods):" << original_graph_count << endl;
   if (best_obtuse_count > 0)
   {
       bool found_steiner_point = true;
       while (found_steiner_point)
       {
           found_steiner_point = false; // aν δεν βρεθεί "καλό σημείο" για να βάλουμε steiner, το loop θα σταματήσει
           // Κάνουμε iterate στα finite faces της τριγωνοποίησης
           for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
           {
               Point a = fit->vertex(0)->point();
               Point b = fit->vertex(1)->point();
               Point c = fit->vertex(2)->point();
               int best_strategy = 0;
               for (int strategy = 0; strategy < 5; strategy++)
               {
                   Point steiner = select_steiner_point(a, b, c, strategy, cdt, convex_hull);
                   if (convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDED_SIDE || convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDARY)
                   {
                       CDT temp_cdt = cdt;
                       temp_cdt.insert(steiner);
                       int new_obtuse_count = count_Obtuse_Angles(temp_cdt);
                       if (new_obtuse_count < best_obtuse_count)
                       {
                           found_steiner_point = true;
                           best_obtuse_count = new_obtuse_count;
                           best_cdt = temp_cdt;
                           best_strategy = strategy;
                           draw(temp_cdt);
                       }
                   }
               }
               if (original_graph_count > best_obtuse_count)
               {
                   cdt = best_cdt;
                   original_graph_count = best_obtuse_count;
                   steiner_counter++;
                   cout << "Best strategy was:" << best_strategy << endl;
                   break;
               }
           }
       }

       cout << "Number of obtuse angles in the triangulation: " << best_obtuse_count << endl;
       cout << "Number of steiner points added in the triangulation: " << steiner_counter << endl;
       draw(cdt);
   }
}
*/