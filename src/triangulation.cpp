#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Delaunay_mesher_2.h>
#include "triangulation.h"
#include "strategies.h"
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
TriangulationResult triangulate(const vector<int> &points_x, const vector<int> &points_y, const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints)
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
    TriangulationResult results;
    results.obtuse_count = best.obtuse_count;
    for (CDT::Finite_edges_iterator eit = best.cdt.finite_edges_begin();
         eit != best.cdt.finite_edges_end(); ++eit) {
        CDT::Face_handle face = eit->first;
        int index = eit->second;
        CDT::Vertex_handle v1 = face->vertex((index+1)%3);
        CDT::Vertex_handle v2 = face->vertex((index+2)%3);
        
        // Find the indices of these points in the original input
        auto it1 = find(points.begin(), points.end(), v1->point());
        auto it2 = find(points.begin(), points.end(), v2->point());
        
        if (it1 != points.end() && it2 != points.end()) {
            int index1 = distance(points.begin(), it1);
            int index2 = distance(points.begin(), it2);
            results.edges.push_back({index1, index2});
        }
    }
 // Collect Steiner points
    set<Point> original_points(points.begin(), points.end());
    for (CDT::Finite_vertices_iterator vit = best.cdt.finite_vertices_begin();
         vit != best.cdt.finite_vertices_end(); ++vit) {
        Point p = vit->point();
        
        // If the point is not in the original set, it's a Steiner point
        if (original_points.find(p) == original_points.end()) {
            results.steiner_points_x.push_back(CGAL::to_double(p.x()));
            results.steiner_points_y.push_back(CGAL::to_double(p.y()));
        }
    }
    return results;
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