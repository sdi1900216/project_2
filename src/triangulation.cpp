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
    vector<Point> steiner_locations;
    vector<int> strategies;

    bool operator==(const State &other) const
    {
        return obtuse_count == other.obtuse_count &&
               steiner_points == other.steiner_points &&
               steiner_locations == other.steiner_locations &&
               strategies == other.strategies;
    }
};

void printStateDetails(State &state)
{
    cout << "Obtuse Count: " << state.obtuse_count << endl;
    cout << "Steiner Points: " << state.steiner_points << endl;
    cout << endl;
}
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
    size_t operator()(const State &s) const
    {
        size_t hash_val = 0;
        // Hash για obtuse_count, steiner_points, και άλλα χαρακτηριστικά
        hash_val ^= hash<int>{}(s.obtuse_count) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        hash_val ^= hash<int>{}(s.steiner_points) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        for (const auto &loc : s.steiner_locations)
        {
            hash_val ^= hash<double>{}(loc.x()) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            hash_val ^= hash<double>{}(loc.y()) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        }
        return hash_val;
    }
};

State bfs_triangulation(CDT &initial_cdt, Polygon_2 &convex_hull, int &best_obtuse, CDT &best_cdt, int max_iterations)
{
    queue<State> queue;
    unordered_set<State, StateHash> visited; // Χρησιμοποιούμε custom hash για State
    // Αρχικοποίηση με την αρχική κατάσταση
    State initial_state = {initial_cdt, count_Obtuse_Angles(initial_cdt), 0, {}, {}};
    State best_state = initial_state;
    queue.push(initial_state);
    visited.insert(initial_state);
    int iteration_count = 0;
    best_cdt = initial_cdt;

    // Εξερεύνηση μέσω BFS
    while (!queue.empty() && iteration_count < max_iterations && best_state.obtuse_count > 0)
    {
        State current_state = queue.front();
        queue.pop();
        // Αν η τρέχουσα κατάσταση είναι βέλτιστη, ενημερώνουμε τη βέλτιστη λύση
        if (current_state.obtuse_count < best_state.obtuse_count)
        {
            best_cdt = current_state.cdt;
            best_state = current_state;
            iteration_count -= 10; // Επαναφορά του μετρητή επαναλήψεων επειδή βελτιώθηκε
        }
        // Αν φτάσουμε στο μέγιστο βάθος ή δεν έχουμε άλλες αμβλείες γωνίες, σταματάμε
        if (iteration_count >= max_iterations || best_state.obtuse_count == 0)
            return best_state;
        // Εξερεύνηση όλων των τριγώνων με αμβλείες γωνίες
        for (auto fit = current_state.cdt.finite_faces_begin(); fit != current_state.cdt.finite_faces_end(); ++fit)
        {
            Point a = fit->vertex(0)->point();
            Point b = fit->vertex(1)->point();
            Point c = fit->vertex(2)->point();

            // if (is_obtuse_angle(a, b, c) || is_obtuse_angle(b, c, a) || is_obtuse_angle(c, a, b))
            //{
            //  Δοκιμή όλων των στρατηγικών
            for (int strategy = 0; strategy < 5; ++strategy)
            {
                Point steiner = select_steiner_point(a, b, c, strategy, current_state.cdt, convex_hull);
                // Έλεγχος αν το σημείο είναι μέσα στο κυρτό περίβλημα
                if (convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDED_SIDE || convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDARY)
                {
                    CDT temp_cdt = current_state.cdt;
                    temp_cdt.insert(steiner);
                    int new_obtuse = count_Obtuse_Angles(temp_cdt);
                    // Δημιουργούμε μια νέα κατάσταση και ελέγχουμε αν υπάρχει ήδη
                    State new_state = {temp_cdt, new_obtuse, current_state.steiner_points + 1, {}, {}};
                    if (new_obtuse <= initial_state.obtuse_count)
                    {
                        // Αν η νέα κατάσταση δεν έχει επισκεφθεί ξανά, την προσθέτουμε
                        if (visited.find(new_state) == visited.end())
                        {
                            queue.push(new_state);
                            visited.insert(new_state);
                        }
                    }
                    else
                    {
                        visited.insert(new_state);
                    }
                }
            }
            //}
            if (current_state.obtuse_count >= best_state.obtuse_count)
            {
                iteration_count++;
            }
        }
    }
    return best_state;
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
    int max_depth = 20000;

    State best = bfs_triangulation(cdt, convex_hull, best_obtuse, best_cdt, max_depth);
    printStateDetails(best);
    draw(best.cdt);
    TriangulationResult results;
    results.obtuse_count = best.obtuse_count;
    for (CDT::Finite_edges_iterator eit = best.cdt.finite_edges_begin();
         eit != best.cdt.finite_edges_end(); ++eit)
    {
        CDT::Face_handle face = eit->first;
        int index = eit->second;
        CDT::Vertex_handle v1 = face->vertex((index + 1) % 3);
        CDT::Vertex_handle v2 = face->vertex((index + 2) % 3);

        // Find the indices of these points in the original input
        auto it1 = find(points.begin(), points.end(), v1->point());
        auto it2 = find(points.begin(), points.end(), v2->point());

        if (it1 != points.end() && it2 != points.end())
        {
            int index1 = distance(points.begin(), it1);
            int index2 = distance(points.begin(), it2);
            results.edges.push_back({index1, index2});
        }
    }
    // Collect Steiner points
    set<Point> original_points(points.begin(), points.end());
    for (CDT::Finite_vertices_iterator vit = best.cdt.finite_vertices_begin();
         vit != best.cdt.finite_vertices_end(); ++vit)
    {
        Point p = vit->point();

        // If the point is not in the original set, it's a Steiner point
        if (original_points.find(p) == original_points.end())
        {
            results.steiner_points_x.push_back(CGAL::to_double(p.x()));
            results.steiner_points_y.push_back(CGAL::to_double(p.y()));
        }
    }
    return results;
}
//////////////////////////////////////////////////////////////////////////
