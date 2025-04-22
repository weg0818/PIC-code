#include "../extra/math.h"
#include "PositiveMeanValueR2.h"
#include "MeanValueR2.hpp"

using namespace gbc;

void gbc::PositiveMeanValueR2::compute(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior(p, b);
    return;
}

void gbc::PositiveMeanValueR2::bcInterior(const VertexR2 &p, std::vector<double> &b)
{
    std::vector<bool> visible;
    std::vector<int> is_new;
    std::vector<VertexR2> new_poly; 
    std::unordered_map<int, std::pair<double, std::pair<int, int>>> new_points;
    VertexR2 pp = p;
    cal_visible_polygon(pp, visible, is_new, new_poly, new_points);

    // mvc
    MeanValueR2 mvc = MeanValueR2(new_poly);
    std::vector<double> bc;
    mvc.compute(pp, bc);

    // pmvc
    for (int i = 0; i < is_new.size(); ++i)
    {
        int index = is_new[i];
        if (index == -1)
        {
            int index1 = new_points[i].second.first;
            int index2 = new_points[i].second.second;
            double t1 =  1 - new_points[i].first;
            double t2 = 1 - t1;
            b[index1] += t1 * bc[i];
            b[index2] += t2 * bc[i];
        }
        else
        {
            b[index] += bc[i];
        }
    }

}

// Find the visible polygon about the point p
// visible: Marks vertices as visible or invisible, 1 is visible, 0 is invisible.
// is_new: Stores the situation of the vertices of the visible polygon, >0 are all under the table of old vertices, -1 means new vertices.
// new_poly: Stores the vertices of the visible polygon
// new_points: Store the added vertices with structure {added vertices: make_pair(t, make_pair(vertex1 subscript, vertex2 subscript))}, added vertices = (1-t) * vertex1 + t * vertex2
void gbc::PositiveMeanValueR2::cal_visible_polygon(VertexR2 & p, std::vector<bool> & visible, 
    std::vector<int> & is_new, std::vector<VertexR2> & new_poly, 
    std::unordered_map<int, std::pair<double, std::pair<int, int>>> & new_points)
{
    int n = n_v;
    visible.clear();
    visible.resize(n, true);
    is_new.clear();

    for (int i = 0; i < n; ++i)
    {
        bool flag_intersect = true;
        for (int j = 0; j < n - 2; ++j)
        {
            // Determine if the line between point p and vertex i intersects other edges, (no more edges containing vertex i).
            int p_1 = (i + j + 1) % n;
            int p_2 = (i + j + 2) % n;
            if (is_visible(p, poly[i], p_1, p_2) == false)
            {
                flag_intersect = false;
                break;
            }
        }
        if (flag_intersect == false)
        {
            // Vertex i Invisible
            visible[i] = false;
        }
    }

    // Start looking for additional visual vertices
    int invisible_begin = 0;
    int invisible_end = -1;
    while (visible[invisible_begin] == 0)
    {
        invisible_begin += 1;
    }
    int i = 0;
    int index1, index2, index3, index4, index5, begin;
    while (i < n)
    {
        index1 = invisible_begin % n;
        if (visible[index1] == 1)
        {
            // vertex visibility
            new_poly.push_back(poly[index1]);
            is_new.push_back(index1);
            invisible_begin += 1;
            i += 1;
        }
        else
        {
            // be unable to see the vertex
            invisible_end = invisible_begin;
            index2 = invisible_end % n;
            // Finding the invisible zone
            while (visible[index2] == 0)
            {
                invisible_end += 1;
                index2 = invisible_end % n;
            }
            // Additional Viewing Points
            invisible_begin -= 1;
            index1 = invisible_begin % n;
            begin = invisible_begin + 1;
            while (begin < invisible_end)
            {
                index3 = begin % n;
                index4 = (begin + 1) % n;
                VertexR2 intersect_point;
                bool state1 = is_intersect_line(p, poly[index1], poly[index3], poly[index4], intersect_point);
                if (state1 == true)
                {
                    // Straight lines have intersections
                    double u;
                    bool state2 = judge_point_on_lines_segment(poly[index3], poly[index4], intersect_point, u);
                    if (state2 != false)
                    {
                        // The intersection point is on the side of the polygon
                        bool flag_intersect;
                        for (int j = 0; j < n - 1; ++j)
                        {
                            int p_1 = (index4 + j) % n;
                            int p_2 = (index4 + j + 1) % n;
                            flag_intersect = true;
                            if (is_visible(p, intersect_point, p_1, p_2) == false)
                            {
                                flag_intersect == false;
                                break;
                            }
                        }
                        if (flag_intersect == true)
                        {
                            new_poly.push_back(intersect_point);
                            index5 = new_poly.size() - 1;
                            new_points[index5] = std::make_pair(u, std::make_pair(index3, index4));
                            is_new.push_back(-1);
                            break;
                        }
                    }
                }
                begin += 1;
            }

            begin = invisible_begin;
            while (begin < invisible_end - 1)
            {
                index3 = begin % n;
                index4 = (begin + 1) % n;
                VertexR2 intersect_point;
                bool state1 = is_intersect_line(p, poly[index2], poly[index3], poly[index4], intersect_point);
                if (state1 == true)
                {
                    // Straight lines have intersections
                    double u;
                    bool state2 = judge_point_on_lines_segment(poly[index3], poly[index4], intersect_point, u);
                    if (state2 != false)
                    {
                        // The intersection point is on the side of the polygon
                        bool flag_intersect;
                        for (int j = 0; j < n - 1; ++j)
                        {
                            int p_1 = (index4 + j) % n;
                            int p_2 = (index4 + j + 1) % n;
                            flag_intersect = true;
                            if (is_visible(p, intersect_point, p_1, p_2) == false)
                            {
                                flag_intersect == false;
                                break;
                            }
                        }
                        if (flag_intersect == true)
                        {
                            new_poly.push_back(intersect_point);
                            index5 = new_poly.size() - 1;
                            new_points[index5] = std::make_pair(u, std::make_pair(index3, index4));
                            is_new.push_back(-1);
                            break;
                        }
                    }
                }
                begin += 1;
            }
            int index = new_poly.size();
            if (index > 1 && same_point(new_poly[index - 1], new_poly[index - 2]))
            {
                new_poly.pop_back();
                new_poly.pop_back();
                new_points.erase(index - 1);
                new_points.erase(index - 2);
                is_new.pop_back();
                is_new.pop_back();
            }
            i = i + invisible_end - invisible_begin - 1;
            invisible_begin = invisible_end;
        }
    }


    // de-duplicate
    std::vector<VertexR2> nnew_poly;
    std::vector<int> is_nnew;
    std::unordered_map<int, std::pair<double, std::pair<int, int>>> nnew_points;
    remove_same(new_poly, is_new, new_points, nnew_poly, is_nnew, nnew_points);
    new_points = nnew_points;
    is_new = is_nnew;
    new_points = nnew_points;
}


 // Determine whether the intersection point and other edges have intersections, do not judge the edges that contain p-points, do not judge the judgement that contains a known intersection point
bool gbc::PositiveMeanValueR2::is_visible(VertexR2 & p, VertexR2 & intersect_point, int p_1, int p_2)
{
    int n = n_v;
    if (is_intersect_segment(p, intersect_point, poly[p_1], poly[p_2]) == true)
    {
        // If intersect, one case is to intersect the endpoints, to determine whether it is through the endpoints, through the endpoints also has two cases
        // If they don't intersect, continue traversing
        VertexR2 ip;
        bool state1 = is_intersect_line(p, intersect_point, poly[p_1], poly[p_2], ip);
        if (same_point(ip, poly[p_1]) == true)
        {
            // Through the vertex.
            if (check_ip_visible(p, intersect_point, poly[(p_1 - 1) % n], poly[p_2]) == true)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else if (same_point(ip, poly[p_2]) == true)
        {
            // Through the vertex.
            if (check_ip_visible(p, intersect_point, poly[p_1], poly[(p_2 + 1) % n]) == true)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    return true;
}

bool gbc::PositiveMeanValueR2::check_ip_visible(VertexR2 & p, VertexR2 & ip, VertexR2 & i_1, VertexR2 & i_2)
{
    VertexR2 a = ip - p;
    VertexR2 b = i_1 - p;
    VertexR2 c = i_2 - p;
    double tmp1 = cross(a, b);
    double tmp2 = cross(a, c);
    double tmp = tmp1 * tmp2;
    if (tmp >= -(_tol))
    {
        return true;
    }        
    return false;
}

void gbc::PositiveMeanValueR2::remove_same(std::vector<VertexR2> & old_poly, std::vector<int> & old_is_new, 
    std::unordered_map<int, std::pair<double, std::pair<int, int>>> & old_points,
    std::vector<VertexR2> & new_poly, std::vector<int> & is_new,
    std::unordered_map<int, std::pair<double, std::pair<int, int>>> & new_points)
{
    new_poly.clear();
    is_new.clear();

    int begin = 0;
    // VertexR2 a = old_poly[old_poly.size() - 1] - old_poly[begin];
    VertexR2 a = old_poly[begin - 1] - old_poly[begin];
    VertexR2 b = old_poly[begin + 1] - old_poly[begin];
    while (abs(cross(a, b)) < _tol)
    {
        begin += 1;
        a[0] = old_poly[begin - 1][0] - old_poly[begin][0];
        a[1] = old_poly[begin - 1][1] - old_poly[begin][1];
        b[0] = old_poly[begin + 1][0] - old_poly[begin][0];
        b[1] = old_poly[begin + 1][1] - old_poly[begin][1];
    }
    int n = old_poly.size();
    for (int i = 0; i < n; ++i)
    {
        int index1 = (begin + i - 1) % n;
        int index2 = (begin + i) % n;
        int index3 = (begin + i + 1) % n;
        a[0] = old_poly[index1][0] - old_poly[index2][0];
        a[1] = old_poly[index1][1] - old_poly[index2][1];
        b[0] = old_poly[index3][0] - old_poly[index2][0];
        b[1] = old_poly[index3][1] - old_poly[index2][1];

        if (cross(a, b) != 0)
        {
            new_poly.push_back(old_poly[index2]);
            is_new.push_back(old_is_new[i]);
            int index = new_poly.size() - 1;
            if (old_is_new[i] == -1)
            {
                new_points[index] = old_points[index];
            }
        }   
    }
}