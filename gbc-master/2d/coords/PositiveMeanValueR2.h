/*** 
 * @Author: Mao jiaping
 * @Date: 2024-10-09 15:29:23
 * @LastEditTime: 2024-11-09 13:41:59
 * @Description: 
 * @FileName: 
 * @FilePath: /Codes/C++/gbc-master/2d/coords/PositiveMeanValueR2.h
 */

#ifndef GBC_POSITIVEMEANVALUER2_H
#define GBC_POSITIVEMEANVALUER2_H

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>
// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"
#include "../extra/math.h"

// Libs.
#include "../extra/libs/eigen/Eigen/Core"
#include "../extra/libs/eigen/Eigen/Dense"

namespace gbc {

    // Mean value coordinates in R3.
    class PositiveMeanValueR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        PositiveMeanValueR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10, int k = 0) : super(v, tol) 
        {
            n_v = _v.size();
            poly = v;
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "PositiveMeanValueR2";
        }
        
        void compute(const VertexR2 &p, std::vector<double> &b);

        // Compute the coordinates at p using the internal storage from the VertexR2 class.
        inline void compute(VertexR2 &p){
            compute(p, p.b());
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb)
        {
            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) 
            {   
                compute(p[i], bb[i]);
            }
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p){

            const size_t numP = p.size();
// #pragma omp paralel for
            for (size_t i = 0; i < numP; ++i) 
            {
                compute(p[i], p[i].b());
            }
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

        // Find the visible polygon about the point p
        // visible: Marks vertices as visible or invisible, 1 is visible, 0 is invisible.
        // is_new: Stores the situation of the vertices of the visible polygon, >0 are all under the table of old vertices, -1 means new vertices.
        // new_poly: Stores the vertices of the visible polygon
        // new_points: Store the added vertices with structure {added vertices: make_pair(t, make_pair(vertex1 subscript, vertex2 subscript))}, added vertices = (1-t) * vertex1 + t * vertex2
        void cal_visible_polygon(VertexR2 & p, std::vector<bool> & visible, std::vector<int> & is_new, 
            std::vector<VertexR2> & new_poly, std::unordered_map<int, std::pair<double, std::pair<int, int>>> & new_points);
    

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        void bcInterior(const VertexR2 &p, std::vector<double> &b);
        
        // Determine whether the intersection point and other edges have intersections, do not judge the edges that contain p-points, do not judge the judgement that contains a known intersection point
        bool is_visible(VertexR2 & p, VertexR2 & intersect_point, int p_1, int p_2);

        bool check_ip_visible(VertexR2 & p, VertexR2 & ip, VertexR2 & i_1, VertexR2 & i_2);

        // de-duplicate
        void remove_same(std::vector<VertexR2> & old_poly, std::vector<int> & old_is_new, std::unordered_map<int, std::pair<double, std::pair<int, int>>> & old_points,
        std::vector<VertexR2> & new_poly, std::vector<int> & is_new, std::unordered_map<int, std::pair<double, std::pair<int, int>>> & new_points);

        std::vector<VertexR2> poly;
    
    };

} // namespace gbc

#endif // GBC_POSITIVEMEANVALUER2_H
