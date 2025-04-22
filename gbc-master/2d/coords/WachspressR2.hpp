/*** 
 * @Author: Mao jiaping
 * @Date: 2025-03-01 12:37:11
 * @LastEditTime: 2025-03-01 13:03:17
 * @Description: 
 * @FileName: 
 * @FilePath: /Codes/C++/gbc-master/2d/coords/WachspressR2.hpp
 */

// README:
/*

    Wachspress coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/

#ifndef GBC_WACHSPRESSR2_HPP
#define GBC_WACHSPRESSR2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"
#include "../extra/math.h"

namespace gbc {

    // Wachspress coordinates in R2.
    class WachspressR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        WachspressR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "WachspressR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // M. S. Floater. Generalized barycentric coordinates and applications.
        // Acta Numerica, 24:161-214, 2015 (see Section 5.2).
        // The formula from the paper works only for convex polygons, however here it works for any simple polygon
        // due to the modification proposed in:
        // http://doc.cgal.org/latest/Barycentric_coordinates_2/index.html#Chapter_2D_Generalized_Barycentric_Coordinates (see Section 4.5).
        // Note that this is a more robust but slower O(n^2) version of mean value coordinates. For the fast O(n) version see below.
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const int n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                int i_1 = (i - 1 + n) % n;
                int i_2 = (i + 1 + n) % n;

                double A1 = std::abs(gbc::cross(_v[i_1] - _v[i], _v[i_2] - _v[i]));
                double A2 = std::abs(gbc::cross(_v[i_1] - p, _v[i] - p));
                double A3 = std::abs(gbc::cross(_v[i] - p, _v[i_2] - p));
                b[i] = A1 / (A2 * A3);
                sum += b[i];
            }
            gbc::div(b, sum);
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p) const {

            const int numP = p.size();
#pragma omp parallel for
            for (int i = 0; i < numP; ++i) 
            {
                compute(p[i], p[i].b());
            }
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        int n; // count  of vertex
    };

} // namespace gbc

#endif // GBC_WACHSPRESSR2_HPP
