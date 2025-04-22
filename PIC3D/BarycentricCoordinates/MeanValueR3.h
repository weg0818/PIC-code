/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-19 10:02:54
 * @LastEditTime: 2025-03-06 21:51:58
 * @Description: writed by Mao jiaping, reference to Floater's paper Mean value coordinates in 3d
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/BarycentricCoordinates/MeanValueR3.h
 */

#ifndef GBC_MEANVALUER3_HPP 
#define GBC_MEANVALUER3_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>
// Local includes.
#include "../Basic/VertexR3.hpp"
#include "../Basic/Face.hpp"
#include "BarycentricCoordinatesR3.hpp"

namespace gbc {

    // Mean value coordinates in R3.
    class MeanValueR3 : public BarycentricCoordinatesR3 {

    public:
        // Constructor.
        MeanValueR3(const std::vector<VertexR3> &v, const std::vector<Face> &f, const double tol = 1.0e-10) : super(v, f, tol) 
        {
            n_v = _v.size();
            n_f = _f.size();
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "MeanValueR3";
        }
        
        void compute(const VertexR3 &p, std::vector<double> &b);

        // Compute the coordinates at p using the internal storage from the VertexR3 class.
        inline void compute(VertexR3 &p){
            compute(p, p.b());
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR3> &p, std::vector<std::vector<double> > &bb)
        {
            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) 
            {   
                compute(p[i], bb[i]);
            }
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR3 class.
        void compute(std::vector<VertexR3> &p){

            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) 
            {
                // std::cout << i << std::endl;
                compute(p[i], p[i].b());
            }
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR3> &p) {
            compute(p);
        }
        inline double squaredArea(double a, double b, double c, double p)
        {
            return p * ( p - a ) * ( p - b ) * ( p - c );
        }
    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR3 super;

        void bcInterior(const VertexR3 &p, std::vector<double> &b);

        void get_p2_information(std::vector<VertexR3>& p1, std::vector<std::vector<double>> &p2_coe, std::vector<Face>& face);

        std::vector<VertexR3> p0;
    };

} // namespace gbc

#endif // GBC_MEANVALUER3_HPP
