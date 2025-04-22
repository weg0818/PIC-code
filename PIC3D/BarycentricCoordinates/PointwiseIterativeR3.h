/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-08 14:59:48
 * @LastEditTime: 2025-02-24 12:55:05
 * @Description: writed by Mao jiaping
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/BarycentricCoordinates/PointwiseIterativeR3.h
 */

#ifndef GBC_POINTWISEITERATIVER3_HPP
#define GBC_POINTWISEITERATIVER3_HPP
 
// STL includes.
#include <vector>
#include <cassert>
#include <cmath>
#include <chrono>

// Local includes.
#include "../Basic/VertexR3.hpp"
#include "../Basic/Face.hpp"
#include "BarycentricCoordinatesR3.hpp"

namespace gbc {

    // Mean value coordinates in R3.
    class PointwiseIterativeR3 : public BarycentricCoordinatesR3 {

    public:
        // Constructor.
        PointwiseIterativeR3(const std::vector<VertexR3> &v, const std::vector<Face> &f, const double tol = 1.0e-10) : super(v, f, tol) 
        {
            n_v = _v.size();
            n_f = _f.size();
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "PointwiseIterativeR3";
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
// #pragma omp parallel for  
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

        // 孤立点
        void set_iso_p(std::vector<VertexR3> & p)
        {
            _iso_p = p;
            _iso_p_n = p.size();
            
        }

        void compute_iso_pic(std::vector<VertexR3> &p)
        {
            const int numP = p.size();
            for (int i = 0; i < numP; ++i) 
            {
                // std::cout << i << std::endl;
                compute_iso_pic(p[i], p[i].b());
            }
        }

        void compute_iso_pic(const VertexR3 &p, std::vector<double> &b);

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR3 super;

        void bcInterior(const VertexR3 &p, std::vector<double> &b);
    
        void update_delta(std::vector<VertexR3>& p1, std::vector<double>& p1_coe, std::vector<VertexR3>& p2, std::vector<std::vector<double>>& p2_coe,
            std::vector<Face>& face, VertexR3& delta, std::vector<double>& delta_coe);
    
        void get_rhoi_information(std::vector<VertexR3>& pi, double alpha,
            VertexR3& delta, std::vector<double>& delta_coe, 
            std::vector<double>& rhoi_coe, std::vector<double>& di_coe);

        void adjustSequenceOfFase(std::vector<VertexR3>& p, std::vector<Face>& face);
        
        bool get_p2_information(std::vector<VertexR3>& p1, std::vector<VertexR3> &p2, 
            std::vector<std::vector<double>> &p2_coe, std::vector<Face>& face);

        void bcInterior_iso(const VertexR3 &p, std::vector<double> &b);

        
        bool isZero(VertexR3 & a);

        inline double f1(double x)
        {
            return x * x * x + 3 * x * x * (1 - x);
        };

        inline double squaredArea(double a, double b, double c, double p)
        {
            return p * ( p - a ) * ( p - b ) * ( p - c );
        }

        std::vector<VertexR3> _iso_p;
        int _iso_p_n = 0;
    };

} // namespace gbc

#endif // GBC_POINTWISEITERATIVER3_HPP
