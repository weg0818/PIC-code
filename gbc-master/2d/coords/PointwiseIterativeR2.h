/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-07 14:22:29
 * @LastEditTime: 2025-02-28 13:55:07
 * @Description: writed by Mao jiaping
 * @FileName: 
 * @FilePath: /Codes/C++/gbc-master/2d/coords/PointwiseIterativeR2.h
 */

#ifndef POINTWISEITERATIVER2_HPP
#define POINTWISEITERATIVER2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // PointwiseIterativeCoordinates in R2.
    class PointwiseIterativeR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        PointwiseIterativeR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) 
        {
            n_v = _v.size();
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "PointwiseIterativeR2";
        }
        
        void compute(const VertexR2 &p, std::vector<double> &b);

        // Compute the coordinates at p using the internal storage from the VertexR2 class.
        inline void compute(VertexR2 &p){
            compute(p, p.b());
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb)
        {
            const int numP = p.size();
            
            bb.resize(numP);
            for (int i = 0; i < numP; ++i) 
            {   
                compute(p[i], bb[i]);
            }
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p){

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

        // isolated point
        void set_iso_p(std::vector<VertexR2> & p, std::vector<std::vector<int>> lines)
        {
            _iso_p = p;
            _iso_p_n = p.size();
            _iso_p_line = lines;
        }

        void compute_iso_pic(std::vector<VertexR2> &p)
        {
            const int numP = p.size();
#pragma omp parallel for           
            for (int i = 0; i < numP; ++i) 
            {
                // std::cout << i << std::endl;
                compute_iso_pic(p[i], p[i].b());
            }
        }

        void compute_iso_pic(const VertexR2 &p, std::vector<double> &b);

        // 对比实验1
        void compute_1(std::vector<VertexR2> &p){
            const int numP = p.size();
#pragma omp parallel for
            for (int i = 0; i < numP; ++i) 
            {
                compute_1(p[i], p[i].b());
            }
        }
        void compute_1(const VertexR2 &p, std::vector<double> &b);
        void bcInterior_1(const VertexR2 &p, std::vector<double> &b);
        // 对比实验2
        void compute_1toN(std::vector<VertexR2> &p){
            const int numP = p.size();
#pragma omp parallel for
            for (int i = 0; i < numP; ++i) 
            {
                compute_1toN(p[i], p[i].b());
            }
        }
        void compute_1toN(const VertexR2 &p, std::vector<double> &b);
        void bcInterior_1toN(const VertexR2 &p, std::vector<double> &b);
        // 对比实验3
        void compute_Nto2N(std::vector<VertexR2> &p){
            const int numP = p.size();
#pragma omp parallel for
            for (int i = 0; i < numP; ++i) 
            {
                compute_Nto2N(p[i], p[i].b());
            }
        }
        void compute_Nto2N(const VertexR2 &p, std::vector<double> &b);
        void bcInterior_Nto2N(const VertexR2 &p, std::vector<double> &b);
        // 对比实验4
        void compute_1to2N(std::vector<VertexR2> &p){
            const int numP = p.size();
#pragma omp parallel for
            for (int i = 0; i < numP; ++i) 
            {
                compute_1to2N(p[i], p[i].b());
            }
        }
        void compute_1to2N(const VertexR2 &p, std::vector<double> &b);
        void bcInterior_1to2N(const VertexR2 &p, std::vector<double> &b);
        // 对比实验5
        void compute_p2(std::vector<VertexR2> &p){
            const int numP = p.size();
#pragma omp parallel for
            for (int i = 0; i < numP; ++i) 
            {
                compute_p2(p[i], p[i].b());
            }
        }
        void compute_p2(const VertexR2 &p, std::vector<double> &b);
        void bcInterior_p2(const VertexR2 &p, std::vector<double> &b);

        // 对比实验6 1/4 1/2 1/4
        void compute_p1p2(std::vector<VertexR2> &p){
            const int numP = p.size();
#pragma omp parallel for
            for (int i = 0; i < numP; ++i) 
            {
                compute_p1p2(p[i], p[i].b());
            }
        }
        void compute_p1p2(const VertexR2 &p, std::vector<double> &b);
        void bcInterior_p1p2(const VertexR2 &p, std::vector<double> &b);

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        void bcInterior(const VertexR2 &p, std::vector<double> &b);
    
        void update_delta(std::vector<VertexR2>& p1, std::vector<double>& p1_coe, std::vector<VertexR2>& p2, std::vector<std::vector<double>>& p2_coe,
            VertexR2& delta, std::vector<double>& delta_coe);
    
        void get_rhoi_information(std::vector<VertexR2>& pi, double alpha,
            VertexR2& delta, std::vector<double>& delta_coe, 
            std::vector<double>& rhoi_coe, std::vector<double>& d_coei);


        void bcInterior_iso(const VertexR2 &p, std::vector<double> &b);
        
        void update_delta_iso(std::vector<VertexR2>& p1, std::vector<double>& p1_coe, 
            std::vector<VertexR2>& p2,std::vector<std::vector<double>> & p2_coe, std::vector<std::vector<int>> & p2_index,
            VertexR2& delta, std::vector<double>& delta_coe);

        inline double f1(double x)
        {
            return x * x * x + 3 * x * x * (1 - x);
        };

        std::vector<VertexR2> _iso_p;
        std::vector<std::vector<int>> _iso_p_line;
        int _iso_p_n = 0;

    };

} // namespace gbc

#endif // GBC_POINTWISEITERATIVER2_HPP
