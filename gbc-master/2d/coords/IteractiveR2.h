/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-19 12:52:55
 * @LastEditTime: 2024-10-10 19:07:10
 * @Description: wirted by Mao jiaping, reference to Chongyang Deng's paper Iterative coordinates
 * @FileName: 
 * @FilePath: /Codes/C++/gbc-master/2d/coords/IteractiveR2.h
 */

#ifndef GBC_ITERACTIVER2_HPP
#define GBC_ITERACTIVER2_HPP

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
    class IteractiveR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        IteractiveR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10, int k = 0) : super(v, tol) 
        {
            n_v = _v.size();
            _k = k;
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "IteractiveR2";
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

#pragma omp parallel for
            for (size_t i = 0; i < numP; ++i) 
            {
                compute(p[i], p[i].b());
            }
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

        void compute_iteractive_times(std::vector<VertexR2> &p){

            const size_t numP = p.size();
#pragma omp parallel for
            for (size_t i = 0; i < numP; ++i) 
            {
                // std::cout << "Find the number of iterations: " << i << std::endl;
                compute_iteractive_times(p[i], p[i].b());
            }
            std::cout << "The number of iterations should be: " << _k << std::endl;
        }

        void compute_iteractive_times(const VertexR2 &p, std::vector<double>  &bb);

        void set_k(int k)
        {
            _k = k;
        }

        void set_estimate_k()
        {
            _k = round(2 / (gbc::PI * gbc::PI) * n_v * n_v * log(n_v + 1));
            std::cout << "The number of iterations is estimated as: " << _k << std::endl;
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        void bcInterior(const VertexR2 &p, std::vector<double> &b);
        void bcInterior_i(const VertexR2 &p);

        void get_M_k(Eigen::MatrixXd & p_i, Eigen::MatrixXd & M_k);
        void MatrixXdToVecVertex(Eigen::MatrixXd & p_i, std::vector<VertexR2> & poly);
        int _k = 0;
    };

} // namespace gbc

#endif // GBC_ITERACTIVER2_HPP
