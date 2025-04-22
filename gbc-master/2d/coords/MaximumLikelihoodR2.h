/*** 
 * @Author: Chang qing jun
 * @Date: 2024-09-12 16:02:17
 * @LastEditTime: 2024-10-10 17:02:26
 * @Description: writed by Chang qing jun
 * @FileName: 
 * @FilePath: /Codes/C++/gbc-master/2d/coords/MaximumLikelihoodR2.h
 */

#ifndef GBC_MAXIMUMLIKEHOOD_H
#define GBC_MAXIMUMLIKEHOOD_H

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
    class MaximumLikelihoodR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        MaximumLikelihoodR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10, int k = 0) : super(v, tol) 
        {
            n_v = _v.size();

            poly = Eigen::MatrixXd::Zero(2, _v.size());
            for (int i = 0; i < _v.size(); ++i)
            {
                poly(0, i) = _v[i][0];
                poly(1, i) = _v[i][1];
            }

        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "MaximumLikelihoodR2";
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


    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        Eigen::VectorXd dmlCoordinates(Eigen::Vector2d v, Eigen::MatrixXd & polygon);

        /// @brief Solve convex optimization problem
        /// @param f the objective function of convex optimization
        /// @param f1 the first derivative of the objective function
        /// @param f2 the Hessian matrix function of objective function
        /// @param Phi extreme point
        /// @param error error
        /// @param rho a parameter in Line Search Method
        /// @param sigma a parameter in Line Search Method
        /// @param M a parameter in Line Search Method
        /// @return true/false
        bool convex_opt(std::function<double(Eigen::Vector2d)> f,
                        std::function<Eigen::Vector2d(Eigen::Vector2d)> f1,
                        std::function<Eigen::Matrix2d(Eigen::Vector2d)> f2,
                        Eigen::Vector2d *Phi, double error = 1E-10, double rho = 0.55,
                        double sigma = 0.4, int M = 20);
        Eigen::MatrixXd poly;
    };



} // namespace gbc

#endif // GBC_MAXIMUMLIKEHOOD_H
