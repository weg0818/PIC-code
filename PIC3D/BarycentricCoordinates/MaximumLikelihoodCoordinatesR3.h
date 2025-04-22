/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-21 15:40:17
 * @LastEditTime: 2024-10-10 20:43:22
 * @Description: From the source https://github.com/DanStroeter/CageModeler
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/BarycentricCoordinates/MaximumLikelihoodCoordinatesR3.h
 */
#ifndef GBC_MAXIMUMLIKELIHOOD_H
#define GBC_MAXIMUMLIKELIHOOD_H

// includes.
#include <vector>
#include <cassert>
#include <cmath>
#include <Eigen/Dense>

// Local includes.
#include "../Basic/VertexR3.hpp"
#include "../Basic/Face.hpp"
#include "BarycentricCoordinatesR3.hpp"

namespace gbc {

    // Mean value coordinates in R3.
    class MaximumLikelihoodCoordinatesR3 : public BarycentricCoordinatesR3 {

    public:
        // Constructor.
        MaximumLikelihoodCoordinatesR3(const std::vector<VertexR3> &v, const std::vector<Face> &f, const double tol = 1.0e-10) : super(v, f, tol) 
        {
            n_v = _v.size();
            n_f = _f.size();
        }



        // Return name of the coordinate function.
        inline std::string name() const {
            return "MaximumLikelihoodCoordinatesR3";
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR3 class.
        void compute(std::vector<VertexR3> &p);

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR3> &p) {
            compute(p);
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR3 super;

        void calculateMaximumLikelihoodCoordinates();

        void computeIntegralUnitNormals(const Eigen::MatrixXd &v,  const Eigen::MatrixXi &f, Eigen::MatrixXd &transMatrix, Eigen::MatrixXd &normalized_integral_outward_allfaces);

        Eigen::VectorXd f_gradient(Eigen::VectorXd x, Eigen::MatrixXd v);
        
        Eigen::MatrixXd f_Hessian(Eigen::VectorXd x, Eigen::MatrixXd v);

        double f(Eigen::VectorXd x, Eigen::MatrixXd v)
        {
            int n = v.cols();
            double z = 0;

            for (int i = 0; i < n; ++i)
            {
                z -= log(n + x.dot(v.col(i)));
            }

            return z;
        }


        double f(unsigned n, const double *x, double *grad, void *data)
        {
            Eigen::MatrixXd* v = (Eigen::MatrixXd *) data;
            Eigen::Map<Eigen::Vector3d> gradient(grad, n);
            Eigen::Map<const Eigen::Vector3d> X(x);
            gradient = - (*v) * (1 / ((((*v).transpose()*X).array() + 1.0*(*v).cols()))).matrix();

            return f(X, *v);
        }

        Eigen::MatrixXd cage_v; // 3 * _n
        Eigen::MatrixXi cage_f; // 3 * _f
        Eigen::MatrixXd model_v;
        Eigen::MatrixXd mlc;
    };

} // namespace gbc

#endif // GBC_MAXIMUMLIKELIHOOD_H
