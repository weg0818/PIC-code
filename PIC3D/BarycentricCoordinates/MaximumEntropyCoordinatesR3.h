/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-21 15:40:17
 * @LastEditTime: 2024-10-10 20:42:30
 * @Description: From the source https://github.com/DanStroeter/CageModeler
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/BarycentricCoordinates/MaximumEntropyCoordinatesR3.h
 */

#ifndef GBC_MAXIMUMEntropy_H
#define GBC_MAXIMUMEntropy_H

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
    class MaximumEntropyCoordinatesR3 : public BarycentricCoordinatesR3 {

    public:
        // Constructor.
        MaximumEntropyCoordinatesR3(const std::vector<VertexR3> &v, const std::vector<Face> &f, const double tol = 1.0e-10) : super(v, f, tol) 
        {
            n_v = _v.size();
            n_f = _f.size();
        }



        // Return name of the coordinate function.
        inline std::string name() const {
            return "MaximumEntropyCoordinatesR3";
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR3 class.
        void compute(std::vector<VertexR3> &p);

        void set_mec_flag(int flag = 1)
        {
            mec_flag = flag;
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR3> &p) {
            compute(p);
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR3 super;

        void calculateMaximumEntropyCoordinates();
        void priorFunctions(const Eigen::Vector3d v, const Eigen::MatrixXd &cage_v, const Eigen::MatrixXi &cage_f, std::vector<std::vector<int>> & adjs, Eigen::VectorXd &priors, int mec_flag);
        double areaOfTriangle(const Eigen::Vector3d & a, const Eigen::Vector3d & b, const Eigen::Vector3d & c);

        Eigen::MatrixXd cage_v; // 3 * _n
        Eigen::MatrixXi cage_f; // 3 * _f
        Eigen::MatrixXd model_v;
        Eigen::MatrixXd mec;
        int mec_flag = 1;
    };

} // namespace gbc

#endif // GBC_MAXIMUMEntropy_H
