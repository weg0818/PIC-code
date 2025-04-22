/*** 
 * @Author: Mao jiaping
 * @Date: 2024-09-20 09:33:38
 * @LastEditTime: 2024-10-10 17:06:41
 * @Description: writed by Mao jiaping
 * @FileName: 
 * @FilePath: /Codes/C++/gbc-master/2d/coords/PointwiseIterativeR2_Scaled.h
 */

#ifndef SPOINTWISEITERATIVER2_H
#define SPOINTWISEITERATIVER2_H

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/MeshR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

#include "../extra/libs/eigen/Eigen/Dense"
#include "../extra/libs/libigl/include/igl/exact_geodesic.h"

namespace gbc {

    // ScaledPointwiseIterativeR2 in R2.
    class ScaledPointwiseIterativeR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        ScaledPointwiseIterativeR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) 
        {
            n_v = _v.size();
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "ScaledPointwiseIterativeR2";
        }
        
        void compute(const VertexR2 &p, std::vector<double> &b, std::vector<double> & cd1_i, std::vector<double> & cd2_i);

        // Compute the coordinates at p using the internal storage from the VertexR2 class.
        inline void compute(VertexR2 &p){
            compute(p, p.b(), _cd1_i, _cd2_i);
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb)
        {
            const int numP = p.size();
            
            bb.resize(numP);
            for (int i = 0; i < numP; ++i) 
            {   
                compute(p[i], bb[i], _cd1[i], _cd2[i]);
            }
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p){

            const int numP = p.size();
#pragma omp parallel for    
            for (int i = 0; i < numP; ++i) 
            {
                compute(p[i], p[i].b(), _cd1[i], _cd2[i]);
            }
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

        void set_scaled_cd(std::vector<std::vector<double>> & cd1, std::vector<std::vector<double>> & cd2)
        {
            _cd1 = cd1;
            _cd2 = cd2;
        }

        void set_f(std::vector<Face> & f)
        {
            _f = f;
        }

        void compute_scald_cd(std::vector<std::vector<double>> & cd1, std::vector<std::vector<double>> & cd2, std::vector<VertexR2> & p);
    
    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        void bcInterior(const VertexR2 &p, std::vector<double> &b, std::vector<double> & cd1_i, std::vector<double> & cd2_i);
    
        void update_delta(std::vector<VertexR2>& p1, std::vector<double>& p1_coe, std::vector<VertexR2>& p2, std::vector<std::vector<double>>& p2_coe,
            VertexR2& delta, std::vector<double>& delta_coe, std::vector<double> & cd1_i, std::vector<double> & cd2_i);
    
        void get_rhoi_information(std::vector<VertexR2>& pi, double alpha,
            VertexR2& delta, std::vector<double>& delta_coe,  
            std::vector<double>& rhoi_coe, std::vector<double>& d_coei);

        void cal_p_bisectorPos(const Eigen::MatrixXd& p, const Eigen::MatrixXd& x, const std::vector<int>& vertex_markers, 
            std::vector<double> & pbP, std::vector<int> & boundryPoints);

        void cal_p_bisectorDis(std::vector<int>& boundryPoints, const std::unordered_map<int, Eigen::VectorXd>& boundryDistance, 
            const std::vector<double>& pbP, int index, Eigen::VectorXd & cd2);


        inline double f1(double x)
        {
            return x * x * x + 3 * x * x * (1 - x);
        };

        std::vector<std::vector<double>> _cd1;
        std::vector<std::vector<double>> _cd2;
        std::vector<double> _cd1_i;
        std::vector<double> _cd2_i;
        std::vector<Face> _f;

    };

} // namespace gbc

#endif // GBC_SPOINTWISEITERATIVER2_HPP
