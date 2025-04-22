/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-21 22:43:48
 * @LastEditTime: 2024-10-10 20:40:46
 * @Description: From the source https://github.com/DanStroeter/CageModeler
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/BarycentricCoordinates/LocalBarycentricCoordinatesR3.h
 */

#ifndef GBC_LOCALBARYCENTRICCOORDINATESR3_H
#define GBC_LOCALBARYCENTRICCOORDINATESR3_H

// includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../Basic/VertexR3.hpp"
#include "../Basic/Face.hpp"
#include "BarycentricCoordinatesR3.hpp"
#include "../ThirdPartyLib/LBC/include/LBC/LBCSolver.h"

namespace gbc {

    // Mean value coordinates in R3.
    class LocalBarycentricCoordinatesR3 : public BarycentricCoordinatesR3 {

    public:
        // Constructor.
        LocalBarycentricCoordinatesR3(const std::vector<VertexR3> &v, const std::vector<Face> &f, const double tol = 1.0e-10) : super(v, f, tol) 
        {
            n_v = _v.size();
            n_f = _f.size();
        }



        // Return name of the coordinate function.
        inline std::string name() const {
            return "LocalBarycentricCoordinatesR3";
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR3 class.
        void compute(std::vector<VertexR3> &p);

        void set_lbc_scheme(int lbc_scheme, int numBBWSteps = 3000)
        {
            _lbc_scheme = lbc_scheme;
            _numBBWSteps = numBBWSteps;
        }


        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR3> &p) {
            compute(p);
        }

        // About the initialisation of some data
        void Init(Eigen::MatrixXd & V, Eigen::MatrixXi & T, Eigen::MatrixXd & C, Eigen::MatrixXi & CF); 

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR3 super;   
        int _lbc_scheme;
        int _numBBWSteps;
        Eigen::MatrixXd sample_points; //tetrahedral sectioning vertex 3*n
        Eigen::MatrixXi cell_vertices; // Tetrahedra in tetrahedral sectioning 4 * n
        LBC::IndexVector control_point_idx;
        Eigen::MatrixXd _V, _V_model, _C;		// C is the set of cage vertices
	    Eigen::MatrixXi _T, _T_model, _CF;
        Eigen::MatrixXi BE; // BE by 2 list of bone edge indices into C
        Eigen::MatrixXi CE; // CE by 2 list of cage edge indices into *P*
        Eigen::VectorXi b; //  b list of boundary indices (indices into V of vertices which have known, fixed values)
        Eigen::MatrixXd bbc; // b by #weights list of known/fixed values for boundary vertices 
        Eigen::MatrixXd W, W_interpolated; //W_interpolated is GBC
    };

} // namespace gbc

#endif // GBC_LOCALBARYCENTRICCOORDINATESR3_H
