/*** 
 * @Author: Mao jiaping
 * @Date: 
 * @LastEditTime: 2024-10-10 20:38:45
 * @Description: From the source https://github.com/DanStroeter/CageModeler
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/BarycentricCoordinates/HarmonicR3.h
 */


#ifndef GBC_HARMONICR3_H 
#define GBC_HARMONICR3_H

// includes.
#include <vector>
#include <cassert>
#include <cmath>
#include <igl/harmonic.h>

// Local includes.
#include "../Basic/VertexR3.hpp"
#include "../Basic/Face.hpp"
#include "BarycentricCoordinatesR3.hpp"

namespace gbc {

    // Mean value coordinates in R3.
    class HarmonicR3 : public BarycentricCoordinatesR3 {

    public:
        // Constructor.
        HarmonicR3(const std::vector<VertexR3> &v, const std::vector<Face> &f, const double tol = 1.0e-10) : super(v, f, tol) 
        {
            n_v = _v.size();
            n_f = _f.size();
        }



        // Return name of the coordinate function.
        inline std::string name() const {
            return "HarmonicR3";
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR3 class.
        void compute(std::vector<VertexR3> &p);


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
        Eigen::MatrixXd _V, _V_model, _C;		//  C is the set of cage vertices
	    Eigen::MatrixXi _T, _T_model, _CF;
        Eigen::MatrixXi BE; // BE by 2 list of bone edge indices into C
        Eigen::MatrixXi CE; // CE by 2 list of cage edge indices into *P*
        Eigen::VectorXi b; //  b list of boundary indices (indices into V of vertices which have known, fixed values)
        Eigen::MatrixXd bbc; // #b by #weights list of known/fixed values for boundary vertices 
        Eigen::MatrixXd W, W_interpolated; // W_interpolated is GBC
    };

} // namespace gbc

#endif // GBC_HARMONICR3_H
