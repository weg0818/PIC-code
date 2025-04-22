/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-23 15:10:59
 * @LastEditTime: 2025-01-21 22:00:49
 * @Description: 
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/BarycentricCoordinates/HarmonicR3.cpp
 */
#include "HarmonicR3.h"
#include "../Basic/GBCmath.h"


#include <igl/boundary_conditions.h>



void gbc::HarmonicR3::compute(std::vector<VertexR3> &p)
{

    Basic::vecVertexToMatXd2(p, _V_model);
    if (!igl::harmonic(_V, _T, b, bbc, 1, W))
    {
        std::cerr << "Failed to compute harmonic weights!\n";
        return;
    }

    std::cout << W.rows() << ' ' << W.cols() << std::endl;

    Eigen::MatrixXd M;
    Basic::interpolateWeightsInEmbedding(_V_model, W, _V, _T, _C.rows(), W_interpolated, M);
    std::cout << "W" << ' ' << W.rows() << ' ' << W.cols() << std::endl;
	std::cout << "W_interpolated" << ' ' << W_interpolated.rows() << ' ' << W_interpolated.cols() << std::endl;
    if (W_interpolated.rows() == p.size())
    {
        for (int i = 0; i < p.size(); ++i)
        {
            std::vector<double> & vecDouble = p[i].b();
            vecDouble.clear();
            vecDouble.resize(n_v, 0.0);
            for (int j = 0; j < n_v; ++j)
            {
                vecDouble[j] = W_interpolated(i, j);
            }
        }
    }
    else
    {
        std::cout << "No calculation of hc" << std::endl;
    }
}

void gbc::HarmonicR3::Init(Eigen::MatrixXd & V, Eigen::MatrixXi & T, Eigen::MatrixXd & C, Eigen::MatrixXi & CF)
{

    _V = V;
    _T = T;
    _C = C;

    Eigen::VectorXi P;
    P.resize(C.rows());
	for (int i = 0; i < C.rows(); ++i) {
		P(i) = i;
	}

    if (!igl::boundary_conditions(V, T, C, P, BE, CE, CF, b, bbc))
	{
		std::cerr << "Failed to extract boundary conditions for cage!\n";
		return;
	}
    else
	{
		std::cout << "Cage Boundary Extraction Successful!" << std::endl;
	}
}