#include "LocalBarycentricCoordinatesR3.h"
#include "../Basic/GBCmath.h"

#include <igl/boundary_conditions.h>


void gbc::LocalBarycentricCoordinatesR3::compute(std::vector<VertexR3> &p)
{

	Basic::vecVertexToMatXd2(p, _V_model);
	// std::cout << _V_model.rows() << " " << _V_model.cols() << std::endl;

    LBC::DataSetup::WeightingScheme scheme = static_cast<LBC::DataSetup::WeightingScheme>(_lbc_scheme);

    if (bbc.rows() == 0)
    {
        return;
    }
    std::vector< LBC::DataSetup::CageBoundaryFacetInfo > boundary_facet_info;
    for (int i = 0; i < n_f; ++i)
	{
		const Eigen::Vector3i tri_indices(_f[i].v[0], _f[i].v[1], _f[i].v[2]);
		std::vector<int> boundary_points;

		for (int j = control_point_idx.size(); j < bbc.rows(); ++j)
		{
			auto const row = bbc.row(j);
			bool contains = true;
			for (int l = 0; l < 3; ++l)
			{
				if (row(tri_indices(l)) == 0)
				{
					contains = false;
					break;
				}
			}
			if (contains)
			{
				boundary_points.push_back(b(j));
			}
		}

		auto boundary_points_vec = LBC::IndexVector(boundary_points.size());
		for (int i = 0; i < boundary_points.size(); ++i)
		{
			boundary_points_vec(i) = boundary_points[i];
		}

		LBC::DataSetup::CageBoundaryFacetInfo info(tri_indices, boundary_points_vec);
		boundary_facet_info.push_back(info);
	}

    LBC::DataSetup ds(sample_points, control_point_idx, cell_vertices, boundary_facet_info, scheme);

	LBC::Param param;
	param.max_iterations = _numBBWSteps;
	param.relaxation_alpha = 1.65;
	param.convergence_check_frequency = 10;
	param.output_frequency_ratio = 10;
	param.rel_primal_eps = 0.00000000001;
	param.rel_dual_eps = 0.00000000001;
	param.penalty_weight = 10;

	std::cout << "Call the solver" << std::endl;
	LBC::LBCSolver solver(param, ds);
	solver.solve();
	W = ds.get_full_coordinate_values(solver.get_coordinates());

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
        std::cout << "No calculation of lbc" << std::endl;
    }
}

// V：tetrahedral sectioning vertex n*3
// T: Tetrahedra in tetrahedral sectioning n*4
void gbc::LocalBarycentricCoordinatesR3::Init(Eigen::MatrixXd & V, Eigen::MatrixXi & T, Eigen::MatrixXd & C, Eigen::MatrixXi & CF)
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
		std::cout << "Cage Boundary Extraction Successful" << std::endl;
	}

    int P_size = P.size();
    control_point_idx.resize(P.size());
	for (int i = 0; i < P.size(); ++i)
	{
		control_point_idx(i) = i;
	}

    sample_points.resize(V.cols(), V.rows());
    for (int i = 0; i < V.rows(); ++i)
	{
		sample_points.col(i) = V.row(i);
	}

    cell_vertices.resize(T.cols(), T.rows());
	for (int i = 0; i < T.rows(); ++i)
	{
		cell_vertices.col(i) = T.row(i);
	}
}


