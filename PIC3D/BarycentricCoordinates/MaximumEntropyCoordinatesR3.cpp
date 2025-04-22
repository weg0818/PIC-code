#include "MaximumEntropyCoordinatesR3.h"
#include "../Basic/GBCmath.h"

#include <igl/vertex_triangle_adjacency.h>

void gbc::MaximumEntropyCoordinatesR3::compute(std::vector<VertexR3> &p)
{
    Basic::vecVertexToMatXd(_v, cage_v);
    Basic::vecFaceToMatXi(_f, cage_f);
    Basic::vecVertexToMatXd(p, model_v);
    calculateMaximumEntropyCoordinates();

    if (mec.cols() == p.size())
    {
        for (int i = 0; i < p.size(); ++i)
        {
            std::vector<double> & vecDouble = p[i].b();
            vecDouble.clear();
            vecDouble.resize(n_v, 0.0);
            for (int j = 0; j < n_v; ++j)
            {
                vecDouble[j] = mec(j, i);
            }
        }
    }
    else
    {
        std::cout << "No calculation of mec" << std::endl;
    }

}

void gbc::MaximumEntropyCoordinatesR3::calculateMaximumEntropyCoordinates()
{
    int nv = n_v;
    int nf = n_f;

    Eigen::MatrixXd shifted_cage_v = cage_v;

    mec.resize(nv, model_v.cols());

    std::vector<std::vector<int>> allAdjFaces;
    std::vector< std::vector< int > > VF;
    igl::vertex_triangle_adjacency(nv, cage_f, VF, allAdjFaces);
    
    for (int ii = 0; ii < model_v.cols(); ++ii)
    {
         // compute the prior function
        Eigen::VectorXd priors = Eigen::VectorXd::Ones(nv);
        priorFunctions(model_v.col(ii), cage_v, cage_f, allAdjFaces, priors, mec_flag);

        
        shifted_cage_v = cage_v.colwise() - model_v.col(ii);
        // damp Newton
        double error = 1e-10;
        double rho = .7;
        double sigma = .3;

        Eigen::Vector3d lambda(0.0, 0.0, 0.0);
        while (true)
        {
            double step_length = 1.0;
            Eigen::VectorXd Zi = Eigen::VectorXd::Zero(nv);
            for (int i = 0; i < nv; ++i)
            {
                Zi[i] = priors[i] * exp(-lambda.dot(shifted_cage_v.col(i)));
            }
            // compute the gradient
            Eigen::Vector3d gradient(0.0, 0.0, 0.0);

            gradient = -shifted_cage_v*Zi;
            double Z = Zi.sum();
            gradient /= Z;

            if (std::isnan(gradient.norm())) 
            {
                std::cerr<< "There is numerical instability in the process of solving optimization problems. (please check the " << ii << "-th point)" << std::endl;
                assert(false);
            }

            if (gradient.norm() < error)
            {
                mec.col(ii) = Zi / Z;
                if ((cage_v * mec.col(ii) - model_v.col(ii)).norm()>1e-7)
                    std::cerr<< "---------------- The maximum entropy coordinates of the " << ii << "-th point may be wrong, as c*lambda-v=" << (cage_v * mec.col(ii) - model_v.col(ii)).norm() << " > 1e-7" << std::endl;
                break;
            }

            // compute the Hessian matrix
            Eigen::Vector3d SUMZivi(0.0, 0.0, 0.0);
            Eigen::Matrix3d Hessian = Eigen::Matrix3d::Zero();
            for (int i = 0; i < nv; ++i)
            {
                Eigen::Vector3d Zivi = Zi[i] * shifted_cage_v.col(i);
                SUMZivi += Zivi;
                Hessian += Zivi * (shifted_cage_v.col(i).transpose()); // 3-by-3
            }
            Hessian *= Z;
            Hessian -= (SUMZivi * SUMZivi.transpose());
            Hessian /= (Z * Z);

            // determine Newton search direction
            SUMZivi = -Hessian.colPivHouseholderQr().solve(gradient);

            if (gradient.norm() > 1e-4) {
                // determine step size by using Armijo linear search
                while (step_length > error)
                {
                    Eigen::Vector3d lambda_temp = lambda + step_length * SUMZivi;
                    for (int i = 0; i < nv; ++i)
                    {
                        Zi[i] = priors[i] * exp(-lambda_temp.dot(shifted_cage_v.col(i)));
                    }
                    if (log(Zi.sum()) < log(Z) + sigma * step_length * gradient.dot(SUMZivi))
                    {
                        break;
                    }
                    step_length *= rho;
                }
            }
            lambda += (step_length*SUMZivi);
        }
    }
}

void gbc::MaximumEntropyCoordinatesR3::priorFunctions
    (const Eigen::Vector3d v, const Eigen::MatrixXd &cage_v, const Eigen::MatrixXi &cage_f, std::vector<std::vector<int>> & adjs, Eigen::VectorXd &priors, int mec_flag)
{
    for (int i = 0; i < cage_v.cols(); ++i)
    {

        if (mec_flag == 1) // MEC-1
        {
            for (std::vector<int>::iterator iter = adjs.at(i).begin(); iter != adjs.at(i).end(); ++iter)
            {
                priors[i] *= areaOfTriangle(v, cage_v.col(cage_f(0, *iter)), cage_v.col(cage_f(1, *iter))) + areaOfTriangle(v, cage_v.col(cage_f(1, *iter)), cage_v.col(cage_f(2, *iter))) +
                             areaOfTriangle(v, cage_v.col(cage_f(2, *iter)), cage_v.col(cage_f(0, *iter))) - areaOfTriangle(cage_v.col(cage_f(0, *iter)), cage_v.col(cage_f(1, *iter)), cage_v.col(cage_f(2, *iter)));
            }
        }
        if (mec_flag == 2) // MEC-2
        {
            for (std::vector<int>::iterator iter = adjs.at(i).begin(); iter != adjs.at(i).end(); ++iter)
            {
                priors[i] *= areaOfTriangle(v, cage_v.col(cage_f(0, *iter)), cage_v.col(cage_f(1, *iter))) * areaOfTriangle(v, cage_v.col(cage_f(1, *iter)), cage_v.col(cage_f(2, *iter))) *
                                 areaOfTriangle(v, cage_v.col(cage_f(2, *iter)), cage_v.col(cage_f(0, *iter))) +
                             (v - cage_v.col(cage_f(0, *iter))).dot(v - cage_v.col(cage_f(1, *iter))) + (v - cage_v.col(cage_f(1, *iter))).dot(v - cage_v.col(cage_f(2, *iter))) + (v - cage_v.col(cage_f(2, *iter))).dot(v - cage_v.col(cage_f(0, *iter)));
            }
        }
        priors[i] = 1.0 / priors[i];
    }
    priors /= priors.sum();
}

double gbc::MaximumEntropyCoordinatesR3::areaOfTriangle(const Eigen::Vector3d & a, const Eigen::Vector3d & b, const Eigen::Vector3d & c)
{
    return 0.5 * (b - a).cross(c - a).norm();
}