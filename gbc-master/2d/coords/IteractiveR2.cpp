#include "../extra/math.h"
#include "IteractiveR2.h"
#include "MeanValueR2.hpp"

using namespace gbc;

void gbc::IteractiveR2::compute(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior(p, b);
    return;
}

void gbc::IteractiveR2::compute_iteractive_times(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior_i(p);
    return;
}

void gbc::IteractiveR2::bcInterior_i(const VertexR2 &p)
{
    std::vector<VertexR2> p1;

    // translation
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(2, n_v);
    for (int i = 0; i < n_v; ++i)
    {
        P(0, i) = _v[i][0] - p[0];
        P(1, i) = _v[i][1] - p[1];
        p1.push_back(_v[i] - p);
    }

    Eigen::MatrixXd M_0 = Eigen::MatrixXd::Zero(n_v, n_v);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n_v, n_v);
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = 1 / p1[i].length();
        M_0(i, i) = tmp;
        M(i, i) = tmp;
    }

    Eigen::MatrixXd p_i = P * M_0; // 2 * n
    int times = 0;
    while (1)
    {
        times++;
        Eigen::MatrixXd M_k = Eigen::MatrixXd::Zero(n_v, n_v);
        get_M_k(p_i, M_k);
        p_i = p_i * M_k;
        M = M * M_k;

        std::vector<double> bc;
        std::vector<VertexR2> poly;
        MatrixXdToVecVertex(p_i, poly);
        MeanValueR2 mvc = MeanValueR2(poly);
        mvc.compute(VertexR2(0.0, 0.0), bc);
        bool flag = false;
        for (int i = 0; i < n_v; ++i)
        {
            if (bc[i] < 0)
            {
                flag = true;
                break;
            }
        }
        if (flag == false || times > 1000)
        {
            break;
        }
    }
    _k = std::max(_k, times);

}

void gbc::IteractiveR2::MatrixXdToVecVertex(Eigen::MatrixXd & p_i, std::vector<VertexR2> & poly)
{
    for (int i = 0; i < n_v; ++i)
    {
        poly.push_back(VertexR2(p_i(0,i), p_i(1,i)));
    }
}

void gbc::IteractiveR2::get_M_k(Eigen::MatrixXd & p_i, Eigen::MatrixXd & M_k)
{
    for (int i = 0; i < n_v; ++i)
    {
        int ip = (i - 1 + n_v) % n_v;
        double tmp = 1 / (VertexR2(p_i(0, i), p_i(1, i)) + VertexR2(p_i(0, ip), p_i(1, ip))).length();
        M_k(ip, i) = tmp;
        M_k(i, i) = tmp;
    }
}


void gbc::IteractiveR2::bcInterior(const VertexR2 &p, std::vector<double> &b)
{
    std::vector<VertexR2> p1;

    // translation
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(2, n_v);
    for (int i = 0; i < n_v; ++i)
    {
        P(0, i) = _v[i][0] - p[0];
        P(1, i) = _v[i][1] - p[1];
        p1.push_back(_v[i] - p);
    }

    Eigen::MatrixXd M_0 = Eigen::MatrixXd::Zero(n_v, n_v);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n_v, n_v);
    for (int i = 0; i < n_v; ++i)
    {
        if (p1[i].length() == 0)
        {
            return;
        }
        double tmp = 1 / p1[i].length();
        M_0(i, i) = tmp;
        M(i, i) = tmp;
    }

    Eigen::MatrixXd p_i = P * M_0; // 2 * n
    int times = _k;
    while (times--)
    {
        Eigen::MatrixXd M_k = Eigen::MatrixXd::Zero(n_v, n_v);
        get_M_k(p_i, M_k);
        p_i = p_i * M_k;
        M = M * M_k;
    }
    std::vector<double> bc;
    std::vector<VertexR2> poly;
    MatrixXdToVecVertex(p_i, poly);
    MeanValueR2 mvc = MeanValueR2(poly);
    mvc.compute(VertexR2(0.0, 0.0), bc);

    Eigen::MatrixXd mat_bc = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(bc.data(), n_v, 1);
    mat_bc = M * mat_bc;
    for (int i = 0; i < n_v; ++i)
    {
        b[i] = mat_bc(i, 0);
    }

    double b_sum = gbc::sum(b);
    gbc::div(b, b_sum);
}
