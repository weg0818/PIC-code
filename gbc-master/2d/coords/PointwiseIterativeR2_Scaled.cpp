#include "PointwiseIterativeR2_Scaled.h"
#include <omp.h>
#include <vector>
#include "../extra/math.h"

using namespace gbc;

void gbc::ScaledPointwiseIterativeR2::compute(const VertexR2 &p, std::vector<double> &b, std::vector<double> & cd1_i, std::vector<double> & cd2_i)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior(p, b, cd1_i, cd2_i);
    return;
}


void gbc::ScaledPointwiseIterativeR2::bcInterior(const VertexR2 &p, std::vector<double> &b, std::vector<double> & cd1_i, std::vector<double> & cd2_i)
{
    std::vector<VertexR2> p1;
    for (int i = 0; i < n_v; ++i)
    {
        p1.push_back(_v[i]);
    }

    // translation
    for (auto &p1_i : p1)
    {
        p1_i -= p;
    }

    // calutate \hat{P}
    std::vector<double> p1_coe(n_v, 0.0);
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        p1[i] *= 1 / tmp;
        p1_coe[i] += 1 / tmp;
    }

    // calutate \tilde{P}
    std::vector<VertexR2> p2(n_v);
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0));
    int ip;
    for (int i = 0; i < n_v; ++i)
    {   
        ip = (i + 1) % n_v;
        p2[i] = VertexR2(p1[i][0] + p1[ip][0], p1[i][1] + p1[ip][1]);
        tmp = p2[i].length();
        p2[i] *= 1 / tmp;
        p2_coe[i][0] += p1_coe[i] / tmp;
        p2_coe[i][1] += p1_coe[ip] / tmp;
    }

    // iteration
    for (int i = 0; i < n_v; ++i)
    {
        // VertexR2 delta = p1[i];
        // std::vector<double> delta_coe(n_v, 0.0);
        // delta_coe[i] = p1_coe[i];
    
        VertexR2 delta = p2[i]; // Initialise \delta
        std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
        delta_coe[i] = p2_coe[i][0];
        delta_coe[(i+1)%n_v] = p2_coe[i][1];
        
        update_delta(p1, p1_coe, p2, p2_coe, delta, delta_coe, cd1_i, cd2_i);
        double sum_delta = gbc::sum(delta_coe);
        gbc::div(delta_coe, sum_delta);
        gbc::add(b, delta_coe);    
    }
    double sum_b = gbc::sum(b);
    gbc::div(b, sum_b);
}

void gbc::ScaledPointwiseIterativeR2::update_delta(std::vector<VertexR2>& p1, std::vector<double>& p1_coe, 
            std::vector<VertexR2>& p2,std::vector<std::vector<double>> & p2_coe,
            VertexR2& delta, std::vector<double>& delta_coe, std::vector<double> & cd1_i, std::vector<double> & cd2_i)
{
    double delta_norm = delta.length();
    if (delta_norm > 1.0)
    {
        delta_norm = 1.0;
    }
    int times = 0;
    while (delta_norm > _tol)
    {
        times += 1;
        double alpha = std::asin(delta_norm);
        // alpha = 0;
        std::vector<double> rho1_coe(n_v, 0.0);
        std::vector<double> rho2_coe(n_v, 0.0);
        std::vector<double> d1_coe(n_v, 0.0);
        std::vector<double> d2_coe(n_v, 0.0);
        std::vector<double> new_delta_coe(n_v, 0.0);

        get_rhoi_information(p1, alpha, delta, delta_coe, rho1_coe, d1_coe);
        get_rhoi_information(p2, alpha, delta, delta_coe, rho2_coe, d2_coe);
        
        double sum = 0.0;
        double sum_delta_rho = 0.0; 
        for (int i = 0; i < n_v; ++i)
        {
            sum += rho1_coe[i] / cd1_i[i] + rho2_coe[i] / cd2_i[i];
            sum_delta_rho += rho1_coe[i] * d1_coe[i] / cd1_i[i] + rho2_coe[i] * d2_coe[i] / cd2_i[i];
        }

        double ip;
        std::vector<VertexR2> d1(n_v);
        std::vector<VertexR2> d2(n_v);
        for (int i = 0; i < n_v; ++i)
        {
            ip = (i - 1 + n_v) % n_v;
            new_delta_coe[i] = sum_delta_rho * delta_coe[i] + rho1_coe[i] * (1 - d1_coe[i]) * p1_coe[i] / cd1_i[i]
                            + rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][0] / cd2_i[i]
                            + rho2_coe[ip] * (1 - d2_coe[ip]) * p2_coe[ip][1] / cd2_i[ip];
            d1[i] = (rho1_coe[i] / cd1_i[i]) * (d1_coe[i] * delta + (1 - d1_coe[i]) * p1[i]);
            d2[i] = (rho2_coe[i] / cd2_i[i]) * (d2_coe[i] * delta + (1 - d2_coe[i]) * p2[i]);
        }

        delta[0] = 0.0;
        delta[1] = 0.0;
        for (int i = 0; i < n_v; ++i)
        {
            delta_coe[i] = new_delta_coe[i] / sum;
            delta += d1[i] + d2[i];
        }
        delta = delta / sum;
        delta_norm = delta.length();

        if (times > 100)
        {
            int a =1;
        }
   }
}

void gbc::ScaledPointwiseIterativeR2::get_rhoi_information(std::vector<VertexR2>& pi, double alpha,
        VertexR2& delta, std::vector<double>& delta_coe,
        std::vector<double>& rhoi_coe, std::vector<double>& di_coe)
{
    double tmp, theta;
    double delta_norm = delta.length();
    for (int i = 0; i < pi.size(); ++i)
    {
        tmp = (pi[i][0] * delta[0] + pi[i][1] * delta[1]) / (pi[i].length() * delta_norm);
        if (tmp > 1.0) 
        {
            tmp = 1.0;
        }
        if (tmp < -1.0)
        {
            tmp = -1.0;
        }
        theta = std::acos(tmp);
        if (theta > PI / 2 - alpha)
        {
            rhoi_coe[i] = f1((theta - PI/2 + alpha) / (PI/2 + alpha));
            VertexR2 del_p = delta - pi[i];
            di_coe[i] = 0.0 - gbc::dot(pi[i], del_p) / gbc::dot(del_p, del_p);
        }
    }
}

void gbc::ScaledPointwiseIterativeR2::compute_scald_cd(std::vector<std::vector<double>> & cd1_in, std::vector<std::vector<double>> & cd2_in, std::vector<VertexR2> & p)
{
    omp_set_num_threads(4);
    int n = n_v;
    int n_x = p.size();
 
    Eigen::MatrixXd polygon;  // n*2
    vecVertexToMatXd2(_v, polygon);
    Eigen::MatrixXd x;
    vecVertexToMatXd2(p, x); // n*2

    Eigen::MatrixXd cds1(n, n_x);
    Eigen::MatrixXd cds2 = Eigen::MatrixXd::Zero(n_x, n);
    std::unordered_map<int, Eigen::VectorXd> boundryDistance;

    Eigen::MatrixXd v(n_x, 3);
    v << x, Eigen::MatrixXd::Zero(n_x, 1); //  Extend x to 3 dimensions
    Eigen::MatrixXi f; // n * 3
    vecFaceToMatXi2(_f, f);

    Eigen::VectorXi vt(n_x);
    for (int i = 0; i < n_x; i++) vt(i) = i;

    std::cout << "Compute the distances of vertices" << std::endl;
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        Eigen::VectorXi vs(1);
        vs(0) = i;
        std::cout << i << std::endl;
        // igl::exact_geodesic function call
        Eigen::VectorXd D;
        igl::exact_geodesic(v, f, vs, Eigen::VectorXi(), vt, Eigen::VectorXi(), D);
        
        cds1.row(i) = D.transpose();
        boundryDistance[i] = D;
    }

    cd1_in.resize(n_x, std::vector<double>(n, 0.0));
    Eigen::MatrixXd cds1_T = cds1.transpose();
    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cd1_in[i][j] = cds1_T(i,j);
        }
    }


    // Is it on the border
    std::vector<int> vertex_markers(n_x);
#pragma omp parallel for   
    for (int i = 0; i < n_x; ++i)
    {
        if (computeBoundaryCoordinates(p[i], p[i].b()))
        {
            vertex_markers[i] = 1; // at the border
        }
        else
        {
            vertex_markers[i] = 0;
        }
    }

    std::cout << "Calculate the  distances of the boundaries" << std::endl;
#pragma omp parallel for 
    for (int i = n; i < n_x; i++) {
        
        // not at the border
        if (vertex_markers[i] == 0) {
            std::cout << i << std::endl;
            std::vector<double> pbP; // Ratio of points about the two nearest points on a side
            std::vector<int> boundryPoints; // Index of the two nearest points on the edge
            cal_p_bisectorPos(polygon.rowwise() - x.row(i), x.rowwise() - x.row(i), vertex_markers, pbP, boundryPoints);

            for (int bp : boundryPoints) {
                if (boundryDistance.find(bp) == boundryDistance.end()) {
                    Eigen::VectorXi vs(1);
                    vs(0) = bp;
                    
                    Eigen::VectorXd D;
                    igl::exact_geodesic(v, f, vs, Eigen::VectorXi(), vt, Eigen::VectorXi(), D);

                    boundryDistance[bp] = D;
                }
            }
            Eigen::VectorXd cd2(boundryPoints.size() / 2);
            cal_p_bisectorDis(boundryPoints, boundryDistance, pbP, i, cd2);
            cds2.row(i) = cd2;
        }
    }

    cd2_in.resize(n_x, std::vector<double>(n, 0.0));
    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cd2_in[i][j] = cds2(i, j);
        }
    }
}

// Calculate the position of the angle bisector
void gbc::ScaledPointwiseIterativeR2::cal_p_bisectorPos(const Eigen::MatrixXd& p, const Eigen::MatrixXd& x, const std::vector<int>& vertex_markers, 
        std::vector<double> & pbP, std::vector<int> & vt) 
{
    Eigen::VectorXd dis = p.rowwise().norm(); // Calculated Distance
    int n = p.rows();
    Eigen::MatrixXd p_bisector(n, 2);

    for (int i = 0; i < n; i++) {
        int ip = (i + 1) % n;
        double c = dis(i) + dis(ip);
        p_bisector.row(i) = (dis(ip) * p.row(i) + dis(i) * p.row(ip)) / c;
    }


    for (int i = 0; i < n; i++) {
        int l_index = i;
        double l_dis = (p_bisector.row(i) - p.row(i)).norm();
        int r_index = (i + 1) % n;
        double r_dis = (p_bisector.row(i) - p.row(r_index)).norm();

        for (int j = 0; j < vertex_markers.size(); j++) {
            if (vertex_markers[j] == 1) { // If on the border
                double dir2 = (p_bisector.row(i) - p.row(i)).dot(p_bisector.row(i) - x.row(j));
                double dir1 = cross(x(j, 0) - p(i, 0), x(j, 1) - p(i, 1), x(j, 0) - p((i + 1) % n, 0), x(j, 1) - p((i + 1) % n, 1));

                if (abs(dir1) < 1e-10 && dir2 > 0) {
                    double dis = (p_bisector.row(i) - x.row(j)).norm();
                    if (dis < l_dis) {
                        l_index = j;
                        l_dis = dis;
                    }
                } else if (abs(dir1) < 1e-10 && dir2 < 0) {
                    double dis = (p_bisector.row(i) - x.row(j)).norm();
                    if (dis < r_dis) {
                        r_index = j;
                        r_dis = dis;
                    }
                }
            }
        }

        // Calculate ratios and store
        pbP.push_back(r_dis / (r_dis + l_dis));
        vt.push_back(l_index);
        vt.push_back(r_index);
    }
}

// Calculate the angle bisector distance
void ScaledPointwiseIterativeR2::cal_p_bisectorDis(std::vector<int>& boundryPoints, const std::unordered_map<int, Eigen::VectorXd>& boundryDistance, const std::vector<double>& pbP, int index, Eigen::VectorXd & cd2) {
    int n = boundryPoints.size() / 2;

    for (int i = 0; i < n; i++) {
        double distance1 = boundryDistance.at(boundryPoints[2 * i])(index);
        double distance2 = boundryDistance.at(boundryPoints[2 * i + 1])(index);
        cd2(i) = distance1 * pbP[i] + distance2 * (1 - pbP[i]);
    }
}
