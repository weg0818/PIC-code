#include <vector>
#include "PointwiseIterativeR2.h"
#include "../extra/math.h"

using namespace gbc;

void gbc::PointwiseIterativeR2::compute(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior(p, b);
    return;
}


void gbc::PointwiseIterativeR2::bcInterior(const VertexR2 &p, std::vector<double> &b)
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
    std::vector<double> p1_coe(n_v, 0.0); // Coefficients of the first projected point: 1/r_i
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        p1[i] *= 1 / tmp; // the first projected point
        p1_coe[i] += 1 / tmp;
    }

    // calutate \tilde{P}
    std::vector<VertexR2> p2(n_v); // the second projected point. 
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0)); // Coefficients of the second projected point
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
        VertexR2 delta = p1[i]; // Initialise \delta
        std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
        delta_coe[i] = p1_coe[i];

        update_delta(p1, p1_coe, p2, p2_coe, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
        double sum_delta = gbc::sum(delta_coe);
        gbc::div(delta_coe, sum_delta); // normalize
        gbc::add(b, delta_coe); 
    }
    double sum_b = gbc::sum(b);
    gbc::div(b, sum_b);
}

void gbc::PointwiseIterativeR2::update_delta(std::vector<VertexR2>& p1, std::vector<double>& p1_coe, 
            std::vector<VertexR2>& p2,std::vector<std::vector<double>> & p2_coe,
            VertexR2& delta, std::vector<double>& delta_coe)
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
        std::vector<double> rho1_coe(n_v, 0.0);
        std::vector<double> rho2_coe(n_v, 0.0);
        std::vector<double> d1_coe(n_v, 0.0);
        std::vector<double> d2_coe(n_v, 0.0);
        std::vector<double> new_delta_coe(n_v, 0.0);

        get_rhoi_information(p1, alpha, delta, delta_coe, rho1_coe, d1_coe); // get the \rho and u
        get_rhoi_information(p2, alpha, delta, delta_coe, rho2_coe, d2_coe); // get the \rho and u
        
        double sum = 0.0; // \sum{\rho}
        double sum_delta_rho = 0.0; // // \sum{\rho u}
        for (int i = 0; i < n_v; ++i)
        {
            sum += rho1_coe[i] + rho2_coe[i];
            sum_delta_rho += rho1_coe[i] * d1_coe[i] + rho2_coe[i] * d2_coe[i];
        }

        double ip;
        std::vector<VertexR2> d1(n_v);
        std::vector<VertexR2> d2(n_v);
        for (int i = 0; i < n_v; ++i)
        {
            ip = (i - 1 + n_v) % n_v;
            new_delta_coe[i] = sum_delta_rho * delta_coe[i] + rho1_coe[i] * (1 - d1_coe[i]) * p1_coe[i] 
                            + rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][0]
                            + rho2_coe[ip] * (1 - d2_coe[ip]) * p2_coe[ip][1];
            d1[i] = rho1_coe[i] * (d1_coe[i] * delta + (1 - d1_coe[i]) * p1[i]);
            d2[i] = rho2_coe[i] * (d2_coe[i] * delta + (1 - d2_coe[i]) * p2[i]);
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
   }
}

void gbc::PointwiseIterativeR2::get_rhoi_information(std::vector<VertexR2>& pi, double alpha,
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




void gbc::PointwiseIterativeR2::compute_iso_pic(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v + _iso_p_n, 0.0);
    // Interior.
    bcInterior_iso(p, b);
    return;
}


void gbc::PointwiseIterativeR2::bcInterior_iso(const VertexR2 &p, std::vector<double> &b)
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

    std::vector<VertexR2> q1 = _iso_p;
    for (auto & item : q1)
    {
        item -= p;
    }
    std::vector<VertexR2> q0 = q1;

    // calutate \hat{P}
    std::vector<double> p1_coe(n_v, 0.0); // Coefficients of the first projected point: 1/r_i
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        tmp = p1[i].length();
        if (tmp < 1e-10)
        {
            b[i] = 1.0;
            return;
        }
        p1[i] *= 1 / tmp; // the first projected point
        p1_coe[i] += 1 / tmp;
    }

    std::vector<double> q1_coe(q1.size(), 0.0);
    for (int i = 0; i < q1.size(); ++i)
    {
        tmp = q1[i].length();
        if (tmp < 1e-10)
        {
            b[n_v + i] = 1.0;
            return;
        }
        q1[i] *= 1 / tmp;
        q1_coe[i] += 1 / tmp;
    }

    // calutate \tilde{P}
    std::vector<VertexR2> p2(n_v); // the second projected point. 
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0)); // Coefficients of the second projected point
    std::vector<std::vector<int>> p2_index(n_v, std::vector<int>(2, 0));
    int ip;
    for (int i = 0; i < n_v; ++i)
    {   
        ip = (i + 1) % n_v;
        p2[i] = VertexR2(p1[i][0] + p1[ip][0], p1[i][1] + p1[ip][1]);
        tmp = p2[i].length();
        if (tmp < 1e-10)
        {
            double a1 = (_v[ip] - p).length();
            double b1 = (_v[i] - p).length();
            b[i] = a1 / (a1 + b1);
            b[ip] = b1 / (a1 + b1);
            return;
        }
        p2[i] *= 1 / tmp;
        p2_coe[i][0] += p1_coe[i] / tmp;
        p2_coe[i][1] += p1_coe[ip] / tmp;
        p2_index[i][0] = i;
        p2_index[i][1] = ip;
    }

    std::vector<VertexR2> q2;
    std::vector<std::vector<double>> q2_coe;
    int pos_q = n_v;
    std::vector<double> qq2_ceo(2);
    std::vector<int> qq2_index(2);
    for (int i = 0; i < _iso_p_line.size(); ++i)
    {   
        int index1 = _iso_p_line[i][0];
        int index2 = _iso_p_line[i][1];
        VertexR2 q_i = VertexR2(q1[index1][0] + q1[index2][0], q1[index1][1] + q1[index2][1]);
        tmp = q_i.length();
        if (tmp < 1e-10)
        {
            double a1 = (_iso_p[index2] - p).length();
            double b1 = (_iso_p[index1] - p).length();
            b[pos_q] = a1 / (a1 + b1);
            b[pos_q + 1] = b1 / (a1 + b1);
            return;
        }
        q_i *= 1 / tmp;
        q2.push_back(q_i);
        qq2_ceo[0] = q1_coe[index1] / tmp;
        qq2_ceo[1] = q1_coe[index2] / tmp;
        q2_coe.push_back(qq2_ceo);
        qq2_index[0] = n_v + index1;
        qq2_index[1] = n_v + index2;
        p2_index.push_back(qq2_index);
    }

    //Merge P
    for (int i = 0; i < q0.size(); ++i)
    {
        p1.push_back(q1[i]);
        p1_coe.push_back(q1_coe[i]);
    }

    for (int i = 0; i < q2.size(); ++i)
    {
        p2.push_back(q2[i]);
        p2_coe.push_back(q2_coe[i]);
    }


    // iteration
    for (int i = 0; i < p1.size(); ++i)
    {
        VertexR2 delta = p1[i];
        std::vector<double> delta_coe(p1.size(), 0.0);
        delta_coe[i] = p1_coe[i];

        update_delta_iso(p1, p1_coe, p2, p2_coe, p2_index, delta, delta_coe);
        double sum_delta = gbc::sum(delta_coe);
        gbc::div(delta_coe, sum_delta);
        gbc::add(b, delta_coe);    
    }
    double sum_b = gbc::sum(b);
    gbc::div(b, sum_b);

}

void gbc::PointwiseIterativeR2::update_delta_iso(std::vector<VertexR2>& p1, std::vector<double>& p1_coe, 
            std::vector<VertexR2>& p2,std::vector<std::vector<double>> & p2_coe, std::vector<std::vector<int>> & p2_index,
            VertexR2& delta, std::vector<double>& delta_coe)
{
    double delta_norm = delta.length();
    if (delta_norm > 1.0)
    {
        delta_norm = 1.0;
    }
    int times = 0;
    int n_vv = p1.size();
    int n_ff = p2.size();
    while (delta_norm > _tol)
    {
        times += 1;
        double alpha = std::asin(delta_norm);
        std::vector<double> rho1_coe(n_vv, 0.0);
        std::vector<double> rho2_coe(n_ff, 0.0);
        std::vector<double> d1_coe(n_vv, 0.0);
        std::vector<double> d2_coe(n_ff, 0.0);
        std::vector<double> new_delta_coe(n_vv, 0.0);

        get_rhoi_information(p1, alpha, delta, delta_coe, rho1_coe, d1_coe);
        get_rhoi_information(p2, alpha, delta, delta_coe, rho2_coe, d2_coe);
        
        double sum = 0.0;
        double sum_delta_rho = 0.0; 
        for (int i = 0; i < n_vv; ++i)
        {
            sum += rho1_coe[i];
            sum_delta_rho += rho1_coe[i] * d1_coe[i];
        }
        for (int i = 0; i < n_ff; ++i)
        {
            sum += rho2_coe[i];
            sum_delta_rho += rho2_coe[i] * d2_coe[i];
        }

        // double ip;
        // std::vector<VertexR2> d1(n_vv);
        // std::vector<VertexR2> d2(n_ff);
        // for (int i = 0; i < n_vv; ++i)
        // {
        //     ip = (i - 1 + n_vv) % n_vv;
        //     new_delta_coe[i] = sum_delta_rho * delta_coe[i] + rho1_coe[i] * (1 - d1_coe[i]) * p1_coe[i];
        //     d1[i] = rho1_coe[i] * (d1_coe[i] * delta + (1 - d1_coe[i]) * p1[i]);
        // }

        // for (int i = 0; i < n_ff; ++i)
        // {
        //     int ip1 = p2_index[i][0];
        //     int ip2 = p2_index[i][1];
        //     new_delta_coe[ip1] += rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][0];
        //     new_delta_coe[ip2] += rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][1];
        //     d2[i] = rho2_coe[i] * (d2_coe[i] * delta + (1 - d2_coe[i]) * p2[i]);
        // }

        // delta[0] = 0.0;
        // delta[1] = 0.0;
        // for (int i = 0; i < n_vv; ++i)
        // {
        //     delta_coe[i] = new_delta_coe[i] / sum;
        //     delta += d1[i]; 
        // }
        // for (int i = 0; i < n_ff; ++i)
        // {
        //     delta += d2[i]; 
        // }
        // delta = delta / sum;
        // delta_norm = delta.length();
        double ip;
        std::vector<VertexR2> d1(n_vv);
        std::vector<VertexR2> d2(n_ff);
        for (int i = 0; i < n_vv; ++i)
        {
            ip = (i - 1 + n_vv) % n_vv;
            new_delta_coe[i] = sum_delta_rho * delta_coe[i] + rho1_coe[i] * (1 - d1_coe[i]) * p1_coe[i];
            d1[i] = rho1_coe[i] * (d1_coe[i] * delta + (1 - d1_coe[i]) * p1[i]);
        }

        for (int i = 0; i < n_ff; ++i)
        {
            int ip1 = p2_index[i][0];
            int ip2 = p2_index[i][1];
            new_delta_coe[ip1] += rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][0];
            new_delta_coe[ip2] += rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][1];
            d2[i] = rho2_coe[i] * (d2_coe[i] * delta + (1 - d2_coe[i]) * p2[i]);
        }

        delta[0] = 0.0;
        delta[1] = 0.0;
        for (int i = 0; i < n_vv; ++i)
        {
            delta_coe[i] = new_delta_coe[i] / sum;
            delta += d1[i]; 
        }
        for (int i = 0; i < n_ff; ++i)
        {
            delta += d2[i]; 
        }
        delta = delta / sum;
        delta_norm = delta.length();
   }
}
/// -------------------------------对比实验-----------------------------------------
/// 对比实验1

void gbc::PointwiseIterativeR2::compute_1(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior_1(p, b);
    return;
}

void gbc::PointwiseIterativeR2::bcInterior_1(const VertexR2 &p, std::vector<double> &b)
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
    std::vector<double> p1_coe(n_v, 0.0); // Coefficients of the first projected point: 1/r_i
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        p1[i] *= 1 / tmp; // the first projected point
        p1_coe[i] += 1 / tmp;
    }

    // calutate \tilde{P}
    std::vector<VertexR2> p2(n_v); // the second projected point. 
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0)); // Coefficients of the second projected point
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
        VertexR2 delta = p1[1]; // Initialise \delta
        std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
        delta_coe[1] = p1_coe[1];

        update_delta(p1, p1_coe, p2, p2_coe, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
        double sum_delta = gbc::sum(delta_coe);
        gbc::div(delta_coe, sum_delta); // normalize
        gbc::add(b, delta_coe); 
    }
    double sum_b = gbc::sum(b);
    gbc::div(b, sum_b);
}

/// 对比实验2
void gbc::PointwiseIterativeR2::compute_1toN(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    // if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior_1toN(p, b);
    return;
}

void gbc::PointwiseIterativeR2::bcInterior_1toN(const VertexR2 &p, std::vector<double> &b)
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
    std::vector<double> p1_coe(n_v, 0.0); // Coefficients of the first projected point: 1/r_i
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        p1[i] *= 1 / tmp; // the first projected point
        p1_coe[i] += 1 / tmp;
    }

    // calutate \tilde{P}
    std::vector<VertexR2> p2(n_v); // the second projected point. 
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0)); // Coefficients of the second projected point
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
    VertexR2 delta = VertexR2(); // Initialise \delta
    std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
    for (int k = 0; k < n_v; ++k)
    {
        delta[0] += p1[k][0] / n_v;
        delta[1] += p1[k][1] / n_v;
        delta_coe[k] = p1_coe[k] / n_v;
    }

    update_delta(p1, p1_coe, p2, p2_coe, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
    double sum_delta = gbc::sum(delta_coe);
    gbc::div(delta_coe, sum_delta); // normalize
    gbc::add(b, delta_coe); 
}

/// 对比实验3
void gbc::PointwiseIterativeR2::compute_Nto2N(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior_Nto2N(p, b);
    return;
}

void gbc::PointwiseIterativeR2::bcInterior_Nto2N(const VertexR2 &p, std::vector<double> &b)
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
    std::vector<double> p1_coe(n_v, 0.0); // Coefficients of the first projected point: 1/r_i
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        p1[i] *= 1 / tmp; // the first projected point
        p1_coe[i] += 1 / tmp;
    }

    // calutate \tilde{P}
    std::vector<VertexR2> p2(n_v); // the second projected point. 
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0)); // Coefficients of the second projected point
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
        VertexR2 delta = VertexR2(); // Initialise \delta
        std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
        for (int k = 0; k < n_v; ++k)
        {
            delta[0] += p2[k][0] / n_v;
            delta[1] += p2[k][1] / n_v;
            delta_coe[k] = (p2_coe[k][0] + p2_coe[(k-1+n_v) % n_v][1]) / n_v;
        }

        update_delta(p1, p1_coe, p2, p2_coe, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
        double sum_delta = gbc::sum(delta_coe);
        gbc::div(delta_coe, sum_delta); // normalize
        gbc::add(b, delta_coe); 
    }
    double sum_b = gbc::sum(b);
    gbc::div(b, sum_b);
}

/// 对比实验4
void gbc::PointwiseIterativeR2::compute_1to2N(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior_1to2N(p, b);
    return;
}

void gbc::PointwiseIterativeR2::bcInterior_1to2N(const VertexR2 &p, std::vector<double> &b)
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
    std::vector<double> p1_coe(n_v, 0.0); // Coefficients of the first projected point: 1/r_i
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        p1[i] *= 1 / tmp; // the first projected point
        p1_coe[i] += 1 / tmp;
    }

    // calutate \tilde{P}
    std::vector<VertexR2> p2(n_v); // the second projected point. 
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0)); // Coefficients of the second projected point
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
        VertexR2 delta = VertexR2(); // Initialise \delta
        std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
        for (int k = 0; k < n_v; ++k)
        {
            delta[0] += p1[k][0] / (n_v * 2);
            delta[1] += p1[k][1] / (n_v * 2);
            delta[0] += p2[k][0] / (n_v * 2);
            delta[1] += p2[k][1] / (n_v * 2);
            delta_coe[k] += (p2_coe[k][0] + p2_coe[(k-1+n_v) % n_v][1]) / (n_v * 2);
            delta_coe[k] += p1_coe[k] / (n_v * 2);
        }

        update_delta(p1, p1_coe, p2, p2_coe, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
        double sum_delta = gbc::sum(delta_coe);
        gbc::div(delta_coe, sum_delta); // normalize
        gbc::add(b, delta_coe); 
    }
    double sum_b = gbc::sum(b);
    gbc::div(b, sum_b);
}

/// 对比实验5
void gbc::PointwiseIterativeR2::compute_p2(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior_p2(p, b);
    return;
}

void gbc::PointwiseIterativeR2::bcInterior_p2(const VertexR2 &p, std::vector<double> &b)
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
    std::vector<double> p1_coe(n_v, 0.0); // Coefficients of the first projected point: 1/r_i
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        p1[i] *= 1 / tmp; // the first projected point
        p1_coe[i] += 1 / tmp;
    }

    // calutate \tilde{P}
    std::vector<VertexR2> p2(n_v); // the second projected point. 
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0)); // Coefficients of the second projected point
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
        VertexR2 delta = p2[i]; // Initialise \delta
        std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
        delta_coe[i] = p2_coe[i][0];
        delta_coe[(i+1)%n_v] = p2_coe[i][1];

        update_delta(p1, p1_coe, p2, p2_coe, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
        double sum_delta = gbc::sum(delta_coe);
        gbc::div(delta_coe, sum_delta); // normalize
        gbc::add(b, delta_coe); 
    }
    double sum_b = gbc::sum(b);
    gbc::div(b, sum_b);
}


/// 对比实验6
void gbc::PointwiseIterativeR2::compute_p1p2(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior_p1p2(p, b);
    return;
}

void gbc::PointwiseIterativeR2::bcInterior_p1p2(const VertexR2 &p, std::vector<double> &b)
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
    std::vector<double> p1_coe(n_v, 0.0); // Coefficients of the first projected point: 1/r_i
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        p1[i] *= 1 / tmp; // the first projected point
        p1_coe[i] += 1 / tmp;
    }

    // calutate \tilde{P}
    std::vector<VertexR2> p2(n_v); // the second projected point. 
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0)); // Coefficients of the second projected point
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
        VertexR2 delta = p2[(i + n_v - 1) % n_v] / 4 + p1[i] / 2 + p2[(i + n_v) % n_v] / 4; // Initialise \delta
        std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
        delta_coe[(i + n_v - 1) % n_v] += p2_coe[(i + n_v - 1) % n_v][0] / 4;
        delta_coe[i] += p2_coe[(i + n_v - 1) % n_v][1] / 4;
        delta_coe[i] += p2_coe[(i + n_v) % n_v][0] / 4;
        delta_coe[(i + n_v + 1) % n_v] += p2_coe[(i + n_v) % n_v][1] / 4;
        delta_coe[i] += p1_coe[i] / 2;

        // VertexR2 delta = p1[(i + n_v) % n_v] / 4 + p2[i] / 2 + p1[(i + n_v + 1) % n_v] / 4; // Initialise \delta
        // std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
        // delta_coe[(i + n_v + 1) % n_v] += p2_coe[i][1] / 2;
        // delta_coe[i] += p2_coe[i][0] / 2;
        // delta_coe[i] += p1_coe[i] / 4;
        // delta_coe[(i + n_v + 1) % n_v] += p1_coe[(i + n_v + 1) % n_v] / 4;

        update_delta(p1, p1_coe, p2, p2_coe, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
        double sum_delta = gbc::sum(delta_coe);
        gbc::div(delta_coe, sum_delta); // normalize
        gbc::add(b, delta_coe); 
    }
    double sum_b = gbc::sum(b);
    gbc::div(b, sum_b);
}

