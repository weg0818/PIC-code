#include <vector>
#include <cmath>
#include "PointwiseIterativeR3.h"
#include "../Basic/GBCmath.h"

void gbc::PointwiseIterativeR3::compute(const VertexR3 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior(p, b);
    return;
}


void gbc::PointwiseIterativeR3::bcInterior(const VertexR3 &p, std::vector<double> &b)
{
    std::vector<VertexR3> p1;
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

    std::vector<Face> face(n_f); // the faces of cage
    adjustSequenceOfFase(p1, face); // Make sure the vertices of the faces are in the same order

    // calutate \tilde{P}
    std::vector<VertexR3> p2(n_f); // the second projected point. 
    std::vector<std::vector<double>> p2_coe(n_f, std::vector<double>(3, 0.0)); // Coefficients of the second projected point
    get_p2_information(p1, p2, p2_coe, face); 

    // iteration
    // for (int i = 0; i < n_v; ++i)
    // {
    //     VertexR3 delta = p1[i]; // Initialise \delta
    //     std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
    //     delta_coe[i] = p1_coe[i];
    //     update_delta(p1, p1_coe, p2, p2_coe, face, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
    //     double sum_delta = Basic::sum(delta_coe); 
    //     Basic::div(delta_coe, sum_delta); // normalize
    //     Basic::add(b, delta_coe);    
    // }

    for (int i = 0; i < p2.size(); ++i)
    {
        VertexR3 delta = p2[i]; // Initialise \delta
        std::vector<double> delta_coe(n_v, 0.0); // Initialise the coefficients of \delta
        delta_coe[face[i].v[0]] += p2_coe[i][0] * p1_coe[face[i].v[0]];
        delta_coe[face[i].v[1]] += p2_coe[i][1] * p1_coe[face[i].v[1]];
        delta_coe[face[i].v[2]] += p2_coe[i][2] * p1_coe[face[i].v[2]];
        update_delta(p1, p1_coe, p2, p2_coe, face, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
        double sum_delta = Basic::sum(delta_coe); 
        Basic::div(delta_coe, sum_delta); // normalize
        Basic::add(b, delta_coe);    
    }
    
    double sum_b = Basic::sum(b); 
    Basic::div(b, sum_b); // get PIC
}

void gbc::PointwiseIterativeR3::update_delta(std::vector<VertexR3>& p1, std::vector<double>& p1_coe, 
            std::vector<VertexR3>& p2,std::vector<std::vector<double>>& p2_coe, std::vector<Face>& face,
            VertexR3& delta, std::vector<double>& delta_coe)
{
    int _n = p1.size();
    int _f = p2.size();
    double delta_norm = delta.length();
    if (delta_norm > 1.0)
    {
        delta_norm = 1.0;
    }
    int k = 0; 
    while (delta_norm > _tol)
    {
        k += 1;
        double alpha = std::asin(delta_norm);
        std::vector<double> rho1_coe(_n, 0.0);
        std::vector<double> rho2_coe(_f, 0.0);
        std::vector<double> d1_coe(_n, 0.0);
        std::vector<double> d2_coe(_f, 0.0);
        std::vector<double> new_delta_coe(_n, 0.0);

        get_rhoi_information(p1, alpha, delta, delta_coe, rho1_coe, d1_coe); // get the \rho and u
        get_rhoi_information(p2, alpha, delta, delta_coe, rho2_coe, d2_coe); // get the \rho and u
        
        double sum = 0.0; // \sum{\rho}
        double sum_delta_rho = 0.0; // \sum{\rho u}
        for (int i = 0; i < _n; ++i)
        {
            sum += rho1_coe[i];
            sum_delta_rho += rho1_coe[i] * d1_coe[i];
        }
        for (int i = 0; i < _f; ++i)
        {
            sum += rho2_coe[i];
            sum_delta_rho += rho2_coe[i] * d2_coe[i];
        }

        std::vector<VertexR3> d1(_n);
        std::vector<VertexR3> d2(_f);
        for (int i = 0; i < _n; ++i)
        {   
            new_delta_coe[i] += sum_delta_rho * delta_coe[i] + rho1_coe[i] * (1 - d1_coe[i]) * p1_coe[i];
            d1[i] = (delta * d1_coe[i] + p1[i] * (1 - d1_coe[i])) * rho1_coe[i];
        }
        for (int i = 0; i < _f; ++i)
        {   
            int i0 = face[i].v[0], i1 = face[i].v[1], i2 = face[i].v[2];
            new_delta_coe[i0] += rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][0] * p1_coe[i0];
            new_delta_coe[i1] += rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][1] * p1_coe[i1];
            new_delta_coe[i2] += rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][2] * p1_coe[i2];
            d2[i] = (delta * d2_coe[i] + p2[i] * (1 - d2_coe[i])) * rho2_coe[i];
        }

        delta[0] = 0.0;
        delta[1] = 0.0;
        delta[2] = 0.0;
        for (int i = 0; i < _n; ++i)
        {
            delta_coe[i] = new_delta_coe[i] / sum;
            delta += d1[i];
        }
       for (int i = 0; i < _f; ++i)
        {
            delta += d2[i];
        }
        delta = delta / sum;
        delta_norm = delta.length();
    }
}


void gbc::PointwiseIterativeR3::get_rhoi_information(std::vector<VertexR3>& pi, double alpha,
        VertexR3& delta, std::vector<double>& delta_coe,
        std::vector<double>& rhoi_coe, std::vector<double>& di_coe)
{
    double tmp, theta;
    double delta_norm = delta.length();
    for (int i = 0; i < pi.size(); ++i)
    {
        tmp = (Basic::dot(pi[i], delta)) / (pi[i].length() * delta_norm);
        if (tmp > 1.0) 
        {
            tmp = 1.0;
        }
        if (tmp < -1.0)
        {
            tmp = -1.0;
        }
        theta = std::acos(tmp);
        if (theta > Basic::PI / 2 - alpha)
        {
            rhoi_coe[i] = f1((theta - Basic::PI/2 + alpha) / (Basic::PI/2 + alpha));
            VertexR3 del_p = delta - pi[i];
            di_coe[i] = 0.0 - Basic::dot(pi[i], del_p) / Basic::dot(del_p, del_p);
        }
    }
}

void gbc::PointwiseIterativeR3::adjustSequenceOfFase(std::vector<VertexR3>& p, std::vector<Face>& face)
{
    for (int i = 0; i < n_f; i++)
    {
        const Face& face_i = _f[i];
        VertexR3 a = p[face_i.v[1]] - p[face_i.v[0]];
        VertexR3 b = p[face_i.v[2]] - p[face_i.v[0]];
        VertexR3 n = Basic::cross(a, b);
        VertexR3 c = p[face_i.v[0]];
        if (Basic::dot(n, c) < 0)
        {
            face[i] = Face(face_i.v[2], face_i.v[1], face_i.v[0]);
        }
        else
        {
            face[i] = Face(face_i.v[0], face_i.v[1], face_i.v[2]);
        }
    }
}

bool gbc::PointwiseIterativeR3::isZero(VertexR3 & a)
{
    return std::fabs(a[0]) < _tol && std::fabs(a[1]) < _tol && std::fabs(a[2]) < _tol;
}

bool gbc::PointwiseIterativeR3::get_p2_information(std::vector<VertexR3>& p1, std::vector<VertexR3> &p2, 
            std::vector<std::vector<double>> &p2_coe, std::vector<Face>& face)
{   
    for (int i = 0; i < n_f; ++i)
    {
        int i0 = face[i].v[0], i1 = face[i].v[1], i2 = face[i].v[2];

        // Calculate the unit normal, and note the direction of the cross product.
        VertexR3 n0, n1, n2;
        Basic::get_normal_face(p1[i1], p1[i2], n0);
        Basic::get_normal_face(p1[i2], p1[i0], n1);
        Basic::get_normal_face(p1[i0], p1[i1], n2);

        // If three vertices of a face are on a straight line, the coefficient is set to zero.
        double ab = (p1[i0] - p1[i1]).length();
        double bc = (p1[i1] - p1[i2]).length();
        double ac = (p1[i0] - p1[i2]).length();
        double abc = (ab + bc + ac) / 2;
        double area = squaredArea(ab, bc, ac, abc);
        if (area < _tol)
        {
            continue;
        }

        // Getting the angle
        double b0, b1, b2;
        b0 = std::acos(Basic::dot(p1[i1], p1[i2]));
        b1 = std::acos(Basic::dot(p1[i2], p1[i0]));
        b2 = std::acos(Basic::dot(p1[i0], p1[i1]));

        // Calculate the coefficients of normal
        double c0 = (b0 + b1 * Basic::dot(n1, n0) + b2 * Basic::dot(n2, n0)) / (2 * Basic::dot(p1[i0], n0));
        double c1 = (b1 + b2 * Basic::dot(n2, n1) + b0 * Basic::dot(n0, n1)) / (2 * Basic::dot(p1[i1], n1));
        double c2 = (b2 + b0 * Basic::dot(n0, n2) + b1 * Basic::dot(n1, n2)) / (2 * Basic::dot(p1[i2], n2));
        if (fabs(c0) < 1e-6 && fabs(c1) < 1e-6 && fabs(c2) < 1e-6) // If all coefficients are zero
        {
            continue;
        }
        if (c0 < 0 || c1 < 0 || c2 < 0)
        {
            std::cout << "Error: p2 is incorrectly calculated, with a negative coefficient." << c0 << ' ' << c1 << ' ' << c2 << std::endl;
        }
        p2_coe[i][0] = c0;
        p2_coe[i][1] = c1;
        p2_coe[i][2] = c2;
        p2[i] = p1[i0] * p2_coe[i][0] + p1[i1] * p2_coe[i][1] +  p1[i2] * p2_coe[i][2];
       
        // unitization
        double p2_norm = p2[i].length();
        p2_coe[i][0] /= p2_norm;
        p2_coe[i][1] /= p2_norm;
        p2_coe[i][2] /= p2_norm;

        p2[i] /= p2_norm;
    }
    return true;

}


void gbc::PointwiseIterativeR3::compute_iso_pic(const VertexR3 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v + _iso_p_n, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior_iso(p, b);
    return;
}

void gbc::PointwiseIterativeR3::bcInterior_iso(const VertexR3 &p, std::vector<double> &b)
{

    std::vector<VertexR3> p1;
    for (int i = 0; i < n_v; ++i)
    {
        p1.push_back(_v[i]);
    }

    // translation
    for (auto &p1_i : p1)
    {
        p1_i -= p;
    }

    std::vector<VertexR3> q1 = _iso_p;
    for (auto & item : q1)
    {
        item -= p;
    }

    // calutate \hat{P}
    std::vector<double> p1_coe(n_v, 0.0); // Coefficients of the first projected point: 1/r_i
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        if (tmp < 1e-13)
        {
            b[i] = 1.0;
            return;
        }
        p1[i] *= 1 / tmp;
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

    std::vector<Face> face(n_f); // the faces of cage
    adjustSequenceOfFase(p1, face); // Make sure the vertices of the faces are in the same order 

    // calutate \tilde{P}
    std::vector<VertexR3> p2(n_f);
    std::vector<std::vector<double>> p2_coe(n_f, std::vector<double>(3, 0.0));
    get_p2_information(p1, p2, p2_coe, face);

    //Merge
    for (int i = 0; i < q1.size(); ++i)
    {
        p1.push_back(q1[i]);
        p1_coe.push_back(q1_coe[i]);
        
    }
    // iteration
    for (int i = 0; i < p1.size(); ++i)
    {
        VertexR3 delta = p1[i]; // Initialise \delta
        std::vector<double> delta_coe(p1.size(), 0.0); // Initialise the coefficients of \delta
        delta_coe[i] = p1_coe[i];
        update_delta(p1, p1_coe, p2, p2_coe, face, delta, delta_coe); // to get the homogeneous coordinates of PIC with respect to the vertex p1_i
        double sum_delta = Basic::sum(delta_coe);
        Basic::div(delta_coe, sum_delta); // normalize
        Basic::add(b, delta_coe);    
    }
    
    double sum_b = Basic::sum(b);
    Basic::div(b, sum_b); // get PIC
}