/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-19 10:03:06
 * @LastEditTime: 2025-03-06 21:55:52
 * @Description: 
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/BarycentricCoordinates/MeanValueR3.cpp
 */
#include <vector>
#include "MeanValueR3.h"
#include "../Basic/GBCmath.h"
 

void gbc::MeanValueR3::compute(const VertexR3 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior(p, b);
    return;
}


void gbc::MeanValueR3::bcInterior(const VertexR3 &p, std::vector<double> &b)
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
        p0.push_back(p1_i);
    }

    std::vector<Face> face(n_f);
    for (int i = 0; i < n_f; i++)
    {
        const Face& face_i = _f[i];
        face[i] = Face(face_i.v[0], face_i.v[1], face_i.v[2]);
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
    std::vector<std::vector<double>> p2_coe(n_f, std::vector<double>(3, 0.0));
    get_p2_information(p1, p2_coe, face);

    // get MVC
    for (int i = 0; i < n_f; ++i)
    {
        int i0 = face[i].v[0], i1 = face[i].v[1], i2 = face[i].v[2];
        b[i0] += p2_coe[i][0] * p1_coe[i0];
        b[i1] += p2_coe[i][1] * p1_coe[i1];
        b[i2] += p2_coe[i][2] * p1_coe[i2];
    }

    double sum_b = Basic::sum(b);
    Basic::div(b, sum_b);

}

void gbc::MeanValueR3::get_p2_information(std::vector<VertexR3>& p1, std::vector<std::vector<double>> &p2_coe, std::vector<Face>& face)
{
    for (int i = 0; i < n_f; ++i)
    {
        int i0 = face[i].v[0], i1 = face[i].v[1], i2 = face[i].v[2];

        // Calculate the unit normal, and note the direction of the cross product.
        VertexR3 n0, n1, n2;
        Basic::get_normal_face(p1[i1], p1[i2], n0);
        Basic::get_normal_face(p1[i2], p1[i0], n1);
        Basic::get_normal_face(p1[i0], p1[i1], n2);

        // Getting the angle
        double b0, b1, b2;
        b0 = std::acos(Basic::dot(p1[i1], p1[i2]));
        b1 = std::acos(Basic::dot(p1[i2], p1[i0]));
        b2 = std::acos(Basic::dot(p1[i0], p1[i1]));

        // Calculate the coefficients of normal
        p2_coe[i][0] = (b0 + b1 * Basic::dot(n1, n0) + b2 * Basic::dot(n2, n0)) / (2 * Basic::dot(p1[i0], n0));
        p2_coe[i][1] = (b1 + b2 * Basic::dot(n2, n1) + b0 * Basic::dot(n0, n1)) / (2 * Basic::dot(p1[i1], n1));
        p2_coe[i][2] = (b2 + b0 * Basic::dot(n0, n2) + b1 * Basic::dot(n1, n2)) / (2 * Basic::dot(p1[i2], n2));
    }
    
}
