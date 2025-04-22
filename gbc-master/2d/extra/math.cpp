#include "math.h"
using namespace gbc;

double gbc::dot(VertexR2 a, VertexR2 b)
{
    return a[0] * b[0] + a[1] * b[1];
}

double gbc::dot(double a, double b, double c, double d)
{
    return a * c + b * d;
}

double gbc::cross(VertexR2 a, VertexR2 b)
{
    return a[0] * b[1] - a[1] * b[0];
}

double gbc::cross(double a, double b, double c, double d)
{
    return a * d - b * c;
}

VertexR2 gbc::diff(VertexR2 a, VertexR2 b)
{
    return VertexR2(a[0] - b[0], a[1] - b[1]);
}

double gbc::sum(std::vector<double>& nums)
{
    double ret = 0.0;
    for (int i = 0; i < nums.size(); ++i)
    {
        ret += nums[i];
    }
    return ret;
}

void gbc::div(std::vector<double>& nums, double s)
{
    for (int i = 0; i < nums.size(); ++i)
    {
        nums[i] /= s;
    }
}

void gbc::add(std::vector<double>& nums1, std::vector<double>& nums2)
{
    for (int i = 0; i < nums1.size(); ++i)
    {
        nums1[i] += nums2[i];
    }
}

// Rotate 90 degrees
void gbc::rotateSomeAngel(VertexR2 & a)
{
    double tmp = a[0];
    a[0] = -a[1];
    a[1] = tmp;
}

void gbc::rotateSomeAngel(VertexR2 & a, double angle)
{
    VertexR2 b;
    b[0] = cos(angle) * a[0] - sin(angle) * a[1];
    b[1] = sin(angle) * a[0] + cos(angle) * a[1];
    a[0] = b[0];
    a[1] = b[1];
}

double gbc::get_unsigned_angle(VertexR2 a, VertexR2 b, VertexR2 c)
{
    return get_unsigned_angle(a - b, c - b);
}


double gbc::get_unsigned_angle(VertexR2 a, VertexR2 b)
{
    double _dot = dot(a, b);
    double magA = a.length();
    double magB = b.length();

    double cosTheta = _dot / (magA * magB);

    // Calculation of angle (radians)
    return std::acos(cosTheta);
}

void gbc::vecVertexToMatXd2(const std::vector<VertexR2> & vecVertex, Eigen::MatrixXd & MatXd)
{
    MatXd = Eigen::MatrixXd::Zero(vecVertex.size(), 2);
    for (int i = 0; i < vecVertex.size(); ++i)
    {
        MatXd(i, 0) = vecVertex[i][0];
        MatXd(i, 1) = vecVertex[i][1];
    }
}

void gbc::vecFaceToMatXi2(const std::vector<Face> & vecFace, Eigen::MatrixXi & MatXi)
{
    MatXi = Eigen::MatrixXi::Zero(vecFace.size(), 3);
    for (int i = 0; i < vecFace.size(); ++i)
    {
        MatXi(i, 0) = vecFace[i].v[0];
        MatXi(i, 1) = vecFace[i].v[1];
        MatXi(i, 2) = vecFace[i].v[2];
    }   
}

bool gbc::is_intersect_segment(VertexR2 & a1, VertexR2 & a2, VertexR2 & b1, VertexR2 & b2)
{
    // Determine whether the line segments are intersecting, if intersecting then return true and return, if not intersecting then return flase
    if (std::min(a1[0], a2[0]) > std::max(b1[0], b2[0]) || 
        std::min(a1[1], a2[1]) > std::max(b1[1], b2[1]) ||
        std::min(b1[0], b2[0]) > std::max(a1[0], a2[0]) || 
        std::min(b1[1], b2[1]) > std::max(a1[1], a2[1]))
    {
        // Failed the rapid rejection test
        return false;
    }
    double tmp1 = cross(b2 - b1, a1 - b1) * cross(b2 - b1, a2 - b1);
    double tmp2 = cross(a2 - a1, b1 - a1) * cross(a2 - a1, b2 - a1);
    if (tmp1 > 0 || tmp2 > 0)
    {
        // Failure to pass the straddle experiment
        return false;
    }
    else if (std::abs(tmp1) < 1e-13 && std::abs(tmp2) < 1e-13)
    {
        return false;
    }
    return true;
}

bool gbc::is_intersect_line(VertexR2 & a1, VertexR2 & a2, VertexR2 & b1, VertexR2 & b2, VertexR2 & point)
{
    // Determine whether two lines intersect
    double tol = 1e-10;
    VertexR2 a = a2 - a1;
    VertexR2 b = b2 - b1;
    double tmp1 = cross(a, b);
    if (std::abs(tmp1) < tol)
    {
        // Two lines are parallel
        return false;
    }
    VertexR2 c = a1 - b2;
    VertexR2 d = a2 - b2;
    double tmp2 = cross(c, d);
    double t; 
    if (tmp2 == 0)
    {
        point[0] = b2[0];
        point[1] = b2[1];
        return true;
    }
    else
    {
        VertexR2 e = a2 - b1;
        VertexR2 f = a1 - b1;
        double tmp3 = cross(e, f);
        t = tmp3 / tmp2;
    }
    double _lambda = t / (1 + t);
    VertexR2 p = (1 - _lambda) * b1 + _lambda * b2;
    point[0] = p[0];
    point[1] = p[1];
    return true;
}

bool gbc::judge_point_on_lines_segment(VertexR2 & A, VertexR2 & B, VertexR2 & P, double & u)
{
    // Determine whether P is on line AB
    double tol = 1e-10;
    if ((std::abs((P -  A)[0]) < tol && std::abs((P -  A)[1]) < tol) ||
        (std::abs((P -  B)[0]) < tol && std::abs((P -  B)[1]) < tol) )
    {
        return false;
    }
    VertexR2 AP = P - A;
    VertexR2 BP = P - B;
    // First determine if the point is on a straight line
    if (std::abs(cross(AP, BP)) < tol)
    {
        // If the point is on the line, then determine if the point is on the segment
        double min_x = std::min(A[0], B[0]);
        double max_x = std::max(A[0], B[0]);
        double min_y = std::min(A[1], B[1]);
        double max_y = std::max(A[1], B[1]);
        if (min_x - tol <= P[0] && P[0] <= max_x + tol &&
            min_y - tol <= P[1] && P[1] <= max_y + tol)
        {
            u = AP.length() / (B - A).length();
            return true;
        }
        else
        {
            return false;
        }
    }
    return false;
}

bool gbc::same_point(VertexR2 & A, VertexR2 & B)
{
    double tol = 1e-10;
    if (std::abs((A - B)[0]) < tol && std::abs((A - B)[1]) < tol)
    {
        return true;
    }
    return false;
}
