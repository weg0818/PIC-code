/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-07 14:31:48
 * @LastEditTime: 2025-03-01 12:50:40
 * @Description: 
 * @FileName: 
 * @FilePath: /Codes/C++/gbc-master/2d/extra/math.h
 */

#ifndef GBC_MATH_HPP
#define GBC_MATH_HPP

#include "VertexR2.hpp"
#include "MeshR2.hpp"
#include "../extra/libs/eigen/Eigen/Dense"
#include <vector>
namespace gbc
{
    const double PI = 3.14159265358979323846; // 自定义的 PI 常量

    double dot(VertexR2 a, VertexR2 b);

    double dot(double a, double b, double c, double d);

    double cross(VertexR2 a, VertexR2 b);

    double cross(double a, double b, double c, double d);

    VertexR2 diff(VertexR2 a, VertexR2 b);

    double sum(std::vector<double>& nums);

    void div(std::vector<double>& nums, double s);

    void add(std::vector<double>& nums1, std::vector<double>& nums2);

    void rotateSomeAngel(VertexR2 & a);

    void rotateSomeAngel(VertexR2 & a, double angle);

    // Calculate the angle abc
    double get_unsigned_angle(VertexR2 a, VertexR2 b, VertexR2 c);
    
    // Calculate the angle between the angle vectors
    double get_unsigned_angle(VertexR2 a, VertexR2 b);

    // Converted to MatrixXd of size n * 2
    void vecVertexToMatXd2(const std::vector<VertexR2> & vecVertex, Eigen::MatrixXd & MatXd);

    // Converted to MatrixXi of size n * 3
    void vecFaceToMatXi2(const std::vector<Face> & vecFace, Eigen::MatrixXi & MatXi);

    // Determine whether the line segments are intersecting, if intersecting then return true and return, if not intersecting then return flase
    bool is_intersect_segment(VertexR2 & a1, VertexR2 & a2, VertexR2 & b1, VertexR2 & b2);

    // Determine if two lines intersect
    // output:: point. If it returns true and there is an intersection, the intersection is point
    bool is_intersect_line(VertexR2 & a1, VertexR2 & a2, VertexR2 & b1, VertexR2 & b2, VertexR2 & point);

    // Determine whether P is on line AB, if so, return true, otherwise return false
    // output:: u, constant ratio splitting, AP/AB
    bool judge_point_on_lines_segment(VertexR2 & A, VertexR2 & B, VertexR2 & P, double & u);

    // Determine if the two points agree
    bool same_point(VertexR2 & A, VertexR2 & B);
}

#endif // GBC_MATH_HPP