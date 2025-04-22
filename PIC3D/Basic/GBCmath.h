/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-08 15:06:13
 * @LastEditTime: 2025-02-21 12:26:17
 * @Description: 
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/Basic/GBCmath.h
 */
#ifndef GBC_MATH_H
#define GBC_MATH_H
#include "VertexR3.hpp"
#include "Face.hpp"

#include <vector>
#include <Eigen/Dense>

namespace Basic
{
    const double PI = 3.14159265358979323846; // 自定义的 PI 常量

    // 3D vetex dot 
    double dot(VertexR3 &a, VertexR3 &b);

    VertexR3 cross(VertexR3 &a, VertexR3 &b);

    // Get the unit normal direction of the face
    //return n
    void get_normal_face(VertexR3 &a, VertexR3 &b, VertexR3 &n);

    VertexR3 getNormal(const VertexR3& p1, const VertexR3& p2, const VertexR3& p3);

    // vector summation
    double sum(std::vector<double> & nums);

    void div(std::vector<double> & nums, double a);

    void add(std::vector<double> & nums1, std::vector<double> & nums2);

    // Converted to MatrixXd of size 3 * n
    void vecVertexToMatXd(const std::vector<VertexR3> & vecVertex, Eigen::MatrixXd & MatXd);
    
    // Converted to MatrixXd of size n * 3
    void vecVertexToMatXd2(const std::vector<VertexR3> & vecVertex, Eigen::MatrixXd & MatXd);

    // Converted to MatrixXi of size 3 * n
    void vecFaceToMatXi(const std::vector<Face> & vecFace, Eigen::MatrixXi & MatXi);

    // Transforms to MatrixXi of size n * 3
    void vecFaceToMatXi2(const std::vector<Face> & vecFace, Eigen::MatrixXi & MatXi);

    // Converted to std::vector<VertexR3> 
    void MatXdToVecVertex(const Eigen::MatrixXd & MatXd, std::vector<VertexR3> & vecVertex);

    // ajust the sequance of index of the Cells ，底面三个点右手系，最后一个点在上面
    void adjustCellIndex(const std::vector<VertexR3> & vecVertex, Eigen::MatrixXi & cells);

    bool isRightHanded(const VertexR3& p1, const VertexR3& p2, const VertexR3& p3);

    void interpolateWeightsInEmbedding(const Eigen::MatrixXd& V, const Eigen::MatrixXd& W, const Eigen::MatrixXd& V_embedding, const Eigen::MatrixXi& T_embedding, int numCageVertices, Eigen::MatrixXd& W_interpolated, Eigen::MatrixXd& M);


}

#endif // GBC_MATH_H