#include "GBCmath.h"
#include <omp.h>
#include <igl/barycentric_coordinates.h>
#include <igl/EPS.h>
#include <igl/lbs_matrix.h>

double Basic::dot(VertexR3 &a, VertexR3 &b)
{
    return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

VertexR3 Basic::cross(VertexR3 &a, VertexR3 &b)
{
    return VertexR3(
        a.y() * b.z() - a.z() * b.y(),
        a.z() * b.x() - a.x() * b.z(),
        a.x() * b.y() - a.y() * b.x()
    );
}

// Get the unit normal direction of the face
//return n
void Basic::get_normal_face(VertexR3 &a, VertexR3 &b, VertexR3 &n)
{
    n = cross(a, b);
    n /= n.length();
}

// 计算三点的法向量
VertexR3 Basic::getNormal(const VertexR3& p1, const VertexR3& p2, const VertexR3& p3) {
    VertexR3 v1(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
    VertexR3 v2(p3.x() - p1.x(), p3.y() - p1.y(), p3.z() - p1.z());
    return Basic::cross(v1, v2);
}

// vector summation
double Basic::sum(std::vector<double> & nums)
{
    double sum_n = 0.0;
    for (int i = 0; i < nums.size(); ++i)
    {
        sum_n += nums[i];
    }
    return sum_n;
}

void Basic::div(std::vector<double> & nums, double a)
{
    for (int i = 0; i < nums.size(); ++i)
    {
        nums[i] /= a;
    }
}

void Basic::add(std::vector<double> & nums1, std::vector<double> & nums2)
{
    for (int i = 0; i < nums1.size(); ++i)
    {
        nums1[i] += nums2[i];
    }
}

void Basic::vecVertexToMatXd(const std::vector<VertexR3> & vecVertex, Eigen::MatrixXd & MatXd)
{
    MatXd = Eigen::MatrixXd::Zero(3, vecVertex.size());
    for (int i = 0; i < vecVertex.size(); ++i)
    {
        MatXd(0, i) = vecVertex[i][0];
        MatXd(1, i) = vecVertex[i][1];
        MatXd(2, i) = vecVertex[i][2];
    }
}

void Basic::vecVertexToMatXd2(const std::vector<VertexR3> & vecVertex, Eigen::MatrixXd & MatXd)
{
    MatXd = Eigen::MatrixXd::Zero(vecVertex.size(), 3);
    for (int i = 0; i < vecVertex.size(); ++i)
    {
        MatXd(i, 0) = vecVertex[i][0];
        MatXd(i, 1) = vecVertex[i][1];
        MatXd(i, 2) = vecVertex[i][2];
    }
}

void Basic::vecFaceToMatXi(const std::vector<Face> & vecFace, Eigen::MatrixXi & MatXi)
{
    MatXi = Eigen::MatrixXi::Zero(3, vecFace.size());
    for (int i = 0; i < vecFace.size(); ++i)
    {
        MatXi(0, i) = vecFace[i].v[0];
        MatXi(1, i) = vecFace[i].v[1];
        MatXi(2, i) = vecFace[i].v[2];
    }
}

void Basic::vecFaceToMatXi2(const std::vector<Face> & vecFace, Eigen::MatrixXi & MatXi)
{
    MatXi = Eigen::MatrixXi::Zero(vecFace.size(), 3);
    for (int i = 0; i < vecFace.size(); ++i)
    {
        MatXi(i, 0) = vecFace[i].v[0];
        MatXi(i, 1) = vecFace[i].v[1];
        MatXi(i, 2) = vecFace[i].v[2];
    }
}


void Basic::MatXdToVecVertex(const Eigen::MatrixXd & MatXd, std::vector<VertexR3> & vecVertex)
{
    int n = MatXd.rows();
    for (int i = 0; i < n; ++i)
    {
        vecVertex.push_back(VertexR3(MatXd(i, 0), MatXd(i, 1), MatXd(i, 2)));
    }
}

// ajust the sequance of index of the Cells ，底面三个点右手系，最后一个点在上面
void Basic::adjustCellIndex(const std::vector<VertexR3> & vecVertex, Eigen::MatrixXi & cells)
{
    int n = cells.rows();
    for (int i = 0; i < n; ++i)
    {
        std::vector<int> cell = {cells(i, 0), cells(i, 1), cells(i, 2), cells(i, 3)};
        // 按照z值排序
        std::sort(cell.begin(), cell.end(), [&vecVertex](int a, int b) {
            return vecVertex[a][2] < vecVertex[a][2];
        });
        // 前三个点是底面，第四个点是顶点
        int p1_idx = cell[0];
        int p2_idx = cell[1];
        int p3_idx = cell[2];
        int p4_idx = cell[3];
        // 确保底面三个点构成右手坐标系
        if (!isRightHanded(vecVertex[p1_idx], vecVertex[p2_idx], vecVertex[p3_idx])) {
            // 如果不符合右手坐标系，交换p2和p3
            std::swap(p2_idx, p3_idx);
        }

        cells(i, 0) = p1_idx;
        cells(i, 1) = p2_idx;
        cells(i, 2) = p3_idx;
        cells(i, 3) = p4_idx;
    }
}

// 判断点是否在右手坐标系中
bool Basic::isRightHanded(const VertexR3& p1, const VertexR3& p2, const VertexR3& p3) {
    VertexR3 normal = getNormal(p1, p2, p3);
    return normal.z() > 0;  // 如果法向量z值为正，则符合右手坐标系
}


void Basic::interpolateWeightsInEmbedding(const Eigen::MatrixXd& V, const Eigen::MatrixXd& W, const Eigen::MatrixXd& V_embedding,
	const Eigen::MatrixXi& T_embedding, int numCageVertices, Eigen::MatrixXd& W_interpolated, Eigen::MatrixXd& M)
{
	W_interpolated.resize(V.rows(), numCageVertices);
#pragma omp parallel for
	for (int i = 0; i < V.rows(); ++i)
	{
		Eigen::RowVector3d vert = V.row(i);
		bool found = false;

		// Find containing tetrahedra
		// TODO - More efficient approach than linear search
		for (int j = 0; j < T_embedding.rows(); ++j)
		{
			const Eigen::Vector4i tet = T_embedding.row(j);
			const Eigen::RowVector3d v_0 = V_embedding.row(tet(0));
			const Eigen::RowVector3d v_1 = V_embedding.row(tet(1));
			const Eigen::RowVector3d v_2 = V_embedding.row(tet(2));
			const Eigen::RowVector3d v_3 = V_embedding.row(tet(3));
			Eigen::RowVector4d barycentrics;
			found = true;

			igl::barycentric_coordinates(vert, v_0, v_1, v_2, v_3, barycentrics);


			for (int k = 0; k < 4u; ++k)
			{
				if (0. - igl::FLOAT_EPS > barycentrics(k) || 1. + igl::FLOAT_EPS < barycentrics(k))
				{
					found = false;
					break;
				}
			}

			if (!found)
			{
				continue;
			}

			W_interpolated.row(i) = W.row(tet(0)) * barycentrics(0) + W.row(tet(1)) * barycentrics(1) + W.row(tet(2)) * barycentrics(2) +
				W.row(tet(3)) * barycentrics(3);
			break;
		}
		assert(found);
		//std::cout << (static_cast<float>(i) / static_cast<float>(V.rows())) * 100.f << "%\n";
	}

	igl::lbs_matrix(V, W_interpolated, M);
}


