#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <algorithm> // For std::replace
#include <Eigen/Dense>
#include "igl/copyleft/tetgen/tetrahedralize.h"
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
// Local includes.
#include "./BarycentricCoordinates/PointwiseIterativeR3.h"
#include "./BarycentricCoordinates/MeanValueR3.h"
#include "./BarycentricCoordinates/MaximumEntropyCoordinatesR3.h"
#include "./BarycentricCoordinates/MaximumLikelihoodCoordinatesR3.h"
#include "./BarycentricCoordinates/LocalBarycentricCoordinatesR3.h"
#include "./BarycentricCoordinates/HarmonicR3.h"
#include "./BarycentricCoordinates/TestBCR3.hpp"
#include "./Basic/VertexR3.hpp"
#include "./Basic/Face.hpp"
#include "./Basic/GBCmath.h"
#include "./Basic/LoadandSave.h"

using namespace gbc;

void deformed(std::vector<VertexR3> & vertices_mesh, std::vector<VertexR3> & vertices_deformed_cage, std::vector<VertexR3> & new_vertices_mesh)
{
    for (int i = 0; i < vertices_mesh.size(); ++i)
    {
        std::vector<double> bc = vertices_mesh[i].b();
        for (int j = 0; j < bc.size(); ++j)
        {
            new_vertices_mesh[i] += vertices_deformed_cage[j] * bc[j];
        }
    }
}

// Examples.
int main(int argc, char* argv[]) {
    std::cout << "Loading data" << std::endl;
    std::string name = "ogre";
    // std::string name = "beast";
    // std::string name = "Cactus";
    // std::string name = "Beast";
    // const std::string path_cage = "/home/mjp/Codes/C++/PIC3D/data/" + name + "_Cage.obj";
    // const std::string path_mesh = "/home/mjp/Codes/C++/PIC3D/data/" + name + ".obj";
    // const std::string save_path = "/home/mjp/Codes/C++/PIC3D/out/OUTPUT/" + name + "/";

    const std::string path_cage = "/home/mjp/Codes/大论文多边形数据/八面体.obj";
    const std::string path_mesh = "/home/mjp/Codes/大论文多边形数据/八面体截面.obj";
    const std::string save_path = "/home/mjp/Codes/大论文多边形数据/";

    std::vector<VertexR3> vertices_cage;
    std::vector<Face> faces_cage; // index: 0-based
    Basic::loadOBJ(path_cage, vertices_cage, faces_cage);

    std::vector<VertexR3> vertices_mesh;
    std::vector<Face> faces_mesh;
    Basic::loadOBJ(path_mesh, vertices_mesh, faces_mesh);

    std::vector<VertexR3> vertices_deformed_cage;
    std::vector<Face> faces_deformed_cage;
    
    // std::string path_deformed_cage = "/home/mjp/Codes/C++/PIC3D/data/" + name + "_Cage_Deformed.obj";
    // Basic::loadOBJ(path_deformed_cage, vertices_deformed_cage, faces_deformed_cage);
    std::cout << "Number of cage vertex: " << vertices_cage.size() << std::endl;
    std::cout << "Number of Cage facet: " << faces_cage.size() << std::endl;
    std::cout << "Number of sampling point: " << vertices_mesh.size() << std::endl;
    std::cout << "Number of Sampling point facet: " << faces_mesh.size() << std::endl;
    auto last = std::chrono::high_resolution_clock::now();
    auto cur = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;

    // for (auto & v : vertices_cage)
    // {
    //     v[0] = v[0] * 1000;
    //     v[1] = v[1] * 1000;
    //     v[2] = v[2] * 1000;
    // }
    // for (auto & v : vertices_mesh)
    // {
    //     v[0] = v[0] * 1000;
    //     v[1] = v[1] * 1000;
    //     v[2] = v[2] * 1000;
    // }


    //PIC
    // std::cout << "Loading data" << std::endl;
    // std::string name = "ConvexPolyhedron";
    // const std::string path_cage = "/home/mjp/Codes/C++/PIC3D/data/" + name + ".obj";
    // const std::string save_path = "/home/mjp/Codes/C++/PIC3D/out/OUTPUT/" + name + "/";

    // std::vector<VertexR3> vertices_cage;
    // std::vector<Face> faces_cage; // index: 0-based
    // Basic::loadOBJ(path_cage, vertices_cage, faces_cage);
    // if (vertices_cage.size() == 0) 
    // {
    //     std::cout << "数据读取失败" << std::endl;
    //     return 0;
    // }

    // auto last = std::chrono::high_resolution_clock::now();
    // auto cur = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed;

    // // tetrahedral section
    // Eigen::MatrixXd C; // vertex matrix  (N x 3)
    // Eigen::MatrixXi CF; // faces matrix  (M x 3) 面索引从0开始
    // Basic::vecVertexToMatXd2(vertices_cage, C);
    // Basic::vecFaceToMatXi2(faces_cage, CF);

    // Eigen::MatrixXd V;  // The dissected vertex matrix
    // Eigen::MatrixXi TF;  // TF by 3 list of triangle face indices
    // Eigen::MatrixXi T; // TT by 4 list of tet face indices
    // //  select TetGen's Parameters
    // std::string mode1 = argv[1]; // Beast:"pq1.07a0.007"
    // // Calling TetGen for tetrahedral dissection using libigl
    // igl::copyleft::tetgen::tetrahedralize(C, CF, mode1, V, T, TF);

    // std::cout << "Successful tetrahedral sectioning" << std::endl;
    // std::cout << " Vertex:  "<<  V.rows() << ",  tetrahedral: " << T.rows() << std::endl;

    // std::vector<VertexR3> vertices_mesh;
    // Basic::MatXdToVecVertex(V, vertices_mesh);

    // Basic::adjustCellIndex(vertices_mesh, T);
    // if (Basic::save_cells(save_path + "Cells.txt", T))
    // {
    //     std::cout << "save cells successfull" << std::endl;
    // }
    // if (Basic::save_points(save_path + "Points.txt", vertices_mesh))
    // {
    //     std::cout << "save points successfull" << std::endl;
    // }

    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR3 pic = PointwiseIterativeR3(vertices_cage, faces_cage);
    // pic.bc(vertices_mesh);
    // if (Basic::save_BC(save_path + "pic_lambda", vertices_mesh))
    // {
    //     std::cout << "pic calculation successfull" << std::endl;
    // }
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "pic-times" << elapsed.count() << " s" << std::endl;

    // TestBCR3 test(vertices_cage, vertices_mesh);


    // //mvc
    last = std::chrono::high_resolution_clock::now();
    MeanValueR3 mvc = MeanValueR3(vertices_cage, faces_cage);
    mvc.bc(vertices_mesh);
    if (Basic::save_BC(save_path + name + "_mvc_lambda", vertices_mesh))
    {
        std::cout << "mvc calculation successfull" << std::endl;
    }
    cur = std::chrono::high_resolution_clock::now();
    elapsed = cur - last; 
    std::cout << "mvc-times: " << elapsed.count() << " s" << std::endl; 
    TestBCR3 test_mvc(vertices_cage, vertices_mesh);

    if (Basic::saveOBJ_contour(save_path + name + "八边形截面等值线.obj", vertices_mesh, faces_mesh, 3))
    {
        std::cout << "mvc contour calculation successfull" << std::endl;
    }

    // std::vector<VertexR3> new_vertices_mesh_mvc(vertices_mesh.size());
    // deformed(vertices_mesh, vertices_deformed_cage, new_vertices_mesh_mvc);

    // if (Basic::saveOBJ(save_path + name + "_mvc.obj", new_vertices_mesh_mvc, faces_mesh))
    // {
    //     std::cout << "mvc saves the deformation result" << std::endl;
    // }

    //mec
    // last = std::chrono::high_resolution_clock::now();
    // MaximumEntropyCoordinatesR3 mec = MaximumEntropyCoordinatesR3(vertices_cage, faces_cage);
    // mec.set_mec_flag(1);
    // mec.bc(vertices_mesh);
    // if (Basic::save_BC(save_path + name + "_mec_lambda", vertices_mesh))
    // {
    //     std::cout << "meccalculation successfull" << std::endl;
    // }
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;  
    // std::cout << "mec-times:  " << elapsed.count() << " s" << std::endl;
    // TestBCR3 test_mec(vertices_cage, vertices_mesh);

    // std::vector<VertexR3> new_vertices_mesh_mec(vertices_mesh.size());
    // deformed(vertices_mesh, vertices_deformed_cage, new_vertices_mesh_mec);

    // if (Basic::saveOBJ(save_path + name + "_mec.obj", new_vertices_mesh_mec, faces_mesh))
    // {
    //     std::cout << "mec saves the deformation result" << std::endl;
    // }

    //mlc
    // last = std::chrono::high_resolution_clock::now();
    // MaximumLikelihoodCoordinatesR3 mlc = MaximumLikelihoodCoordinatesR3(vertices_cage, faces_cage);
    // mlc.bc(vertices_mesh);
    // if (Basic::save_BC(save_path + name + "_mlc_lambda", vertices_mesh))
    // {
    //     std::cout << "mlccalculation successfull" << std::endl;
    // }
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;  
    // std::cout << "mlc-times: " << elapsed.count() << " s" << std::endl; 
    // TestBCR3 test_mlc(vertices_cage, vertices_mesh);

    // std::vector<VertexR3> new_vertices_mesh_mlc(vertices_mesh.size());
    // deformed(vertices_mesh, vertices_deformed_cage, new_vertices_mesh_mlc);

    // if (Basic::saveOBJ(save_path + name + "_mlc.obj", new_vertices_mesh_mlc, faces_mesh))
    // {
    //     std::cout << "mlc saves the deformation result" << std::endl;
    // }


    // // LBC
    // tetrahedral section
    // Eigen::MatrixXd C; // vertex matrix  (N x 3)
    // Eigen::MatrixXi CF; // faces matrix  (M x 3) 面索引从0开始
    // Basic::vecVertexToMatXd2(vertices_cage, C);
    // Basic::vecFaceToMatXi2(faces_cage, CF);

    // Eigen::MatrixXd V;  // The dissected vertex matrix
    // Eigen::MatrixXi TF;  // TF by 3 list of triangle face indices
    // Eigen::MatrixXi T; // TT by 4 list of tet face indices
    // //  select TetGen's Parameters
    // // std::string mode1 = "pq1.07a0.07"; // Ogre:"pq1.07a0.07"
    // // std::string mode1 = "pq1.07a0.0007"; // Cactus:"pq1.07a0.0007"
    // std::string mode1 = argv[1]; // Beast:"pq1.07a0.007"
    // // std::string mode1 = "pq2a0.005"; 
    // // Calling TetGen for tetrahedral dissection using libigl
    // igl::copyleft::tetgen::tetrahedralize(C, CF, mode1, V, T, TF);

    // std::cout << "Successful tetrahedral sectioning" << std::endl;
    // std::cout << " Vertex:  "<<  V.rows() << ",  tetrahedral: " << T.rows() << std::endl;

    // last = std::chrono::high_resolution_clock::now();
    // LocalBarycentricCoordinatesR3 lbc = LocalBarycentricCoordinatesR3(vertices_cage, faces_cage);
    // lbc.Init(V, T, C, CF);
    // lbc.set_lbc_scheme(2, 9000);
    // lbc.bc(vertices_mesh);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;  
    // std::cout << "lbc-times:  " << elapsed.count() << " s" << std::endl;
    // if (Basic::save_BC(save_path +  name + "_lbc_lambda_" + mode1, vertices_mesh))
    // {
    //     std::cout << "lbc calculation successfull" << std::endl;
    // }
    // TestBCR3 test_lbc(vertices_cage, vertices_mesh);

    // std::vector<VertexR3> new_vertices_mesh_lbc(vertices_mesh.size());
    // deformed(vertices_mesh, vertices_deformed_cage, new_vertices_mesh_lbc);

    // if (Basic::saveOBJ(save_path + name + "_lbc.obj", new_vertices_mesh_lbc, faces_mesh))
    // {
    //     std::cout << "lbc saves the deformation result" << std::endl;
    // }

    // HC
    // last = std::chrono::high_resolution_clock::now();
    // HarmonicR3 hc = HarmonicR3(vertices_cage, faces_cage);
    // hc.Init(V, T, C, CF);
    // hc.bc(vertices_mesh);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last; 
    // std::cout << "hc-times: " << elapsed.count() << " s" << std::endl;

    // if (Basic::save_BC(save_path +  name + "_hc_lambda_" + mode1, vertices_mesh))
    // {
    //     std::cout << "hc calculation successfull" << std::endl;
    // }
    // TestBCR3 test_hc(vertices_cage, vertices_mesh);

    // std::vector<VertexR3> new_vertices_mesh_hc(vertices_mesh.size());
    // deformed(vertices_mesh, vertices_deformed_cage, new_vertices_mesh_hc);

    // if (Basic::saveOBJ(save_path + name + "_hc.obj", new_vertices_mesh_hc, faces_mesh))
    // {
    //     std::cout << "hc saves the deformation result" << std::endl;
    // }


    // //PIC  Isolated
    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR3 pic = PointwiseIterativeR3(vertices_cage, faces_cage);
    // std::vector<VertexR3> iso(1);
    // iso[0][0] = 6.06463243486742698;
    // iso[0][1] = 1.354801144931828905;
    // iso[0][2] = -944.581787283142944;
    // pic.set_iso_p(iso);
    // pic.compute_iso_pic(vertices_mesh);
    // if (Basic::save_BC(save_path + "_pic_lambda_111", vertices_mesh))
    // {
    //     std::cout << "pic-iso calculation successfull" << std::endl;
    // }
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "pic-iso-times: " << elapsed.count() << " s" << std::endl;

    // Save the obj of the cross-section
    // Basic::saveOBJ_contour(save_path + "duck_47.obj", vertices_mesh, faces_mesh, 47);
    // Basic::saveOBJ_contour(save_path + "duck_57.obj", vertices_mesh, faces_mesh, 57);
    // Basic::saveOBJ_contour(save_path + "duck_64.obj", vertices_mesh, faces_mesh, 64);
    // Basic::saveOBJ_contour(save_path + "duck_ip_111.obj", vertices_mesh, faces_mesh, vertices_cage.size());
   

}
