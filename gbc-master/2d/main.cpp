#include <string>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <chrono>
// Local includes.
#include "./coords/MeanValueR2.hpp"
#include "./coords/HarmonicR2.hpp"
#include "./coords/PointwiseIterativeR2.h"
#include "./coords/IteractiveR2.h"
#include "./coords/MaximumLikelihoodR2.h"
#include "./coords/PointwiseIterativeR2_Scaled.h"
#include "./coords/PositiveMeanValueR2.h"
#include "./coords/WachspressR2.hpp"

#include "./extra/MeshR2.hpp"
#include "./extra/VertexR2.hpp"
#include "./extra/TriangulatorR2.hpp"
#include "./extra/TestCoordinatesR2.hpp"
#include "./extra/LoadandSave.h"
using namespace gbc;
using namespace LoadandSave;


// Examples.
int main(int argc, char* argv[]) {
    std::cout << "Loading data" << std::endl;
    std::string name = argv[1];
    const std::string path = "/home/mjp/Codes/C++/gbc-master/data/" + name;
    // const std::string path = "/home/mjp/Codes/C++/test-vp/test1";
    const std::string save_path = "/home/mjp/Codes/C++/gbc-master/out/" + name + "/";
    std::vector<VertexR2> poly;
    load_imformation(path, poly);

    char* end;
    const double edgeLength = std::strtod(argv[2], &end);
    // const double edgeLength = 0.08;
    // Refine the polygon to create regular mesh.
    std::vector<VertexR2> refined;
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        refined.push_back(poly[i]);
        const size_t ip = (i + 1) % n;
        const size_t numS = ceil((poly[ip] - poly[i]).length() / edgeLength);
        for (size_t j = 1; j < numS; ++j) {
            VertexR2 vert = poly[i] + (double(j) / double(numS)) * (poly[ip] - poly[i]);
            refined.push_back(vert);
        }
    }

    // Create mesh. 
    TriangulatorR2 tri(refined, edgeLength, true);
    tri.setPlanarGraph(true);

    std::vector<VertexR2> queries;
    std::vector<Face> faces;

    tri.generate(queries, faces);
    MeshR2 mesh;
    mesh.initialize(queries, faces);
    // Clean mesh from the polygon vertices.
    std::vector<VertexR2> cleaned;
    for (size_t i = 0; i < mesh.numVertices(); ++i) {
        if (mesh.vertices()[i].type == INTERIOR || mesh.vertices()[i].type == FLAT)
            cleaned.push_back(mesh.vertices()[i]);
    }
    std::cout << "Save Sampling Points" << ' ' << queries.size() << std::endl;
    // Saving the grid data
    if (save_imformation(save_path + "x", queries))
    {
        std::cout << "Save Sampling Points" << ' ' << queries.size() << std::endl;
    }
    else
    {
        std::cout << "Failed to save sampling point" << std::endl;
    }
    if (save_imformation(save_path + "faces", faces))
    {
        std::cout << "save faces" << std::endl;
    }

    if (saveOBJ(save_path + name + ".obj", queries, faces))
    {
        std::cout << "save obj" << std::endl;
    }
        
    // When calculating 'guaishou', a coordinate value that is too large will result in an illegal value.
    for (int i = 0; i < queries.size(); ++i)
    {
        queries[i] *= 0.001;
    }

    for (int i = 0; i < poly.size(); ++i)
    {
        poly[i] *= 0.001;
    }

    auto last = std::chrono::high_resolution_clock::now();
    auto cur = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;

    // Wachspress
    // last = std::chrono::high_resolution_clock::now();
    // WachspressR2 wachspress = WachspressR2(poly);
    // wachspress.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate wachspress time: " << elapsed.count() << " s" << std::endl;
    // if (save_BC(save_path + "wachspressc_lambda", queries))
    // {
    //     std::cout << "Successful calculation of wachspress" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_mvc = BarycentricPropertiesR2(poly, 1e-8, true);
    // test_mvc.check_linearReproduction(queries);
    // test_mvc.check_partitionOfUnity(queries);


    // // MVC
    // last = std::chrono::high_resolution_clock::now();
    // MeanValueR2 mvc = MeanValueR2(poly);
    // mvc.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate MVC time: " << elapsed.count() << " s" << std::endl;
    // if (save_BC(save_path + "mvc_lambda", queries))
    // {
    //     std::cout << "Successful calculation of MVC" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_mvc = BarycentricPropertiesR2(poly, 1e-8, true);
    // test_mvc.check_linearReproduction(queries);
    // test_mvc.check_partitionOfUnity(queries);

    // PMVC
    // last = std::chrono::high_resolution_clock::now();
    // PositiveMeanValueR2 pmvc = PositiveMeanValueR2(poly);
    
    // pmvc.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate PMVC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "pmvc_lambda", queries))
    // {
    //     std::cout << "Successful calculation of PMVC" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_mvc = BarycentricPropertiesR2(poly, 1e-8, true);
    // test_mvc.check_linearReproduction(queries);
    // test_mvc.check_partitionOfUnity(queries);

    // MEC
    // last = std::chrono::high_resolution_clock::now();
    // MaximumEntropyR2 mec = MaximumEntropyR2(poly);
    // mec.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate MEC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "mec_lambda", queries))
    // {
    //     std::cout << "Successful calculation of MEC" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_mec = BarycentricPropertiesR2(poly, 1e-8, true);
    // test_mec.check_linearReproduction(queries);
    // test_mec.check_partitionOfUnity(queries);

    // MLC
    // last = std::chrono::high_resolution_clock::now();
    // MaximumLikelihoodR2 mlc = MaximumLikelihoodR2(poly);
    // mlc.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate MLC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "mlc_lambda", queries))
    // {
    //     std::cout << "Successful calculation of MLC" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_mlc = BarycentricPropertiesR2(poly, 1e-8, true);
    // test_mlc.check_linearReproduction(queries);
    // test_mlc.check_partitionOfUnity(queries);

    // HC
    // last = std::chrono::high_resolution_clock::now();
    // HarmonicR2 hmc = HarmonicR2(poly);
    // hmc.setMesh(queries, faces);
    // hmc.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate HC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "hmc_lambda", queries))
    // {
    //     std::cout << "Successful calculation of HC" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_hmc = BarycentricPropertiesR2(poly, 1e-8, true);
    // test_hmc.check_linearReproduction(queries);
    // test_hmc.check_partitionOfUnity(queries);

    // LBC
    // last = std::chrono::high_resolution_clock::now();
    // LocalR2 lbc = LocalR2(poly);
    // lbc.setEdgeLength(edgeLength);
    // lbc.setMesh(queries, faces);
    // lbc.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate LBC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "lbc_lambda", queries))
    // {
    //     std::cout << "Successful calculation of LBC" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_lbc = BarycentricPropertiesR2(poly, 1e-8, true);
    // test_lbc.check_linearReproduction(queries);
    // test_lbc.check_partitionOfUnity(queries);

    // IC
    // last = std::chrono::high_resolution_clock::now();
    // IteractiveR2 ic = IteractiveR2(poly);
    // ic.compute_iteractive_times(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate IC time: " << elapsed.count() << "s" << std::endl;
    // ic.set_k(11);
    // last = std::chrono::high_resolution_clock::now();
    // ic.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate IC time: " << elapsed.count() << " 秒" << std::endl;
    // if (save_BC(save_path + "ic_lambda_11", queries))
    // {
    //     std::cout << "Successful calculation of IC" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_ic = BarycentricPropertiesR2(poly, 1e-8, true);
    // test_ic.check_linearReproduction(queries);
    // test_ic.check_partitionOfUnity(queries);
    

    // PIC
    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR2 pic = PointwiseIterativeR2(poly);
    // pic.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate PIC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "pic_lambda", queries))
    // {
    //     std::cout << "Successful calculation of PIC" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_pic = BarycentricPropertiesR2(poly, 1e-6, true);
    // test_pic.check_linearReproduction(queries);
    // test_pic.check_partitionOfUnity(queries);

    // // PIC 对比1
    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR2 pic_1 = PointwiseIterativeR2(poly);
    // pic_1.compute_1(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate PIC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "pic_lambda_1", queries))
    // {
    //     std::cout << "Successful calculation of PIC 对比1" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_pic_1 = BarycentricPropertiesR2(poly, 1e-6, true);
    // test_pic_1.check_linearReproduction(queries);
    // test_pic_1.check_partitionOfUnity(queries);

    // // PIC 对比2
    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR2 pic_1toN = PointwiseIterativeR2(poly);
    // pic_1toN.compute_1toN(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate PIC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "pic_lambda_1toN", queries))
    // {
    //     std::cout << "Successful calculation of PIC 对比2" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_pic_1toN = BarycentricPropertiesR2(poly, 1e-6, true);
    // test_pic_1toN.check_linearReproduction(queries);
    // test_pic_1toN.check_partitionOfUnity(queries);

    // // // PIC 对比3
    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR2 pic_Nto2N = PointwiseIterativeR2(poly);
    // pic_Nto2N.compute_Nto2N(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate PIC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "pic_lambda_Nto2N", queries))
    // {
    //     std::cout << "Successful calculation of PIC 对比3" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_pic_Nto2N = BarycentricPropertiesR2(poly, 1e-6, true);
    // test_pic_Nto2N.check_linearReproduction(queries);
    // test_pic_Nto2N.check_partitionOfUnity(queries);

    // // PIC 对比4
    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR2 pic_1to2N = PointwiseIterativeR2(poly);
    // pic_1to2N.compute_1to2N(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate PIC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "pic_lambda_1to2N", queries))
    // {
    //     std::cout << "Successful calculation of PIC 对比4" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_pic_1to2N = BarycentricPropertiesR2(poly, 1e-6, true);
    // test_pic_1to2N.check_linearReproduction(queries);
    // test_pic_1to2N.check_partitionOfUnity(queries);

    // PIC 对比5
    last = std::chrono::high_resolution_clock::now();
    PointwiseIterativeR2 pic_p2 = PointwiseIterativeR2(poly);
    pic_p2.compute_p2(queries);
    cur = std::chrono::high_resolution_clock::now();
    elapsed = cur - last;
    std::cout << "Calculate PIC time: " << elapsed.count() << "s" << std::endl;
    if (save_BC(save_path + "pic_lambda_p2", queries))
    {
        std::cout << "Successful calculation of PIC 对比5" << std::endl;
    }
    gbc::BarycentricPropertiesR2 test_pic_p2 = BarycentricPropertiesR2(poly, 1e-6, true);
    test_pic_p2.check_linearReproduction(queries);
    test_pic_p2.check_partitionOfUnity(queries);
    
    // PIC 对比6 1/4 1/2 1/4
    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR2 pic_p1p2 = PointwiseIterativeR2(poly);
    // pic_p1p2.compute_p1p2(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate PIC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "pic_lambda_p1p2", queries))
    // {
    //     std::cout << "Successful calculation of PIC 对比6" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_pic_p1p2 = BarycentricPropertiesR2(poly, 1e-6, true);
    // test_pic_p1p2.check_linearReproduction(queries);
    // test_pic_p1p2.check_partitionOfUnity(queries);

    // Isolated PIC
    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR2 iso_pic = PointwiseIterativeR2(poly);
    // // std::vector<std::vector<VertexR2>> iso_p(2);
    // // iso_p[0].push_back(VertexR2(1.5, 1.5));
    // // iso_p[1].push_back(VertexR2(-1.5, 1.5));
    // //non-manifold
    // std::vector<VertexR2> iso_p;
    // iso_p.push_back(VertexR2(0.1, 0.9));
    // iso_p.push_back(VertexR2(0.3, 0.5));
    // iso_p.push_back(VertexR2(0.7, 0.5));
    // iso_p.push_back(VertexR2(0.9, 0.9));
    // iso_p.push_back(VertexR2(0.1, 0.1));
    // iso_p.push_back(VertexR2(0.9, 0.1));
    // std::vector<std::vector<int>> iso_p_line;
    // iso_p_line.push_back(std::vector({0, 1}));
    // iso_p_line.push_back(std::vector({1, 2}));
    // iso_p_line.push_back(std::vector({2, 3}));
    // iso_p_line.push_back(std::vector({1, 4}));
    // iso_p_line.push_back(std::vector({2, 5}));
    // iso_pic.set_iso_p(iso_p, iso_p_line);
    // iso_pic.compute_iso_pic(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate iso_PIC time:" << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "iso_pic_lambda", queries))
    // {
    //     std::cout << "Successful calculation of iso_PIC" << std::endl;
    // }


    // // SPIC 
    // load cd1 cd2
    // std::vector<VertexR2> queries;
    // std::vector<std::vector<double>> cd1;
    // std::vector<std::vector<double>> cd2;
    // // LoadandSave::load_x_txt("/home/mjp/Codes/C++/data/M_cd/x", queries);
    // // LoadandSave::load_cd_txt("/home/mjp/Codes/C++/gbc-master/out/M/cd1", cd1);
    // // LoadandSave::load_cd_txt("/home/mjp/Codes/C++/gbc-master/out/M/cd2", cd2);
    // last = std::chrono::high_resolution_clock::now();
    // ScaledPointwiseIterativeR2 spic = ScaledPointwiseIterativeR2(poly);
    // spic.set_f(faces);
    // spic.compute_scald_cd(cd1, cd2, queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate the distance. " << elapsed.count() << "s" << std::endl;
    // spic.set_scaled_cd(cd1, cd2);
    // spic.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // elapsed = cur - last;
    // std::cout << "Calculate SPIC time: " << elapsed.count() << "s" << std::endl;
    // if (save_BC(save_path + "spic_lambda", queries))
    // {
    //     std::cout << "Successful calculation of SPIC" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test_pic = BarycentricPropertiesR2(poly, 1e-8, true);
    // test_pic.check_linearReproduction(queries);
    // test_pic.check_partitionOfUnity(queries);

    // //save cd1，cd2
    // LoadandSave::saveVecVecDouble("/home/mjp/Codes/C++/gbc-master/out/M/cd1", cd1);
    // LoadandSave::saveVecVecDouble("/home/mjp/Codes/C++/gbc-master/out/M/cd2", cd2);
}
