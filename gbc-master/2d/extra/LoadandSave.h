/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-25 13:31:53
 * @LastEditTime: 2025-02-22 15:03:39
 * @Description: 
 * @FileName: 
 * @FilePath: /Codes/C++/gbc-master/2d/extra/LoadandSave.h
 */
#ifndef LOADANDSAVE_H
#define LOADANDSAVE_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "VertexR2.hpp"
#include "MeshR2.hpp"

using namespace gbc;

namespace LoadandSave
{
    void saveVecVecDouble(const std::string & path, std::vector<std::vector<double>> & vecvecD);
    
    bool load_x_txt(const std::string& path, std::vector<VertexR2>& vertices);

    bool load_cd_txt(const std::string& path, std::vector<std::vector<double>> & cd);

    bool readData(const std::string& path, std::vector<VertexR2> & poly);

    bool readData(const std::string& path, std::vector<Face> & faces);

    bool readData(const std::string& path, std::vector<bool> & poly);

    bool readData(const std::string& path, std::vector<bool> & poly);

    void load_imformation(const std::string path, std::vector<VertexR2> & poly, 
        std::vector<VertexR2> & x, std::vector<bool> & vertex_markers,
        std::vector<Face> & faces);

    void load_imformation(const std::string path, std::vector<VertexR2> & poly);


    bool save_BC(const std::string& path, std::vector<VertexR2> & x);

    bool save_imformation(const std::string& path, std::vector<Face> & faces);

    bool save_imformation(const std::string& path, std::vector<VertexR2> & x);

    bool saveOBJ(const std::string& path, std::vector<VertexR2> & vertex, std::vector<Face> & faces);

}

#endif // LOADANDSAVE_H