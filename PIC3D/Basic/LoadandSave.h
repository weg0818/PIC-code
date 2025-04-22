/*** 
 * @Author: Mao jiaping
 * @Date: 2024-08-22 10:47:57
 * @LastEditTime: 2025-02-21 11:47:11
 * @Description: 
 * @FileName: 
 * @FilePath: /Codes/C++/PIC3D/Basic/LoadandSave.h
 */

#ifndef GBC_LOADANDSAVE_H
#define GBC_LOADANDSAVE_H
#include "VertexR3.hpp"
#include "Face.hpp"

#include <vector>
#include <Eigen/Dense>

namespace Basic
{
    // Read obj file
    void loadOBJ(const std::string& path, std::vector<VertexR3>& vertices, std::vector<Face>& face);

    // Save obj file
    bool saveOBJ(const std::string& path, std::vector<VertexR3> & vertex, std::vector<Face> & faces);

    // Used to save sections, with texture coordinates.
    bool saveOBJ_contour(const std::string& path, std::vector<VertexR3> & vertex, std::vector<Face> & faces, int index = 0);

    // Save centre of gravity coordinates
    bool save_BC(const std::string& path, std::vector<VertexR3> & x);

    // Save points
    bool save_points(const std::string& path, std::vector<VertexR3> & x);

    // Save cells
    bool save_cells(const std::string& path, Eigen::MatrixXi & T);

    
}

#endif // GBC_LOADANDSAVE_H