#include <fstream>
#include <sstream>
#include <Eigen/Core>
#include <igl/readOBJ.h>
#include "LoadandSave.h"

void Basic::loadOBJ(const std::string& path, std::vector<VertexR3>& vertices, std::vector<Face>& face)
{
    
    std::string line;
    double x, y, z;
    int a, b, c;
    std::vector<std::vector<double>> V;
    std::vector<std::vector<int>> F;
    if (!igl::readOBJ(path, V, F))
    {
        std::cerr << "Error: Unable to open file " << path << std::endl;
        return;
    }

    for (int i = 0; i < V.size(); ++i)
    {
        vertices.push_back(VertexR3(V[i][0], V[i][1], V[i][2]));
    }
     for (int i = 0; i < F.size(); ++i)
    {
        face.push_back(Face(F[i][0], F[i][1], F[i][2]));
    }
    
}


bool Basic::save_BC(const std::string& path, std::vector<VertexR3> & x) {
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // write data
    for (auto& item : x) {
        std::vector<double> & bc = item.b(); 
        for (auto& bc_i : bc)
        {
            outputFile << bc_i << ' ';
        }
        outputFile << std::endl;
    }
    // Close file
    outputFile.close();
    return true;
}

bool Basic::saveOBJ(const std::string& path, std::vector<VertexR3> & vertex, std::vector<Face> & faces)
{
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // write data
    for (auto& v : vertex) {
        outputFile << "v" << ' ' << v.x() << ' ' << v.y() << ' ' << v.z() <<std::endl;
    }

    // write data
    for (auto& f : faces) {
        outputFile << "f" << ' ' << f.v[0]+1 << ' ' << f.v[1]+1 << ' ' << f.v[2]+1 <<std::endl;
    }
    // Close file
    outputFile.close();
    return true;
}

bool Basic::saveOBJ_contour(const std::string& path, std::vector<VertexR3> & vertex, std::vector<Face> & faces, int index)
{
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // write data
    outputFile << "mtllib contour.mtl" << std::endl;
    for (auto& v : vertex) {
        outputFile << "v" << ' ' << v.x() << ' ' << v.y() << ' ' << v.z() <<std::endl;
    }

    for (auto& v : vertex) {
        if (v.b()[index] < 0.05)
        {
             outputFile << "vt" << ' ' << 0.05 << ' ' << 0.5 << std::endl;
        }
        else
        {
            outputFile << "vt" << ' ' << v.b()[index] << ' ' << 0.5 << std::endl;
        }
    }
    outputFile << "usemtl texture1.png" << std::endl;

    // write data
    for (auto& f : faces) {
        outputFile << "f" << ' ' << f.v[0]+1 << "/" << f.v[0]+1 << ' ' << f.v[1]+1 << "/" << f.v[1]+1 << ' ' << f.v[2]+1 << "/" << f.v[2]+1<<std::endl;
    }
    // Close file
    outputFile.close();
    return true;
}


// Save points
bool Basic::save_points(const std::string& path, std::vector<VertexR3> & x)
{
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // write data
    for (auto& v : x) {
        outputFile << v.x() << ' ' << v.y() << ' ' << v.z() << std::endl;
    }
    // Close file
    outputFile.close();
    return true;
}

// Save cells
bool Basic::save_cells(const std::string& path, Eigen::MatrixXi & T)
{
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    int n = T.rows();
    // write data
    for (int i = 0; i < n; ++i) {
        outputFile << 4 << ' ' << T(i, 0) << ' ' << T(i, 1) << ' ' << T(i, 2) << ' ' << T(i, 3) << std::endl;
    }
    // Close file
    outputFile.close();
    return true;
}


