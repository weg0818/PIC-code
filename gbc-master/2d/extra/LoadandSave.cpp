#include <sstream>
#include "LoadandSave.h"

void LoadandSave::saveVecVecDouble(const std::string & path, std::vector<std::vector<double>> & vecvecD)
{
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return;
    }
    // write data
    for (auto& item : vecvecD) {
        for (auto& bc_i : item)
        {
            outputFile << bc_i << ' ';
        }
        outputFile << std::endl;
    }
    // Close file
    outputFile.close();
    return;
}


bool LoadandSave::load_x_txt(const std::string& path, std::vector<VertexR2>& vertices)
{
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return false;
    }
    double x, y;

    while (ifs >> x >> y)
    {
        vertices.push_back(VertexR2(x, y));
    }
    ifs.close();
    return true;
}

bool LoadandSave::load_cd_txt(const std::string& path, std::vector<std::vector<double>> & cd)
{
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return false;
    }
    std::string line;
    std::size_t m = 0;

    while (getline(ifs, line))
    {
        std::istringstream iss(line);
        std::vector<double> row;
        double value;

        // Read each row of data one by one and store it in a vector
        while (iss >> value) {
            row.push_back(value);
        }

        // If the current row is not empty and the number of columns m has not been determined, update the value of m.
        if (!row.empty()) {
            if (m == 0) {
                m = row.size();  // The number of elements in the first row is the number of columns m
            } else if (row.size() != m) {
                std::cerr << "Inconsistent number of columns per row of data, file format error!" << std::endl;
                return false;
            }
            cd.push_back(row);  //  Add the row to the data
        }
    }
    ifs.close();
    return true;
}

bool LoadandSave::readData(const std::string& path, std::vector<VertexR2> & poly) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return false;
    }
    int rows, cols;
    ifs >> rows >> cols;
    // Resize the matrix to the appropriate size
    double x, y;
    for (int i = 0; i < rows; ++i) {
        ifs >> x >> y;
        poly.push_back(VertexR2(x, y));
    }
    ifs.close();
    return true;
}

bool LoadandSave::readData(const std::string& path, std::vector<Face> & faces) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return false;
    }
    int rows, cols;
    ifs >> rows >> cols;
    // Resize the matrix to the appropriate size
    double v1, v2, v3;
    for (int i = 0; i < rows; ++i) {
        ifs >> v1 >> v2 >> v3;
        faces.push_back(Face(v1, v2, v3));
    }
    ifs.close();
    return true;
}

bool LoadandSave::readData(const std::string& path, std::vector<bool> & poly) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return false;
    }
    int rows;
    ifs >> rows;
    // Resize the matrix to the appropriate size
    int x;
    for (int i = 0; i < rows; ++i) {
        ifs >> x;
        poly.push_back(x);
    }
    ifs.close();
    return true;
}

void LoadandSave::load_imformation(const std::string path, std::vector<VertexR2> & poly, 
        std::vector<VertexR2> & x, std::vector<bool> & vertex_markers,
        std::vector<Face> & faces)      
{
    std::string path_poly = path + "polygon";
    // load polygon
    if (readData(path_poly, poly))
    {
        std::cout << "Successfully loaded polygon" << std::endl;
    }

    std::string path_x = path + "x";
    // load x
    if (readData(path_x, x))
    {
        std::cout << "Successfully loaded sample point" << std::endl;
    }

    std::string path_vertex_markers = path + "vertex_markers";
    // load vertex_markers
    if (readData(path_vertex_markers, vertex_markers))
    {
        std::cout << "Successfully loaded vertex_markers" << std::endl;
    }

    std::string path_faces = path + "triangles";
    // load faces
    if (readData(path_faces, faces))
    {
        std::cout << "Successfully loaded faces" << std::endl;
    }
} 

void LoadandSave::load_imformation(const std::string path, std::vector<VertexR2> & poly)      
{
    std::string path_poly = path;
    // load polygon
    if (readData(path_poly, poly))
    {
        std::cout << "Successfully loaded polygon" << std::endl;
    }
} 

bool LoadandSave::save_BC(const std::string& path, std::vector<VertexR2> & x) {
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // write data
    // std::cout << x[4800].b()[0] << std::endl;
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

bool LoadandSave::save_imformation(const std::string& path, std::vector<Face> & faces) {
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // write data
    for (auto& item : faces) {
        outputFile << item.v[0] << ' ' << item.v[1] << ' ' << item.v[2] << std::endl;
    }
    // Close file
    outputFile.close();
    return true;
}

bool LoadandSave::save_imformation(const std::string& path, std::vector<VertexR2> & x) {
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // write data
    for (auto& item : x) {
        outputFile << item.x() << ' ' << item.y() << std::endl;
    }
    //Close file
    outputFile.close();
    return true;
}


bool LoadandSave::saveOBJ(const std::string& path, std::vector<VertexR2> & vertex, std::vector<Face> & faces)
{
    std::ofstream outputFile(path, std::ios::out);  // Create a file and open it to write data
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // write data
    for (auto& v : vertex) {
        outputFile << "v" << ' ' << v.x() << ' ' << v.y() << ' ' << 0 << std::endl;
    }

    // write data
    for (auto& f : faces) {
        outputFile << "f" << ' ' << f.v[0]+1 << ' ' << f.v[1]+1 << ' ' << f.v[2]+1 <<std::endl;
    }
    // Close file
    outputFile.close();
    return true;
}