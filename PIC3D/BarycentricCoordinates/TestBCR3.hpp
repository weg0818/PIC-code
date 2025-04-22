
#ifndef GBC_TESTBC_HPP
#define GBC_TESTBC_HPP

// STL includes.
#include <vector>

// Local includes.
#include "../Basic/VertexR3.hpp"
#include "../Basic/Face.hpp"
#include "PointwiseIterativeR3.h"

namespace gbc {

    // Generalized barycentric coordinates in R3.
    class TestBCR3 {

    public:
        TestBCR3()
        {

        }
        TestBCR3(std::vector<VertexR3> & vertices, std::vector<VertexR3> & x, double tol = 1e-9) : _tol(tol)
        {
            std::cout << "Start Test" << std::endl;
            test_nonegativity(x);
            test_partitionOfUnity(x);
            test_linearReproduction(vertices, x);
        }

    private:

        void test_nonegativity(std::vector<VertexR3> & x)
        {
            int out_i = 0;
            for (auto & x_i : x)
            {
                for (auto & bc : x_i.b())
                {
                    if (bc < 0)
                    {
                        std::cout << "ERROR: Not satisfied with no-negativity:" << out_i << std::endl;
                        return;
                    }
                }
                out_i++;
            }
            std::cout << "SUCCESS: satisfied with no-negativity" << std::endl;    
        }

        void test_partitionOfUnity(std::vector<VertexR3> & x)
        {
            double sum = 0.0;
            int i = 0;
            for (auto & x_i : x)
            {
                sum = 0.0;
                for (auto & bc : x_i.b())
                {
                    sum += bc;
                }
                if (fabs(sum - 1.0) > _tol)
                {
                    std::cout << "ERROR: Not satisfied with Partition of unity:" << i << std::endl;
                    return;
                }
                i++;
            }
            std::cout << "SUCCESS: Partition of unity" << std::endl;    
        }

        void test_linearReproduction(std::vector<VertexR3> & vertices, std::vector<VertexR3> & x)
        {
            const size_t n = vertices.size();
            int out_i = 0;
            for (auto & x_i : x)
            {
                VertexR3 sum(0.0, 0.0, 0.0);
                std::vector<double> bc = x_i.b();
                for (int i = 0; i < n; ++i)
                {
                    sum += vertices[i] * bc[i];
                }
                const VertexR3 diff = sum - x_i;
                if (fabs(diff.x()) > _tol || fabs(diff.y()) > _tol || abs(diff.z()) > _tol) {
                    std::cout << "ERROR: Linear reproduction:" << out_i << ' '<<fabs(diff.x()) << ' '<< fabs(diff.y()) << ' '<< fabs(diff.z()) << std::endl;
                }
                out_i++;
            }
            std::cout << "SUCCESS: Linear reproduction" << std::endl;
        }

    private: 
        double _tol;
    };

} // namespace gbc

#endif // GBC_TESTBC_HPP
