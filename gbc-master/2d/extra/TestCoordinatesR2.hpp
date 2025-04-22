// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    This is a test class for testing generalized barycentric coordinates.

*/

#ifndef GBC_TESTCOORDINATESR2_HPP
#define GBC_TESTCOORDINATESR2_HPP

// STL includes.
#include <vector>
#include <iostream>
#include <cmath>

// Local includes.
#include "Face.hpp"
#include "VertexR2.hpp"
#include "TriangulatorR2.hpp"
#include "BarycentricCoordinatesR2.hpp"

#include "../coords/MeanValueR2.hpp"
#include "../coords/MaximumEntropyR2.hpp"
#include "../coords/HarmonicR2.hpp"
#include "../coords/LocalR2.hpp"

namespace gbc {

    // Different theoretical properties of generalized barycentric coordinates.
    class BarycentricPropertiesR2 {

    public:
        // Constructor.
        BarycentricPropertiesR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10, const bool verbose = true)
                : _v(v), _tol(tol), _verbose(verbose) { }

        // Properties.
        bool partitionOfUnity(const std::vector<double> &b) const {

            const int n = _v.size();

            double sum = 0.0;
            for (int i = 0; i < n; ++i) sum += b[i];

            const double diff = fabs(sum - 1.0);

            if (diff > _tol) {
                if (_verbose)
                    std::cerr << "ERROR: Partition of unity difference is " << diff << ";" << std::endl;
                return false;
            }

            // if (_verbose) std::cout << "Partition of unity: SUCCESS" << std::endl;
            return true;
        }

        bool linearReproduction(const VertexR2 &p, const std::vector<double> &b) const {

            const int n = _v.size();

            VertexR2 sum;
            for (int i = 0; i < n; ++i) sum += b[i] * _v[i];

            const VertexR2 diff = sum - p;

            if (fabs(diff.x()) > _tol || fabs(diff.y()) > _tol) {
                if (_verbose)
                    std::cerr << "ERROR: Linear reproduction difference is " << diff << ";" << std::endl;
                return false;
            }

            // if (_verbose) std::cout << "Linear reproduction: SUCCESS" << std::endl;
            return true;
        }

        void check_partitionOfUnity(const std::vector<VertexR2> &x)
        {
            for (int i = 0; i < x.size(); ++i)
            {
                if (!partitionOfUnity(x[i].b()))
                {
                    std::cout << "ERROR: Partition of unity difference is " << i << ";" <<  std::endl;
                    // return;
                }
            }
            std::cout << "SUCCESS: Partition of unity" << std::endl;
        }

        void check_linearReproduction(const std::vector<VertexR2> &x)
        {
            for (int i = 0; i < x.size(); ++i)
            {
                if (!linearReproduction(x[i], x[i].b()))
                {
                    std::cout << "ERROR: Linear reproduction difference is " << i << ";" <<  std::endl;
                    // return;
                }
            }
            std::cout << "SUCCESS: Linear reproduction" << std::endl;
        }

        bool boundedness(const std::vector<double> &b) const {

            const int n = _v.size();

            for (int i = 0; i < n; ++i) {
                if (b[i] < 0.0 || b[i] > 1.0) {
                    if (_verbose)
                        std::cerr << "ERROR: Value out of range [0,1] with index " << i << ": " << b[i] << ";" << std::endl;
                    return false;
                }
            }

            if (_verbose) std::cout << "Boundedness: SUCCESS" << std::endl;
            return true;
        }


    private:
        // Vertices of the polygon.
        const std::vector<VertexR2> &_v;

        // Tolerance.
        const double _tol;

        // Print data.
        const bool _verbose;
    };

    // Class that verifies different properties and behaviour of generalzied barycentric coordinates in R2.
    class TestCoordinatesR2 {

    public:
        // Constructor.
        TestCoordinatesR2() : _tol(1.0e-10) { }
    private:
        // Tolerance.
        const double _tol;

        // Create triangle mesh for a given polygon.
        void createMesh(const std::vector<VertexR2> &poly,
                        std::vector<VertexR2> &tp, 
                        std::vector<Face> &tf, 
                        const double edgeLength) {

            // Refine the polygon to create regular mesh.
            std::vector<VertexR2> refined;
            const int n = poly.size();

            for (int i = 0; i < n; ++i) {
                refined.push_back(poly[i]);

                const int ip = (i + 1) % n;
                const int numS = ceil((poly[ip] - poly[i]).length() / edgeLength);

                for (int j = 1; j < numS; ++j) {

                    VertexR2 vert = poly[i] + (double(j) / double(numS)) * (poly[ip] - poly[i]);
                    refined.push_back(vert);
                }
            }

            // Create mesh.
            TriangulatorR2 tri(refined, edgeLength, true);
            tri.setPlanarGraph(true);

            tp.clear();
            tf.clear();

            tri.generate(tp, tf);
        }

    };

} // namespace gbc

#endif // GBC_TESTCOORDINATESR2_HPP
