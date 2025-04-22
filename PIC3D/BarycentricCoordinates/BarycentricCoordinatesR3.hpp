
#ifndef GBC_BARYCENTRICCOORDINATESR3_HPP
#define GBC_BARYCENTRICCOORDINATESR3_HPP

// STL includes.
#include <vector>
#include <utility>
#include <omp.h>
#include <chrono>

// Local includes.
#include "../Basic/VertexR3.hpp"
#include "../Basic/Face.hpp"
#include "../Basic/GBCmath.h"

namespace gbc {

    // Generalized barycentric coordinates in R3.
    class BarycentricCoordinatesR3 {

    public:
        // Constructor.
        BarycentricCoordinatesR3(const std::vector<VertexR3> &v, const std::vector<Face> &f, const double tol) 
        : _v(v), _tol(tol), _f(f){
            
            assert(!_v.empty());
            assert(!_f.empty());
        }

        // This is a virtual function that is used to compute all coordinate classes at the same time.
        virtual void bc(std::vector<VertexR3> &p) = 0; 

        // Return name of the current coordinate function.
        virtual std::string name() const = 0;

        void bc_triangle(const VertexR3 & a, const VertexR3 & b, const VertexR3 & c, const VertexR3 & p, double & v, double & w, double & u)
        {
            VertexR3 v0 = c - a;
            VertexR3 v1 = b - a;
            VertexR3 v2 = p - a;

            double d00 = Basic::dot(v0, v0);
            double d01 = Basic::dot(v0, v1);
            double d11 = Basic::dot(v1, v1);
            double d20 = Basic::dot(v2, v0);
            double d21 = Basic::dot(v2, v1);

            double dm = d00 * d11 - d01 * d01;
            v = (d11 * d20 - d01 * d21) / dm;
            w = (d00 * d21 - d01 * d20) / dm;
            u = 1.0 - v - w;
        }

    protected:
        // Vertices of the polygon.
        const std::vector<VertexR3> &_v;
        // Faces of the polygon.
        const std::vector<Face> &_f;

        // Tolerance.
        const double _tol;

        int n_v = 0;
        int n_f = 0;

        // Compute boundary coordinates (linearly interpolate along the polygon's edges).
        bool computeBoundaryCoordinates(const VertexR3 &p, std::vector<double> &bb) {
            for (int i = 0; i < _f.size(); ++i)
            {
                int i0 = _f[i].v[0], i1 = _f[i].v[1], i2 = _f[i].v[2];
                VertexR3 a = _v[i1] - _v[i0];
                VertexR3 b = _v[i2] - _v[i0];
                VertexR3 c = p - _v[i0];
                VertexR3 acosb = Basic::cross(a, b);
                double A = acosb.length();
                double tmp = - Basic::dot(acosb, c);
                double tmp1 = tmp / ( A / 2);
                if (-(1e-8) < tmp1 && tmp1 < 1e-8)
                {
                    //  点在面上，计算三角形重心坐标
                    VertexR3 v0 = _v[i0] - p;
                    VertexR3 v1 = _v[i1] - p;
                    VertexR3 v2 = _v[i2] - p;
                    double A0 = (Basic::cross(v1, v2)).length();
                    double A1 = (Basic::cross(v2, v0)).length();
                    double A2 = (Basic::cross(v0, v1)).length();
                    
                    if (-_tol < A0 + A1 + A2 - A && A0 + A1 + A2 - A < _tol)
                    {
                        bb[i0] = A0 / A;
                        bb[i1] = A1 / A;
                        bb[i2] = A2 / A;
                        return true;
                    }   
                }
            }
            return false;
        }

        // Compute boundary coordinates at p using the internal storage from the VertexR3 class.
        inline void computeBoundaryCoordinates(VertexR3 &p) {
            computeBoundaryCoordinates(p, p.b());
        }

    };

} // namespace gbc

#endif // GBC_BARYCENTRICCOORDINATESR3_HPP
