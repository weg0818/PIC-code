#ifndef VERTEXR3_HPP
#define VERTEXR3_HPP
#include <cmath>
#include <vector>
#include <iostream>
#include <cassert>
class VertexR3
{
    private:
        double _x;
        double _y;
        double _z;

        double _tol;
        std::vector<double> _b;
    public:

        // Constructors.
        VertexR3() :  _x(0.0), _y(0.0), _z(0.0), _b(), _tol(1.0e-10) { }
        
        VertexR3(const std::size_t bSize) :  _x(0.0), _y(0.0), _z(0.0), _b(bSize, 0.0), _tol(1.0e-10) { }
    
        VertexR3(const VertexR3 &v) : _x(v.x()), _y(v.y()), _z(v.z()), _b(v.b()), _tol(1.0e-10) { }

        VertexR3(const double x, const double y, const double z) :  _x(x), _y(y), _z(z), _b(), _tol(1.0e-10) { }

        // Return size of the barycentric coordinate vector.
        size_t size() const {
            return _b.size();
        }

        // Return x coordinate.
        inline double &x() {
            return _x;
        }

        // Return const x coordinate.
        inline double x() const {
            return _x;
        }

        // Return y coordinate.
        inline double &y() {
            return _y;
        }

        // Return const y coordinate.
        inline double y() const {
            return _y;
        }

        // Return z coordinate.
        inline double &z() {
            return _z;
        }

        // Return const z coordinate.
        inline double z() const {
            return _z;
        }

        // Return barycentric coordinates attached to the vertex.
        inline std::vector<double> &b() {
            return _b;
        }

        // Return const barycentric coordinates attached to the vertex.
        inline const std::vector<double> &b() const {
            return _b;
        }

        // Return barycentric coordinate value.
        inline double b(const int i) const {
            assert(i >= 0 && i < (int) _b.size());
            return _b[i];
        }

        // Overload some basic operators.

        // Square brakets operator by value.
        inline double operator[](const int i) const {
            
            assert(i >= 0 && i < 3);
            return (&_x)[i];
        }

        // Square brakets operator by reference.
        inline double &operator[](const int i) {
            
            assert(i >= 0 && i < 3);
            return (&_x)[i];
        }

        // Addition of two vertices without creating a new vertex.
        void operator+=(const VertexR3 &v) {
            const std::size_t bSize = _b.size();

            _x += v.x();
            _y += v.y();
            _z += v.z();

            for (size_t i = 0; i < bSize; ++i) _b[i] += v.b()[i];
        }

        // Subtraction of two vertices without creating a new vertex.
        void operator-=(const VertexR3 &v) {
            const std::size_t bSize = _b.size();

            _x -= v.x();
            _y -= v.y();
            _z -= v.z();

            for (size_t i = 0; i < bSize; ++i) _b[i] -= v.b()[i];
        }

        // Multiplication by a constant from the right without creating a new vertex.
        void operator*=(const double scalar) {
            const std::size_t bSize = _b.size();

            _x *= scalar;
            _y *= scalar;
            _z *= scalar;

            for (size_t i = 0; i < bSize; ++i) _b[i] *= scalar;
        }

        // Equal equal operator.
        inline bool operator==(const VertexR3 &v) const {
            return fabs(_x - v.x()) < _tol && fabs(_y - v.y()) < _tol && fabs(_z - v.z()) < _tol;
        }

        // Not equal operator.
        inline bool operator!=(const VertexR3 &v) const {
            return !(this->operator==(v));
        }

        // Squared length.
        inline double squaredLength() const {
            return x() * x() + y() * y() + z() * z();
        }

        // Length - Euclidean 2-norm.
        inline double length() const {
            return sqrt(squaredLength());
        }

         // Addition of two vertices creating a new vertex.
        VertexR3 operator+(const VertexR3 &v) const {
            return VertexR3(_x + v.x(), _y + v.y(), _z + v.z());
        }

         // Subtraction of two vertices creating a new vertex.
        VertexR3 operator-(const VertexR3 &v) const {
            return VertexR3(_x - v.x(), _y - v.y(), _z - v.z());
        }

         // Multiplication of two vertices creating a new vertex.
        VertexR3 operator*(const double scalar) const {
            return VertexR3(_x * scalar, _y * scalar, _z * scalar);
        }

        // Division of a vertex by a scalar
        VertexR3 operator/(double scalar) const {
            assert(scalar != 0 && "error : scalar is 0");
            return VertexR3(_x / scalar, _y / scalar, _z / scalar);
        }

        // Optional: Compound assignment operator for scalar division
        VertexR3& operator/=(double scalar) {
            assert(scalar != 0 && "error : scalar is 0");
            _x /= scalar;
            _y /= scalar;
            _z /= scalar;
            return *this;
        }
};



#endif // VERTEXR3_HPP