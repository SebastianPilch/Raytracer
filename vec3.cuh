#ifndef VEC3_CUH
#define VEC3_CUH

#include <iostream>
#include <fstream>
#include <string>

using std::sqrt;

__host__ __device__ class vec3 {
public:
    double e[3];

    __host__ __device__ vec3();
    __host__ __device__ vec3(double e0, double e1, double e2);

    __host__ __device__ double x() const;
    __host__ __device__ double y() const;
    __host__ __device__ double z() const;

    __host__ __device__ vec3 operator-() const;
    __host__ __device__ double operator[](int i) const;
    __host__ __device__ double& operator[](int i);

    __host__ __device__ vec3& operator+=(const vec3& v);
    __host__ __device__ vec3& operator*=(double t);
    __host__ __device__  vec3& operator/=(double t);

    __host__ __device__ double length() const;
    __host__ __device__ double length_squared() const;
};

using point3 = vec3;
using color = vec3;

__host__ __device__ std::ostream& operator<<(std::ostream& out, const vec3& v);
__host__ __device__   vec3 operator+(const vec3& u, const vec3& v);
__host__ __device__   vec3 operator-(const vec3& u, const vec3& v);
__host__ __device__   vec3 operator*(const vec3& u, const vec3& v);
__host__ __device__   vec3 operator*(double t, const vec3& v);
__host__ __device__   vec3 operator*(const vec3& v, double t);
__host__ __device__   vec3 operator/(const vec3& v, double t);
__host__ __device__   double dot(const vec3& u, const vec3& v);
__host__ __device__   vec3 cross(const vec3& u, const vec3& v);
__host__ __device__   vec3 unit_vector(const vec3& v);
__host__ __device__   vec3 crossProduct_(const vec3& a, const vec3& b);
__host__ __device__   double dotProduct_(const vec3& a, const vec3& b);

struct Plane {
    float A, B, C, D;

    __host__ __device__ Plane();
    __host__ __device__ Plane(double _A, double _B, double _C, double _D);
    __host__ __device__ Plane(const vec3& normal, const point3& point_on_face);
    __host__ __device__ Plane(const Plane& other);
    __host__ __device__ Plane(Plane&& other) noexcept;
    __host__ __device__ Plane& operator=(const Plane& other);
};

__host__ __device__ std::ostream& operator<<(std::ostream& out, const Plane& v);


class Vector {
public:
    float x, y, z;
    __host__ __device__ Vector(float _x, float _y, float _z);
    __host__ __device__ Vector();
};

__global__ void MyKernel(Vector* vectors, int size, vec3* z);
__host__ void printVectors(Vector* vectors, int size);

#endif