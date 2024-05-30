#include "vec3.cuh"
#include <iostream>

__host__ __device__ vec3::vec3() : e{ 0, 0, 0 } {}

__host__ __device__ vec3::vec3(double e0, double e1, double e2) : e{ e0, e1, e2 } {}

__host__ __device__ double vec3::x() const { return e[0]; }
__host__ __device__ double vec3::y() const { return e[1]; }
__host__ __device__ double vec3::z() const { return e[2]; }

__host__ __device__ vec3 vec3::operator-() const { return vec3(-e[0], -e[1], -e[2]); }
__host__ __device__ double vec3::operator[](int i) const { return e[i]; }
__host__ __device__ double& vec3::operator[](int i) { return e[i]; }

__host__ __device__ vec3& vec3::operator+=(const vec3& v) {
    e[0] += v.e[0];
    e[1] += v.e[1];
    e[2] += v.e[2];
    return *this;
}

__host__ __device__ vec3& vec3::operator*=(double t) {
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;
    return *this;
}

__host__ __device__ vec3& vec3::operator/=(double t) {
    return *this *= 1 / t;
}

__host__ __device__ double vec3::length() const {
    return sqrt(length_squared());
}

__host__ __device__ double vec3::length_squared() const {
    return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
}

__host__ __device__ std::ostream& operator<<(std::ostream& out, const vec3& v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

__host__ __device__ inline vec3 operator+(const vec3& u, const vec3& v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

__host__ __device__ inline vec3 operator-(const vec3& u, const vec3& v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

__host__ __device__ inline vec3 operator*(const vec3& u, const vec3& v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

__host__ __device__ inline vec3 operator*(double t, const vec3& v) {
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

__host__ __device__ inline vec3 operator*(const vec3& v, double t) {
    return t * v;
}

__host__ __device__ inline vec3 operator/(const vec3& v, double t) {
    return (1 / t) * v;
}

__host__ __device__ inline double dot(const vec3& u, const vec3& v) {
    return u.e[0] * v.e[0]
        + u.e[1] * v.e[1]
        + u.e[2] * v.e[2];
}

__host__ __device__ inline vec3 cross(const vec3& u, const vec3& v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
        u.e[2] * v.e[0] - u.e[0] * v.e[2],
        u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

__host__ __device__ inline vec3 unit_vector(const vec3& v) {
    return v / v.length();
}

__host__ __device__ inline vec3 crossProduct_(const vec3& a, const vec3& b) {
    return vec3(a.y() * b.z() - a.z() * b.y(), a.z() * b.x() - a.x() * b.z(), a.x() * b.y() - a.y() * b.x());
}

__host__ __device__ inline double dotProduct_(const vec3& a, const vec3& b) {
    return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

__host__ __device__ Plane::Plane() : A(0), B(0), C(0), D(0) {}

__host__ __device__ Plane::Plane(double _A, double _B, double _C, double _D) : A(_A), B(_B), C(_C), D(_D) {}

__host__ __device__ Plane::Plane(const vec3& normal, const point3& point_on_face) {
    A = normal.x();
    B = normal.y();
    C = normal.z();
    D = -(normal.x() * point_on_face.x() + normal.y() * point_on_face.y() + normal.z() * point_on_face.z());
}

__host__ __device__ Plane::Plane(const Plane& other) : A(other.A), B(other.B), C(other.C), D(other.D) {}

__host__ __device__ Plane::Plane(Plane&& other) noexcept : A(std::exchange(other.A, 0)),
B(std::exchange(other.B, 0)),
C(std::exchange(other.C, 0)),
D(std::exchange(other.D, 0)) {}

__host__ __device__ Plane& Plane::operator=(const Plane& other) {
    if (this != &other) {
        A = other.A;
        B = other.B;
        C = other.C;
        D = other.D;
    }
    return *this;
}

__host__ __device__ std::ostream& operator<<(std::ostream& out, const Plane& v) {
    if (abs(v.A) != 1 && abs(v.A) != 0) { out << v.A; }
    if (abs(v.A) != 0 && v.A > 0) { out << "x"; }
    if (abs(v.A) != 0 && v.A < 0) { out << "-x"; }
    if (abs(v.A) != 0 && abs(v.B) != 0 && v.B > 0) { out << "+"; }
    if (abs(v.B) != 0 && v.B < 0) { out << ""; }
    if (abs(v.B) != 1 && abs(v.B) != 0) { out << v.B; }
    if (abs(v.B) != 0 && v.B > 0) { out << "y"; }
    if (abs(v.B) != 0 && v.B < 0) { out << "-y"; }
    if (abs(v.A) != 0 && abs(v.B) != 0 && abs(v.C) != 0 && v.C > 0) { out << "+"; }
    if (abs(v.C) != 0 && v.C < 0) { out << ""; }
    if (abs(v.C) != 1 && abs(v.C) != 0) { out << v.C; }
    if (abs(v.C) != 0 && v.C > 0) { out << "z"; }
    if (abs(v.C) != 0 && v.C < 0) { out << "-z"; }
    if (abs(v.D) != 0 && v.D > 0) { out << "+"; }
    if (abs(v.D) != 0 && v.D < 0) { out << ""; }
    if (abs(v.D) != 0) { out << v.D << " = 0"; }
    return out;
}

__host__ __device__ Vector::Vector(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
__host__ __device__ Vector::Vector() : x(0.0f), y(0.0f), z(0.0f) {}

__global__ void MyKernel(Vector* vectors, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < size) {
        vectors[idx] = Vector(1.0f, 2.0f, 3.0f); 
    }
    vec3 x = vec3();
}
__host__ void printVectors(Vector* vectors, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << "x: " << vectors[i].x << ", y: " << vectors[i].y << ", z: " << vectors[i].z << std::endl;
    }
}