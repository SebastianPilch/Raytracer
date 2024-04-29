#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

using std::sqrt;

class vec3 {
public:
    double e[3];

    vec3() : e{ 0,0,0 } {}
    vec3(double e0, double e1, double e2) : e{ e0, e1, e2 } {}

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }

    vec3& operator+=(const vec3& v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    vec3& operator*=(double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    vec3& operator/=(double t) {
        return *this *= 1 / t;
    }

    double length() const {
        return sqrt(length_squared());
    }

    double length_squared() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

};

// point3 is just an alias for vec3, but useful for geometric clarity in the code.
using point3 = vec3;
using color = vec3;


// Vector Utility Functions

inline std::ostream& operator<<(std::ostream& out, const vec3& v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3& u, const vec3& v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3& u, const vec3& v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3& u, const vec3& v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3& v) {
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline vec3 operator*(const vec3& v, double t) {
    return t * v;
}

inline vec3 operator/(const vec3& v, double t) {
    return (1 / t) * v;
}

inline double dot(const vec3& u, const vec3& v) {
    return u.e[0] * v.e[0]
        + u.e[1] * v.e[1]
        + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3& u, const vec3& v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
        u.e[2] * v.e[0] - u.e[0] * v.e[2],
        u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(const vec3& v) {
    return v / v.length();
}

inline vec3 crossProduct_(const vec3& a, const vec3& b) {
    return vec3(a.y() * b.z() - a.z() * b.y(), a.z() * b.x() - a.x() * b.z(), a.x() * b.y() - a.y() * b.x());
}
inline double dotProduct_(const vec3& a, const vec3& b) {
    return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

struct Plane {
    double A, B, C, D;

    Plane() : A(0), B(0), C(0), D(0) {};
    Plane(double _A, double _B, double _C, double _D) : A(_A), B(_B), C(_C), D(_D) {}

    Plane(const vec3& normal, const point3& point_on_face)
    {
        A = normal.x();
        B = normal.y();
        C = normal.z();
        D = -(normal.x() * point_on_face.x() + normal.y() * point_on_face.y() + normal.z() * point_on_face.z());
    }
    Plane(const Plane& other) : A(other.A), B(other.B), C(other.C), D(other.D) {}

    Plane(Plane&& other) noexcept : A(std::exchange(other.A, 0)),
        B(std::exchange(other.B, 0)),
        C(std::exchange(other.C, 0)),
        D(std::exchange(other.D, 0)) {}

    Plane& operator=(const Plane& other) {
        if (this != &other) {
            A = other.A;
            B = other.B;
            C = other.C;
            D = other.D;
        }
        return *this;
    }

};

inline std::ostream& operator <<(std::ostream& out, const Plane& v) {
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
    if (abs(v.D) != 0){ out << v.D << " = 0"; }
    return out ;
}


#endif