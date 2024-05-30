#include <map>
#include <vector>
#include "vec3.cuh"
#include "Ray.cuh"


using namespace std;
double positiveInfinity = INFINITY;
double negativeInfinity = -INFINITY;
__host__ __device__ ray::ray() {};
__host__ __device__ ray::ray(const point3& origin, const vec3& direction) : orig(origin), dir(direction) {};
__host__ __device__ point3 ray:: at(double t) const { return this->dir*t + this->orig; }
__host__ __device__ point3 ray::findIntersection(const Plane& plane) const
{
double t = -(plane.A * this->orig.x() + plane.B * this->orig.y() + plane.C * this->orig.z() + plane.D) /
    (plane.A * this->dir.x() + plane.B * this->dir.y() + plane.C * this->dir.z());
    point3 intersection(at(t));
    return intersection;
}

bool Face_hit(float* pl, const ray& r, int polygon_langht ,int* polygon, float** vertices_coords)
{
    point3 intersection = r.findIntersection(pl);
    if (intersection[0] == INFINITY || intersection[0] == -INFINITY || intersection[1] == INFINITY ||
        intersection[1] == -INFINITY == -INFINITY || intersection[2] == INFINITY || intersection[2] == -INFINITY) {
        return false;
    }
    vec3 edge;
    for (size_t i = 0; i < polygon_langht; ++i) {
        size_t next_index = (i + 1) % polygon_langht;
        edge = vertices_coords[polygon[next_index]] - vertices_coords[polygon[i]];
        vec3 vp = intersection - vertices_coords[polygon[i]];
        vec3 n = crossProduct_(edge, vp);
        vec3 normal = vec3(pl.A, pl.B, pl.C);
        if (dotProduct_(n, normal) < 0) {
            return false;
        }
    }

////    return true;

};