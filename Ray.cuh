#ifndef RAY_H
#define RAY_H
#include <map>
#include <vector>
#include "vec3.cuh"
using namespace std;


struct inter_ 
{
    double t;
    point3 intersection;

    __host__ __device__ inter_() : t(0), intersection() {}
    __host__ __device__ inter_(double t, point3 intersection) : t(t), intersection(intersection) {}

};




class ray {
public:
    __host__ __device__ ray();
    __host__ __device__ ray(const point3& origin, const vec3& direction);
    __host__ __device__ const point3& origin() { return orig; };
    __host__ __device__ const vec3& direction() { return dir; };
    __host__ __device__ point3 at(double t) const;
    __host__ __device__ inter_ findIntersection(const Plane& plane) const;
    point3 orig;
    vec3 dir;
};
__host__ __device__ std::ostream& operator<<(std::ostream& out, const ray& r);

//bool Face_hit(const Plane& pl, const ray& r,const vector<int>& polygon, map<int, point3> vertices_coords);
__global__ void Face_hit(float* pl, const ray& r, int* polygon_langht, int* polygon, float* vertices_coords, int Face_number, int* start_index_per_face);


#endif