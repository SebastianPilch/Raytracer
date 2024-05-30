#ifndef RAY_H
#define RAY_H
#include <map>
#include <vector>
#include "vec3.cuh"
using namespace std;

class ray {
public:
    __host__ __device__ ray();
    __host__ __device__ ray(const point3& origin, const vec3& direction);
    __host__ __device__ const point3& origin() { return orig; };
    __host__ __device__ const vec3& direction() { return dir; };
    __host__ __device__ point3 at(double t) const;
    __host__ __device__ point3 findIntersection(const Plane& plane) const;
private:
    point3 orig;
    vec3 dir;
};

//bool Face_hit(const Plane& pl, const ray& r,const vector<int>& polygon, map<int, point3> vertices_coords);



#endif