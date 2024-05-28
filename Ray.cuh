#ifndef RAY_H
#define RAY_H
#include <map>
#include <vector>
#include "vec3.cuh"
using namespace std;

class ray {
public:
     ray();
     ray(const point3& origin, const vec3& direction);
     const point3& origin() { return orig; };
     const vec3& direction() { return dir; };
     point3 at(double t) const;
     point3 findIntersection(const Plane& plane) const;
private:
    point3 orig;
    vec3 dir;
};

//bool Face_hit(const Plane& pl, const ray& r,const vector<int>& polygon, map<int, point3> vertices_coords);



#endif