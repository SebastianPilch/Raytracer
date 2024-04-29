#include <map>
#include <vector>
#include "vec3.h"
#include "Ray.h"
using namespace std;

ray::ray() {};
ray::ray(const point3& origin, const vec3& direction) : orig(origin), dir(direction) {};
point3 ray::at(double t) const { return orig + t * dir; }
point3 ray::findIntersection(const Plane& plane) const
{
    double t = -(plane.A * orig.x() + plane.B * orig.y() + plane.C * orig.z() + plane.D) /
        (plane.A * dir.x() + plane.B * dir.y() + plane.C * dir.z());
    point3 intersection(at(t));
    return intersection;
}

bool Face_hit(const Plane& pl, const ray& r, const vector<int>& polygon, map<int, point3> vertices_coords)
{
    point3 intersection = r.findIntersection(pl);
    vec3 edge;
    for (size_t i = 0; i < polygon.size()-1; ++i) {

        edge = vertices_coords[polygon[i+1]] - vertices_coords[polygon[i]];
        vec3 vp = intersection - vertices_coords[polygon[i]];
        vec3 n = crossProduct_(edge, vp);
        vec3 normal = vec3(pl.A,pl.B,pl.C);
        if (dotProduct_(n, normal) < 0) {
            return false;
        }
    }
    return true;

};

