#include <map>
#include <vector>
#include "vec3.h"
#include "Ray.h"
using namespace std;
double positiveInfinity = INFINITY;
double negativeInfinity = -INFINITY;
ray::ray() {};
ray::ray(const point3& origin, const vec3& direction) : orig(origin), dir(direction) {};
point3 ray::at(double t) const { return this->dir*t + this->orig; }
point3 ray::findIntersection(const Plane& plane) const
{
    cout << this->dir << endl;
    double t = -(plane.A * this->orig.x() + plane.B * this->orig.y() + plane.C * this->orig.z() + plane.D) /
        (plane.A * this->dir.x() + plane.B * this->dir.y() + plane.C * this->dir.z());

    point3 intersection(at(t));
    cout << intersection << endl;

    return intersection;
}

bool Face_hit(const Plane& pl, const ray& r, const vector<int>& polygon, map<int, point3> vertices_coords)
{
    point3 intersection = r.findIntersection(pl);
    if (intersection[0] == INFINITY || intersection[0] == -INFINITY || intersection[1] == INFINITY ||
        intersection[1] == -INFINITY == -INFINITY || intersection[2] == INFINITY || intersection[2] == -INFINITY) {
        return false;
    }
    vec3 edge;
    for (size_t i = 0; i < polygon.size(); ++i) {

        edge = vertices_coords[polygon[i + 1]] - vertices_coords[polygon[i]];
        vec3 vp = intersection - vertices_coords[polygon[i]];
        vec3 n = crossProduct_(edge, vp);
        vec3 normal = vec3(pl.A, pl.B, pl.C);
        if (dotProduct_(n, normal) < 0) {
            return false;
        }
    }
    return true;

};