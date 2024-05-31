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

__host__ __device__ std::ostream& operator<<(std::ostream& out, const ray& r) {
    return out <<'[' << r.dir <<"]*t + [" << r.orig <<"]";
}


//bool Face_hit(float* pl, const ray& r, int polygon_langht ,int* polygon, float** vertices_coords)
//{
//    Plane plane = Plane((double)pl[0], (double)pl[1], (double)pl[2], (double)pl[3]);
//    point3 intersection = r.findIntersection(plane);
//    if (intersection[0] == INFINITY || intersection[0] == -INFINITY || intersection[1] == INFINITY ||
//        intersection[1] == -INFINITY == -INFINITY || intersection[2] == INFINITY || intersection[2] == -INFINITY) {
//        return false;
//    }
//    vec3 edge;
//    for (size_t i = start_index_per_face; i < polygon_langht; ++i) {
//        size_t next_index = (i + 1) % polygon_langht;
//        point3 current_vertex = point3((double)vertices_coords[polygon[i] - 1][0], (double)vertices_coords[polygon[i] - 1][1], (double)vertices_coords[polygon[i] - 1][2]);
//        point3 next_vertex = point3((double)vertices_coords[polygon[next_index] - 1][0], (double)vertices_coords[polygon[next_index] - 1][1], (double)vertices_coords[polygon[next_index] - 1][2]);
//        edge = next_vertex - current_vertex;
//        vec3 vp = intersection - current_vertex;
//        vec3 n = crossProduct_(edge, vp);
//        vec3 normal = vec3(plane.A, plane.B, plane.C);
//        if (dotProduct_(n, normal) < 0) {
//            return false;
//        }
//    }
//
//
//};

__global__ void Face_hit(float* pl, const ray& r, int* polygon_langht, int* polygon, float* vertices_coords, int Face_number, int* start_index_per_face) {

    extern __shared__ float distances[];

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int tid = threadIdx.x;
    for (int i = tid; i < Face_number; i += blockDim.x) {
        distances[i] = 0.0f;
    }
    __syncthreads();


    if (idx >= Face_number) return;

    Plane plane = Plane((double)pl[0], (double)pl[1], (double)pl[2], (double)pl[3]);
    point3 intersection = r.findIntersection(plane);
    if (intersection[0] == INFINITY || intersection[0] == -INFINITY || intersection[1] == INFINITY ||
        intersection[1] == -INFINITY == -INFINITY || intersection[2] == INFINITY || intersection[2] == -INFINITY) {
    }

    vec3 edge;
    for (size_t i = start_index_per_face[idx]; i < start_index_per_face[idx] + polygon_langht[idx]; ++i) {


        size_t next_index;
        if (i + 1 > start_index_per_face[idx] + polygon_langht[idx])
        {
            next_index = i + 1;

        }
        else
        {
            next_index = start_index_per_face[idx];

        }

        int vertex_index = polygon[i] - 1;
        int vertex_next = polygon[next_index] - 1;

        point3 next_vertex = point3(
            (double)vertices_coords[3 * vertex_index],
            (double)vertices_coords[3 * vertex_index + 1],
            (double)vertices_coords[3 * vertex_index + 2]
        );

        point3 current_vertex = point3(
            (double)vertices_coords[3 * vertex_next],
            (double)vertices_coords[3 * vertex_next + 1],
            (double)vertices_coords[3 * vertex_next + 2]
        );


        edge = next_vertex - current_vertex;
        vec3 vp = intersection - current_vertex;
        vec3 n = crossProduct_(edge, vp);
        vec3 normal = vec3(plane.A, plane.B, plane.C);
        if (dotProduct_(n, normal) < 0) {
            return;
        }
    }

}