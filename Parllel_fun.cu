#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include "vec3.cuh"
#include "Parllel_fun.cuh"
#include "Ray.cuh"

using namespace std;
#ifndef __CUDACC__ 
#define __CUDACC__
#endif

__global__ void Generate_rays(ray* viewport_rays, double focal_length, point3* camera_center,
    point3* camera_focal, int* d_normal_index_to_face, int* d_number_of_vertices_in_one_face,
    int* d_Faces, float* d_Vertices, float* d_Normals, float* d_Planes, int* start_face_at_index,
    int Face_NUM, int Vertex_NUM, int Normal_NUM, float* d_distances, float* d_closest_angles) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int f = blockIdx.z * blockDim.z + threadIdx.z;
    if (f >= Face_NUM || i >= WIDTH || j >= HEIGHT) return;

    vec3 viewport_upper_left = *camera_center - vec3(0, 0, focal_length) - VIEWPORT_U / 2 - VIEWPORT_V / 2;
    vec3 pixel00_loc = viewport_upper_left + 0.5 * (DELTA_U + DELTA_V);

    vec3 pixel_center = pixel00_loc + (i * DELTA_U) + (j * DELTA_V);
    vec3 ray_direction = pixel_center - *camera_center;
    ray UV_ray = ray(*camera_center, ray_direction);
    viewport_rays[j * HEIGHT + i] = ray(*camera_center, ray_direction);

    __syncthreads();

    Plane plane = Plane((double)d_Planes[f * 4], (double)d_Planes[f * 4 + 1], (double)d_Planes[f * 4 + 2], (double)d_Planes[f * 4 + 3]);
    point3 intersection = UV_ray.findIntersection(plane);
    bool Is_Hitten_correct = true;
    if (intersection[0] == INFINITY || intersection[0] == -INFINITY || intersection[1] == INFINITY ||
        intersection[1] == -INFINITY || intersection[2] == INFINITY || intersection[2] == -INFINITY)
    {
        Is_Hitten_correct = false;
    }

    vec3 edge;
    size_t next_index;

    for (size_t idx = start_face_at_index[f]; idx < start_face_at_index[f] + d_number_of_vertices_in_one_face[f]; ++idx)
    {
        size_t next_index;
        if (idx + 1 > start_face_at_index[f] + d_number_of_vertices_in_one_face[f])
        {
            next_index = idx + 1;
        }
        else
        {
            next_index = start_face_at_index[f];
        }
        int vertex_index = d_Faces[idx] - 1;
        int vertex_next = d_Faces[next_index] - 1;

        point3 next_vertex = point3(
            (double)d_Vertices[3 * vertex_index],
            (double)d_Vertices[3 * vertex_index + 1],
            (double)d_Vertices[3 * vertex_index + 2]
        );

        point3 current_vertex = point3(
            (double)d_Vertices[3 * vertex_next],
            (double)d_Vertices[3 * vertex_next + 1],
            (double)d_Vertices[3 * vertex_next + 2]
        );

        edge = next_vertex - current_vertex;
        vec3 vp = intersection - current_vertex;
        vec3 n = crossProduct_(edge, vp);
        vec3 normal = vec3(plane.A, plane.B, plane.C);
        if (dotProduct_(n, normal) < 0) {
            Is_Hitten_correct = false;
        }
    }
    vec3 dis = *camera_center - intersection;
    if (!Is_Hitten_correct)
    {
        d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] = INFINITY;
        d_closest_angles[(j * WIDTH + i) * 3 + 0] = 0.0f;
        d_closest_angles[(j * WIDTH + i) * 3 + 1] = 0.0f;
        d_closest_angles[(j * WIDTH + i) * 3 + 2] = 0.0f;
    }
    else
    {
        float distance = sqrt(dis[0] * dis[0] + dis[1] * dis[1] + dis[2] * dis[2]);
        d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] = distance;

        vec3 normal = vec3((float)plane.A, (float)plane.B, (float)plane.C);
        float norm_length = sqrt(normal.x() * normal.x() + normal.y() * normal.y() + normal.z() * normal.z());
        normal = normal / norm_length;

        d_closest_angles[(j * WIDTH + i) * 3 + 0] = (normal.x() + 1.0f) * 0.5f;
        d_closest_angles[(j * WIDTH + i) * 3 + 1] = (normal.y() + 1.0f) * 0.5f;
        d_closest_angles[(j * WIDTH + i) * 3 + 2] = (normal.z() + 1.0f) * 0.5f;
    }
    __syncthreads();
}
