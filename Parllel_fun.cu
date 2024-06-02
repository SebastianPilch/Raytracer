
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include "vec3.cuh"
#include "Parllel_fun.cuh"
#include "Ray.cuh"

using namespace std;
#ifndef __CUDACC__ 
#define __CUDACC__
#endif




//std::vector<std::vector<float>> img(img_height, std::vector<float>(img_width));
//
//for (int j = 0; j < img_height; j++) {
//    for (int i = 0; i < img_width; i++) {
//        auto pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
//        auto ray_direction = pixel_center - camera_center;
//        ray newRay = ray(camera_center, ray_direction);
//
//        float hit_anything = 0.0;
//        for (const auto& pair : Planes_to_faces) {
//            bool hit = Face_hit(pair.second, newRay, vertices_to_faces[pair.first], vertices_coors);
//            if (hit) {
//                hit_anything += 0.2;
//            }
//
//        }
//
//        img[j][i] = hit_anything;
//    }
//}

//void saveAsBMP(const std::vector<std::vector<float>>& img, int width, int height, const std::string& filename) {
//    std::ofstream file(filename, std::ios::out | std::ios::binary);
//
//    if (!file) {
//        std::cerr << "Cannot open file: " << filename << std::endl;
//        return;
//    }
//
//    int paddingSize = (4 - (width * 3) % 4) % 4; // Padding wymagany przez format BMP
//
//    // Nag³ówek BMP
//    int filesize = 54 + (3 * width + paddingSize) * height;
//    char fileHeader[54] = { 'B', 'M', 0,0,0,0, 0,0, 0,0, 54,0,0,0, 40,0,0,0, static_cast<char>(width), static_cast<char>(width >> 8), static_cast<char>(width >> 16), static_cast<char>(width >> 24), static_cast<char>(height), static_cast<char>(height >> 8), static_cast<char>(height >> 16), static_cast<char>(height >> 24), 1,0, 24,0, 0,0,0,0, static_cast<char>(filesize), static_cast<char>(filesize >> 8), static_cast<char>(filesize >> 16), static_cast<char>(filesize >> 24), 0,0,0,0, 0,0,0,0 };
//
//    // Zapisanie nag³ówka
//    file.write(fileHeader, 54);
//
//    // Zapisanie danych pikseli
//    for (int i = height - 1; i >= 0; i--) {
//        for (int j = 0; j < width; j++) {
//            unsigned char color = static_cast<unsigned char>(img[i][j] * 255); // Skalowanie wartoœci z [0, 1] do [0, 255]
//            file.put(color);
//            file.put(color);
//            file.put(color);
//        }
//        // Dodanie paddingu
//        for (int k = 0; k < paddingSize; k++) {
//            file.put(0);
//        }
//    }
//
//    file.close();
//}


__global__ void Generate_rays(ray* viewport_rays,double focal_length, point3* camera_center,
    point3* camera_focal, int* d_normal_index_to_face, int* d_number_of_vertices_in_one_face,
    int* d_Faces, float* d_Vertices, float* d_Normals, float* d_Planes,int* start_face_at_index 
    ,int Face_NUM, int Vertex_NUM, int Normal_NUM, float* d_distances)
{
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

    
    Plane plane = Plane((double)d_Planes[f*4], (double)d_Planes[f * 4+1], (double)d_Planes[f * 4+2], (double)d_Planes[f * 4+3]);
    point3 intersection = UV_ray.findIntersection(plane);
    bool Is_Hitten_correct = true;
 
    if (intersection[0] == INFINITY || intersection[0] == -INFINITY || intersection[1] == INFINITY ||
        intersection[1] == -INFINITY == -INFINITY || intersection[2] == INFINITY || intersection[2] == -INFINITY) 
    {
        Is_Hitten_correct = false;
    }
    
    vec3 edge;
    size_t next_index;
    for (size_t i = start_face_at_index[f]; i < start_face_at_index[f] + d_number_of_vertices_in_one_face[f]; ++i) 
    {
        size_t next_index;
        if (i + 1 > start_face_at_index[f] + d_number_of_vertices_in_one_face[f])
        {
            next_index = i + 1;
        }
        else
        {
            next_index = start_face_at_index[f];
        }
        int vertex_index = d_Faces[i] - 1;
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
    }
    else
    {
        d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] =(float)plane.A;//sqrt(dis[0] * dis[0] + dis[1] * dis[1] + dis[2] * dis[2]);
    }
    __syncthreads();
//


        //Face_hit << <Cuda_Blocks, Threads_in_one_block >> > (d_Planes, UV_ray, d_number_of_vertices_in_one_face, d_Faces, d_Vertices, Face_NUM, start_face_at_index);




}


