
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
    ,int Face_NUM, int Vertex_NUM, int Normal_NUM)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    vec3 viewport_upper_left = *camera_center - vec3(0, 0, focal_length) - VIEWPORT_U / 2 - VIEWPORT_V / 2;
    vec3 pixel00_loc = viewport_upper_left + 0.5 * (DELTA_U + DELTA_V);

    vec3 pixel_center = pixel00_loc + (i * DELTA_U) + (j * DELTA_V);
    vec3 ray_direction = pixel_center - *camera_center;
    viewport_rays[j*HEIGHT + i] = ray(*camera_center , ray_direction);



}




__global__ void Face_hit()
{






};