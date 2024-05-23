#include "cuda_runtime.h"
#include "device_launch_parameters.h"
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

#define WIDTH 400
#define HEIGHT 400
#define VIEWPORT_HEIGHT  2.0
#define VIEWPORT_WIDTH VIEWPORT_HEIGHT * (double)(WIDTH / HEIGHT)
#define VIEWPORT_U vec3(VIEWPORT_WIDTH, 0, 0)
#define VIEWPORT_V vec3(0, -VIEWPORT_HEIGHT, 0)
#define DELTA_U VIEWPORT_U / WIDTH
#define DELTA_V VIEWPORT_V / HEIGHT

//__global__ void Matrixadd(int m, int n, float* d_x, float* d_y, float* d_z)
//{
//    int i = blockIdx.x * blockDim.x + threadIdx.x;
//    if (i < m * n)
//    {
//        d_z[i] = d_x[i] + d_y[i];
//    }
//
//}


//void MatrixAddition(int N, int M, float** x, float** y, float** z)
//{
//    x[0] = new float[M * N];
//    for (size_t i = 1; i < M; i++) x[i] = x[0] + i * N;
//    y[0] = new float[M * N];
//    for (size_t i = 1; i < M; i++) y[i] = y[0] + i * N;
//    z[0] = new float[M * N];
//    for (size_t i = 1; i < M; i++) z[i] = z[0] + i * N;
//    float* d_x, *d_y, *d_z;
//    cudaMalloc(&d_x, N * M * sizeof(float));
//    cudaMalloc(&d_y, N * M * sizeof(float));
//    cudaMalloc(&d_z, N * M * sizeof(float));
//    cudaMemcpy(d_x, x[0], N * M * sizeof(float), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_y, y[0], N * M * sizeof(float), cudaMemcpyHostToDevice);
//    Matrixadd << <M, N >> > (M, N, d_x, d_y, d_z);
//    cudaMemcpy(z[0], d_z, N * M * sizeof(float), cudaMemcpyDeviceToHost);
//    cudaFree(d_x);
//    cudaFree(d_y);
//    cudaFree(d_z);
//    delete[] x[0];
//    delete[] x;
//    delete[] y[0];
//    delete[] y;
//    delete[] z[0];
//    delete[] z;
//}

//int img_height = 400;
//int img_width = 400;
//
//auto focal_length = 1.0;
//auto viewport_height = 2.0;
//auto viewport_width = viewport_height * (double(img_width) / img_height);
//auto camera_center = point3(0, 0, 5);
//auto viewport_u = vec3(viewport_width, 0, 0);
//auto viewport_v = vec3(0, -viewport_height, 0);
//
//auto pixel_delta_u = viewport_u / img_width;
//auto pixel_delta_v = viewport_v / img_height;
//auto viewport_upper_left = camera_center - vec3(0, 0, focal_length) - viewport_u / 2 - viewport_v / 2;
//auto pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);
//
//
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



//__global__ void Image_matrix(int** img,double focal_length, point3 camera_center, point3 camera_focal)
//{
//    int i = threadIdx.x + blockIdx.x * blockDim.x;
//    int j = threadIdx.y + blockIdx.y * blockDim.y;

    //auto viewport_upper_left = camera_center - vec3(0, 0, focal_length) - VIEWPORT_U / 2 - VIEWPORT_V / 2;
    //auto pixel00_loc = viewport_upper_left + 0.5 * (DELTA_U + DELTA_V);

    //auto pixel_center = pixel00_loc + (i * DELTA_U) + (j * DELTA_V);
    //auto ray_direction = pixel_center - camera_center;
    //ray newRay = ray(camera_center, ray_direction);
    //float hit_anything = 0.0;

    //bool hit = Face_hit(pair.second, newRay, vertices_to_faces[pair.first], vertices_coors);
    //if (hit) 
    //{
    //   hit_anything += 0.2;
    //}
    //img[j][i] = hit_anything;

//}
//dim3 Image_size(width, height);
