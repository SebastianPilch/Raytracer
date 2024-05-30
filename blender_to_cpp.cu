
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <filesystem>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "vec3.cuh"
#include "ImportObj.cuh"
#include "Ray.cuh"
#include "Parllel_fun.cuh"
//



int main() {


    int vert_num1 = 3;
    int face_num1 = 1;
    int normal_num1 = 1;

    Pointer_storage Kostka = GetDataFromObj(vert_num1,
        face_num1, normal_num1, "../../../helpers/untitled.obj");

    int vert_num2 = 3;
    int face_num2 = 1;
    int normal_num2 = 1;

    Pointer_storage Trojkatna_Kostka = GetDataFromObj(vert_num2,
        face_num2, normal_num2, "../../../helpers/Trojkatny_szescian.obj");



    int vert_num3 = 3;
    int face_num3 = 1;
    int normal_num3 = 1;

    Pointer_storage pociety_walec = GetDataFromObj(vert_num3,
        face_num3, normal_num3, "../../../helpers/walec_ale_kanciastyXD.obj");


    cout << endl << endl << "Kostka" << endl << endl;

    Print_Import_data(Kostka, vert_num1, normal_num1, face_num1);

    cout << endl << endl << "Sześcian pokrojony" << endl << endl;
    Print_Import_data(Trojkatna_Kostka, vert_num2, normal_num2, face_num2);

    cout << endl << endl << "zlosliwy przyklad kanciastego walca" << endl << endl;

    Print_Import_data(pociety_walec, vert_num3, normal_num3, face_num3);

    //saveAsBMP(img, width, height, "result_image.bmp");

    cout << "XDD" << endl;
    cout << endl << endl << " Linia przed cuda";


    const int size = 10;
    Vector* d_vectors;
    vec3* d_z = (vec3*)malloc(size * sizeof(vec3));
    vec3* z = (vec3*)malloc(size * sizeof(vec3));

    cudaMalloc(&d_vectors, size * sizeof(Vector));
    cudaMalloc(&d_z, size * sizeof(vec3));
    int threadsPerBlock = 512;
    int blocksPerGrid = (size + threadsPerBlock - 1) / threadsPerBlock;
    
    for (int i = 0; i < size; i++) {
        z[i] = vec3();
    }

    MyKernel <<<blocksPerGrid, threadsPerBlock >> > (d_vectors, size, d_z);
    Vector h_vectors[size];
    cudaMemcpy(h_vectors, d_vectors, size * sizeof(Vector), cudaMemcpyDeviceToHost);
    cudaMemcpy(z, d_z, size * sizeof(vec3), cudaMemcpyDeviceToHost);
    printVectors(h_vectors, size);
    for (int i = 0; i < size; i++) {
        cout << z[i] << endl;
    }
    cudaFree(d_vectors);
    cudaFree(d_z);



    cout << endl << endl << "Testowanie promieni" << WIDTH << HEIGHT;


    double focal_length = 10;
    point3 h_camera_center(5, 5, 5);
    point3 h_camera_focal(-5, -5, -5);
    point3* d_camera_center;
    point3* d_camera_focal;
    ray** h_ray;
    ray* d_ray;
    h_ray = (ray**)malloc(HEIGHT * sizeof(ray*));
    h_ray[0] = (ray*)malloc(HEIGHT * WIDTH * sizeof(ray));
    for (int i = 1; i < HEIGHT; i++) {
        h_ray[i] = h_ray[0] + i * WIDTH;
    }
    cudaMalloc(&d_ray, WIDTH * HEIGHT * sizeof(ray));
    cudaMalloc((void**)&d_camera_center, sizeof(point3));
    cudaMalloc((void**)&d_camera_focal, sizeof(point3));
    cudaMemcpy(d_camera_center, &h_camera_center, sizeof(point3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_camera_focal, &h_camera_focal, sizeof(point3), cudaMemcpyHostToDevice);
    dim3 threadsPerBlock2(16, 16);
    dim3 numBlocks((WIDTH + threadsPerBlock2.x - 1) / threadsPerBlock2.x, (HEIGHT + threadsPerBlock2.y - 1) / threadsPerBlock2.y);
    Generate_rays<<<numBlocks, threadsPerBlock2>>> (d_ray, focal_length, d_camera_center, d_camera_focal);
    cudaMemcpy(h_ray[0], d_ray, WIDTH * HEIGHT * sizeof(ray), cudaMemcpyDeviceToHost);
    cudaFree(d_ray);
    cudaFree(d_camera_center);
    cudaFree(d_camera_focal);
    free(h_ray[0]);
    free(h_ray);

    cout << endl << endl << "Bicie sciany" ;


    return 0;

}


