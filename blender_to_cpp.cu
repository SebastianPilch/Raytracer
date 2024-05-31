
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


    cout << endl << endl << "Testowanie promieni" ;












    //wybór obiektu
    int Vert_NUM = vert_num1;
    int Face_NUM = face_num1;
    int Normal_NUM = normal_num1;
    Pointer_storage liczony_objekt = Kostka;
    // // // //


    int** Faces = liczony_objekt.Faces;
    float** Verticies = liczony_objekt.Vertices;
    float** Normals = liczony_objekt.Normals;
    float** Planes = liczony_objekt.Planes;
    float* Distances = new float[WIDTH*HEIGHT*Face_NUM];
    // (face_idx * WIDTH * HEIGHT) + (h * WIDTH) + w;


    int* number_of_vertices_in_one_face = liczony_objekt.Face_size;
    int* normal_index_to_face = liczony_objekt.Face_to_Normal;

    int* start_face_at_index = new int[Face_NUM];

    start_face_at_index[0] = 0;
    int Length_to_Allocate_Faces = 0;
    for (int i = 0; i < Face_NUM; i++) { Length_to_Allocate_Faces += number_of_vertices_in_one_face[i]; }

    int* Faces_d = new int[Length_to_Allocate_Faces];
    float* Vertices_d, Normals_d, Planes_d;
    int current_index = 0;
    for (int i = 1; i < Face_NUM; i++)
    {
        Faces[i] = Faces[current_index] + number_of_vertices_in_one_face[i - 1];
        current_index += number_of_vertices_in_one_face[i - 1];
        start_face_at_index[i] = current_index;
        Planes[i] = Planes[0] + i * 4;

    }
    for (int i = 1; i < Vert_NUM; i++)
    {
        Verticies[i] = Verticies[0] + i * 3;
    }
    for (int i = 1; i < Normal_NUM; i++)
    {
        Normals[i] = Normals[0] + i * 3;
    }
    int* d_Faces;
    int* d_number_of_vertices_in_one_face;
    int* d_normal_index_to_face;
    int* d_start_face_at_index;
    float* d_distances;
    float* d_Vertices;
    float* d_Normals;
    float* d_Planes;


    cudaMalloc(&d_Faces, Length_to_Allocate_Faces * sizeof(int));
    cudaMalloc(&d_Planes, 4 * Face_NUM * sizeof(float));
    cudaMalloc(&d_Normals, 3 * Normal_NUM * sizeof(float));
    cudaMalloc(&d_Vertices, 3 * Vert_NUM * sizeof(float));
    cudaMalloc(&d_number_of_vertices_in_one_face, Face_NUM * sizeof(int));
    cudaMalloc(&d_normal_index_to_face, Face_NUM * sizeof(int));
    cudaMalloc(&d_start_face_at_index, Face_NUM * sizeof(int));
    cudaMalloc(&d_distances, WIDTH*HEIGHT*Face_NUM * sizeof(float));



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
    

    dim3 blockDim(512, 512, 512); 
    dim3 gridDim((WIDTH + blockDim.x - 1) / blockDim.x,
        (HEIGHT + blockDim.y - 1) / blockDim.y,
        (Face_NUM + blockDim.z - 1) / blockDim.z);

    Generate_rays<<< gridDim, blockDim >>> (d_ray, focal_length, d_camera_center, d_camera_focal, normal_index_to_face,number_of_vertices_in_one_face,
    Faces[0], Verticies[0], Normals[0], Planes[0], d_start_face_at_index, Face_NUM, Vert_NUM, Normal_NUM, d_distances);


    cudaMemcpy(h_ray[0], d_ray, WIDTH * HEIGHT * sizeof(ray), cudaMemcpyDeviceToHost);
    cudaMemcpy(Distances, d_distances, WIDTH * HEIGHT * sizeof(ray), cudaMemcpyDeviceToHost);

    cudaFree(d_ray);
    cudaFree(d_distances);
    cudaFree(d_camera_center);
    cudaFree(d_camera_focal);


    cout << endl << endl << "Bicie sciany" ;

    for (int i = 0; i < WIDTH; i++) 
    {
        for (int j = 0; j < HEIGHT; j++) 
        {
            for (int f = 0; f < Face_NUM; f++) 
            {
                if (Distances[(f * WIDTH * HEIGHT) + (j * WIDTH) + i] > 0)
                {
                    cout << "  " << Distances[(f * WIDTH * HEIGHT) + (j * WIDTH) + i] << "  ";
                }
            }
        
        }

    
    }




    free(h_ray[0]);
    free(h_ray);

    cout << endl<<Length_to_Allocate_Faces;























    return 0;

}


