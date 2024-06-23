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
#include "Material.cuh"
#include "SaveAsBMP.cuh"
#include<cmath>


int main() {

    ///////////////////////////////////////////////
    //
    //wczytanie obiektów .obj
    //
    //////////////////////////////////////////////////


    //int vert_num1 = 3;
    //int face_num1 = 1;
    //int normal_num1 = 1;

    //Pointer_storage Kostka = GetDataFromObj(vert_num1, face_num1, normal_num1, "../../../helpers/untitled.obj");

    //int vert_num2 = 3;
    //int face_num2 = 1;
    //int normal_num2 = 1;

    //Pointer_storage Trojkatna_Kostka = GetDataFromObj(vert_num2, face_num2, normal_num2, "../../../helpers/Trojkatny_szescian.obj");
    int object_couter = 0;
    int vert_num3 = 3;
    int face_num3 = 1;
    int normal_num3 = 1;

    //Pointer_storage pociety_walec = GetDataFromObj(vert_num3, face_num3, normal_num3, "../../../helpers/walec_ale_kanciastyXD.obj");
    Pointer_storage pociety_walec = GetDataFromObj(vert_num3, face_num3, normal_num3, object_couter, "../../../helpers/scena_jeszce_raz.obj");

    //cout << endl << endl << "Kostka" << endl << endl;

    //cout << endl << endl << "Kostka" << endl << endl;
    //cout << endl << endl << "Kostka" << endl << endl;






    //Print_Import_data(Kostka, vert_num1, normal_num1, face_num1);

    //cout << endl << endl << "Sześcian pokrojony" << endl << endl;
    //Print_Import_data(Trojkatna_Kostka, vert_num2, normal_num2, face_num2);

    //cout << endl << endl << "zlosliwy przyklad kanciastego walca" << endl << endl;
    //Print_Import_data(pociety_walec, vert_num3, normal_num3, face_num3);

    //cout << endl << endl << "Testowanie promieni";



    ///////////////////////////////////////////////
    //
    //wybór obiektu
    //
    //////////////////////////////////////////////////
    int Vert_NUM = vert_num3;
    int Face_NUM = face_num3;
    int Normal_NUM = normal_num3;


    cout << object_couter << endl;

    Pointer_storage liczony_objekt = pociety_walec;
    ///////////////////////////////////////////////
    //
    //przepisanie wskaźników
    //
    //////////////////////////////////////////////////
    float** Planes = new float* [Face_NUM];
    Planes[0] = new float[Face_NUM * 4];

    float** Verticies = new float* [Vert_NUM];
    Verticies[0] = new float[Vert_NUM * 3];

    float** Normals = new float* [Normal_NUM];
    Normals[0] = new float[Normal_NUM * 3];

    float* Distances = new float[WIDTH * HEIGHT * Face_NUM];
    float* Colors = new float[WIDTH * HEIGHT * 3];
    float* shadows = new float[WIDTH * HEIGHT * 3];

    int* number_of_vertices_in_one_face = liczony_objekt.Face_size;
    int* normal_index_to_face = liczony_objekt.Face_to_Normal;
    int* Object_to_Face = liczony_objekt.Object_to_Face;
    int* Object_to_Vertex = liczony_objekt.Object_to_Vertex;
    int* start_face_at_index = new int[Face_NUM];
    float* Intersections = new float[WIDTH * HEIGHT * 3];
    start_face_at_index[0] = 0;
    int Length_to_Allocate_Faces = 0;
    for (int i = 0; i < Face_NUM; i++) {
        Length_to_Allocate_Faces += number_of_vertices_in_one_face[i];
    }

    int** Faces = new int* [Face_NUM];
    Faces[0] = new int[Length_to_Allocate_Faces];

    int current_index = 0;
    for (int i = 1; i < Face_NUM + 1; i++) {
        Faces[i] = Faces[0] + current_index + number_of_vertices_in_one_face[i - 1];
        start_face_at_index[i - 1] = current_index;
        current_index += number_of_vertices_in_one_face[i - 1];
        Planes[i] = Planes[0] + i * 4;


        //cout << i-1 << "  :  " << start_face_at_index[i-1] << endl;
    }

    for (int i = 1; i < Vert_NUM; i++) {
        Verticies[i] = Verticies[0] + i * 3;
    }
    for (int i = 1; i < Normal_NUM; i++) {
        Normals[i] = Normals[0] + i * 3;
    }
    for (int i = 0; i < Face_NUM; i++) {
        for (int j = 0; j < 4; j++)
        {
            Planes[i][j] = liczony_objekt.Planes[i][j];
        }
    }
    for (int i = 0; i < Vert_NUM; i++) {
        for (int j = 0; j < 3; j++)
            Verticies[i][j] = liczony_objekt.Vertices[i][j];
    }
    for (int i = 0; i < Normal_NUM; i++) {
        for (int j = 0; j < 3; j++)
            Normals[i][j] = liczony_objekt.Normals[i][j];
    }
    for (int i = 0; i < Face_NUM; i++) {
        for (int j = 0; j < liczony_objekt.Face_size[i]; j++) {
            Faces[i][j] = liczony_objekt.Faces[i][j];
        }
    }

    ///////////////////////////////////////////////
    //
    //utworzenie materiałów
    //
    //////////////////////////////////////////////////
    Material* Materials = new Material[4];
    Materials[0] = Material(0.8f, 0.1f, 0.1f,   // diffuse (red)
        0.9f, 0.6f, 0.6f,   // specular (light red)
        0.3f, 0.1f, 0.1f,   // ambient (dark red)
        1.0f,               // alpha
        4.0f,  	        // shininess  
        0.5f);             // reflectivity

    Materials[1] = Material(0.1f, 0.8f, 0.1f,   // diffuse (green)
        0.6f, 0.9f, 0.6f,   // specular (light green)
        0.1f, 0.3f, 0.1f,   // ambient (dark green)
        1.0f,               // alpha
        4.0f,  	        // shininess  
        0.5f);             // reflectivity

    Materials[2] = Material(0.1f, 0.1f, 0.8f,   // diffuse (blue)
        0.9f, 0.9f, 0.9f,   // specular (light blue)
        0.1f, 0.1f, 0.3f,   // ambient (dark blue)
        1.0f,               // alpha
        4.0f,  	        // shininess  
        0.5f);             // reflectivity

    Materials[3] = Material(0.8f, 0.8f, 0.1f,   // diffuse (yellow)
        0.9f, 0.9f, 0.6f,   // specular (light yellow)
        0.3f, 0.3f, 0.1f,   // ambient (dark yellow)
        1.0f,               // alpha
        4.0f,  	        // shininess  
        0.5f);             // reflectivity

    ///////////////////////////////////////////////
    //
    //Alokacja i przekazanie danych do karty
    //
    //
    //
    //
    //////////////////////////////////////////////////

    int* d_Faces;
    int* d_number_of_vertices_in_one_face;
    int* d_normal_index_to_face;
    int* d_start_face_at_index;
    int* d_Object_to_Vertex;
    int* d_Object_to_Face;
    float* d_distances;
    float* d_Vertices;
    float* d_Normals;
    float* d_Planes;
    float* d_closest_interesections;
    float* d_shadows;
    Material* d_Materials;

    cudaMalloc(&d_Faces, Length_to_Allocate_Faces * sizeof(int));
    cudaMalloc(&d_Planes, 4 * Face_NUM * sizeof(float));
    cudaMalloc(&d_Normals, 3 * Normal_NUM * sizeof(float));
    cudaMalloc(&d_Vertices, 3 * Vert_NUM * sizeof(float));
    cudaMalloc(&d_number_of_vertices_in_one_face, Face_NUM * sizeof(int));
    cudaMalloc(&d_normal_index_to_face, Face_NUM * sizeof(int));
    cudaMalloc(&d_start_face_at_index, Face_NUM * sizeof(int));
    cudaMalloc(&d_distances, WIDTH * HEIGHT * Face_NUM * sizeof(float));
    cudaMalloc(&d_closest_interesections, WIDTH * HEIGHT * 3 * sizeof(float));
    cudaMalloc(&d_Object_to_Vertex, Vert_NUM * sizeof(int));
    cudaMalloc(&d_Object_to_Face, Face_NUM * sizeof(int));
    cudaMalloc(&d_shadows, WIDTH* HEIGHT * 3 * sizeof(float));


    cudaMemcpy(d_Faces, Faces[0], Length_to_Allocate_Faces * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Planes, Planes[0], 4 * Face_NUM * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Normals, Normals[0], 3 * Face_NUM * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Vertices, Verticies[0], 3 * Vert_NUM * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_number_of_vertices_in_one_face, number_of_vertices_in_one_face, Face_NUM * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_normal_index_to_face, normal_index_to_face, Face_NUM * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_start_face_at_index, start_face_at_index, Face_NUM * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Object_to_Vertex, Object_to_Vertex, Vert_NUM * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Object_to_Face, Object_to_Face, Face_NUM * sizeof(int), cudaMemcpyHostToDevice);

    ///////////////////////////////////////////////
    //
    //Tranformacja obiektu rotacja/skala/przesunięcie
    //
    //////////////////////////////////////////////////

    int threadsPerBlock = 256;
    int blocksPerGrid = (Vert_NUM + threadsPerBlock - 1) / threadsPerBlock;
    float TranslateX = 0.0f;
    float TranslateY = 0.0f;
    float TranslateZ = 0.0f;
    float rotateX = 0.0f;
    float rotateY = 0.0f;
    float rotateZ = 0.0f;
    float scaleX = 1.0f;
    float scaleY = -1.0f;
    float scaleZ = -1.0f;

    int index = 0;
    Transform << <blocksPerGrid, threadsPerBlock >> > (d_Vertices, Vert_NUM, d_Object_to_Vertex, index, TranslateX, TranslateY, TranslateZ, rotateX, rotateY, rotateZ, scaleX, scaleY, scaleZ);
    cudaDeviceSynchronize();
    index = 1;
    TranslateZ = 5.0f;
    TranslateY = 0.3f;


    Transform << <blocksPerGrid, threadsPerBlock >> > (d_Vertices, Vert_NUM, d_Object_to_Vertex, index, TranslateX, TranslateY, TranslateZ, rotateX, rotateY, rotateZ, scaleX, scaleY, scaleZ);
    cudaDeviceSynchronize();
    index = 2;
    TranslateZ = 0.0f;

    Transform <<< blocksPerGrid, threadsPerBlock >> > (d_Vertices, Vert_NUM, d_Object_to_Vertex, index, TranslateX, TranslateY, TranslateZ, rotateX, rotateY, rotateZ, scaleX, scaleY, scaleZ);
    cudaDeviceSynchronize();
    index = 3;
    TranslateY = -1.0f;

    Transform << <blocksPerGrid, threadsPerBlock >> > (d_Vertices, Vert_NUM, d_Object_to_Vertex, index, TranslateX, TranslateY, TranslateZ, rotateX, rotateY, rotateZ, scaleX, scaleY, scaleZ);
    cudaDeviceSynchronize();

    threadsPerBlock = 256;
    blocksPerGrid = (Face_NUM + threadsPerBlock - 1) / threadsPerBlock;
    Update_normals_and_Planes << <blocksPerGrid, threadsPerBlock >> > (d_Vertices, d_Faces, d_Object_to_Face, index, d_Normals, d_Planes, d_number_of_vertices_in_one_face, d_normal_index_to_face, d_start_face_at_index, Face_NUM, Normal_NUM);
    cudaDeviceSynchronize();

    ///////////////////////////////////////////////
    //
    //  Wyznaczanie Promieni, uderzenia i dystanse
    //
    //////////////////////////////////////////////////
    int reflecions = 0;
    double focal_length = 10;
    point3 h_camera_center(120.0, -60.0, -50.0);
    point3  h_camera_focal(60.0 / 2, -25.0 / 2, -25.0 / 2);
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



    dim3 blockDim(8, 8, 8);
    dim3 gridDim((WIDTH + blockDim.x - 1) / blockDim.x,
        (HEIGHT + blockDim.y - 1) / blockDim.y,
        (Face_NUM + blockDim.z - 1) / blockDim.z);

    ///////////////////////////////////////////////
    //
    // wyznaczenie promieni
    //
    ///////////////////////////////////////////////
    cudaMemcpy(d_ray, h_ray[0], WIDTH * HEIGHT * sizeof(ray), cudaMemcpyHostToDevice);


    Generate_rays <<< gridDim, blockDim >>> (d_ray, focal_length, d_camera_center, d_camera_focal);
    cudaDeviceSynchronize();

    ///////////////////////////////////////////////
    //
    // wyznaczenie wektora dystansów
    //
    ///////////////////////////////////////////////

    Generate_distances << <gridDim, blockDim >> > (d_ray, d_camera_center, d_closest_interesections, d_normal_index_to_face, d_number_of_vertices_in_one_face,
        d_Faces, d_Vertices, d_Normals, d_Planes, d_start_face_at_index, Face_NUM, Vert_NUM, Normal_NUM, d_distances, reflecions);
    cudaDeviceSynchronize();
    reflecions += 1;




    ///////////////////////////////////////////////
    //
    // dobór najbliższej ściany z wektora dystansów
    //
    ////////////////////////////////////////////////


    int* d_close_indexes;
    float* d_colors;

    cudaMalloc(&d_Materials, object_couter * sizeof(Material));
    cudaMalloc(&d_close_indexes, WIDTH * HEIGHT * sizeof(int));
    cudaMalloc(&d_colors, WIDTH * HEIGHT * 3 * sizeof(float));


    cudaMemcpy(d_Materials, Materials, object_couter * sizeof(Material), cudaMemcpyHostToDevice);

    Choose_closest <<< gridDim, blockDim >>> (d_distances, Face_NUM, d_colors, d_Planes, d_Materials, d_ray, d_Object_to_Face, d_close_indexes, d_closest_interesections);
    cudaDeviceSynchronize();

    cudaMemcpy(Colors, d_colors, WIDTH * HEIGHT * 3 * sizeof(float), cudaMemcpyDeviceToHost);

    saveAsBMP(Colors, WIDTH, HEIGHT, "normals_image.bmp");
    ///////////////////////////////////////////////
    //
    // Dodanie cienia
    //
    ////////////////////////////////////////////////


    Add_shadows <<< gridDim, blockDim >> > (d_closest_interesections, d_shadows, d_normal_index_to_face, d_number_of_vertices_in_one_face, d_Faces, d_Vertices, d_Normals, d_Planes, d_start_face_at_index, Face_NUM, Vert_NUM, Normal_NUM);
    cudaDeviceSynchronize();

    cudaMemcpy(shadows, d_shadows, WIDTH * HEIGHT * 3 * sizeof(float), cudaMemcpyDeviceToHost);


    saveAsBMP(shadows, WIDTH, HEIGHT, "shadowes_image.bmp");

    for (int i = 0; i < WIDTH * HEIGHT * 3; i++)
    {
        if (shadows[i] >= 1)
        {
            Colors[i] -= 0.1f;
            if (Colors[i] < 0) { Colors[i] = 0.0f; }
        }

    }
    saveAsBMP(Colors, WIDTH, HEIGHT, "shadowed_scene_image.bmp");


    ///////////////////////////////////////////////
    //
    // druga iteracja promieni
    //
    ////////////////////////////////////////////////
    Update_rays << < gridDim, blockDim >>> (d_ray, d_closest_interesections, d_close_indexes, d_normal_index_to_face, d_Normals);
    cudaDeviceSynchronize();


    Generate_distances << < gridDim, blockDim >> > (d_ray, d_camera_center, d_closest_interesections, d_normal_index_to_face, d_number_of_vertices_in_one_face,
        d_Faces, d_Vertices, d_Normals, d_Planes, d_start_face_at_index, Face_NUM, Vert_NUM, Normal_NUM, d_distances, reflecions);
    cudaDeviceSynchronize();
    reflecions += 1;


    cudaMemcpy(d_Materials, Materials, object_couter * sizeof(Material), cudaMemcpyHostToDevice);

    Choose_closest << <gridDim, blockDim >> > (d_distances, Face_NUM, d_colors, d_Planes, d_Materials, d_ray, d_Object_to_Face, d_close_indexes, d_closest_interesections);
    cudaDeviceSynchronize();


    cudaMemcpy(Colors, d_colors, WIDTH * HEIGHT * 3 * sizeof(float), cudaMemcpyDeviceToHost);
    saveAsBMP(Colors, WIDTH, HEIGHT, "reflected_image.bmp");


    cudaFree(d_ray);
    cudaFree(d_distances);
    cudaFree(d_camera_center);
    cudaFree(d_camera_focal);
    cudaFree(d_closest_interesections);

    free(h_ray[0]);
    free(h_ray);
    delete[] Colors;

    return 0;
}