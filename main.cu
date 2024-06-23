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
#include "Parallel_functions.cuh"
#include "Material.cuh"
#include "SaveAsBMP.cuh"
#include<cmath>


int main() {

    ///////////////////////////////////////////////
    //
    //wczytanie obiektów .obj
    //
    //////////////////////////////////////////////////

    int object_couter = 0;
    int vert_num3 = 3;
    int face_num3 = 1;
    int normal_num3 = 1;

    Pointer_storage pociety_walec = GetDataFromObj(vert_num3, face_num3, normal_num3, object_couter, "../../../helpers/Studia_z_budownictwa_budujemy_mosty_2.obj");


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



    int* number_of_vertices_in_one_face = liczony_objekt.Face_size;
    int* normal_index_to_face = liczony_objekt.Face_to_Normal;
    int* Object_to_Face = liczony_objekt.Object_to_Face;
    int* Object_to_Vertex = liczony_objekt.Object_to_Vertex;

    float* Distances = new float[WIDTH * HEIGHT * Face_NUM];
    float* Colors = new float[WIDTH * HEIGHT * 3];
    float* shadows = new float[WIDTH * HEIGHT * 3];

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
    Material* Materials = new Material[object_couter];
    //trawa
    Materials[5] = Material(0.6f, 0.8f, 0.3f,   // diffuse (bright green)
        0.7f, 0.9f, 0.4f,   // specular (light green)
        0.3f, 0.4f, 0.2f,   // ambient (dark green)
        1.0f,               // alpha
        4.0f,              // shininess  
        0.0f);             // reflectivity

    //woda

    Materials[2] = Material(0.3f, 0.6f, 0.8f,   // diffuse (light blue)
        0.5f, 0.7f, 0.9f,   // specular (light blue)
        0.2f, 0.3f, 0.4f,   // ambient (dark blue)
        1.0f,               // alpha
        4.0f,              // shininess  
        0.3f);             // reflectivity

    //pnie

    Materials[4] = Material(0.5f, 0.3f, 0.1f,   // diffuse (brown)
        0.6f, 0.4f, 0.2f,   // specular (light brown)
        0.3f, 0.2f, 0.1f,   // ambient (dark brown)
        1.0f,               // alpha
        4.0f,              // shininess  
        0.0f);             // reflectivity

    //drzewa

    Materials[3] = Material(0.1f, 0.4f, 0.1f,   // diffuse (dark green)
        0.2f, 0.5f, 0.2f,   // specular (light green)
        0.05f, 0.2f, 0.05f, // ambient (very dark green)
        1.0f,               // alpha
        4.0f,              // shininess  
        0.0f);             // reflectivity

    //most

    Materials[1] = Material(0.6f, 0.6f, 0.6f,   // diffuse (metallic grey)
        0.8f, 0.8f, 0.8f,   // specular (bright grey)
        0.3f, 0.3f, 0.3f,   // ambient (dark grey)
        1.0f,               // alpha
        4.0f,              // shininess  
        0.15f);             // reflectivity
    //skała
    Materials[0] = Material(0.8f, 0.8f, 0.8f,   // diffuse (metallic grey)
        0.8f, 0.8f, 0.8f,   // specular (bright grey)
        0.3f, 0.3f, 0.3f,   // ambient (dark grey)
        1.0f,               // alpha
        1.0f,              // shininess  
        0.0f);             // reflectivity
    ///////////////////////////////////////////////
    //
    //Alokacja i przekazanie danych do karty
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
    float TranslateX = 13.0f;
    float TranslateY = 0.0f;
    float TranslateZ = -6.0f;
    float rotateX = 0.0f;
    float rotateY = 135.0f;
    float rotateZ = 0.0f;
    float scaleX = 1.5f;
    float scaleY = -1.5f;
    float scaleZ = -1.5f;

    int index = 0;
    Transform << <blocksPerGrid, threadsPerBlock >> > (d_Vertices, Vert_NUM, d_Object_to_Vertex, index, TranslateX, TranslateY, TranslateZ, rotateX, rotateY, rotateZ, scaleX, scaleY, scaleZ);
    cudaDeviceSynchronize();
    index = 1;
    Transform << <blocksPerGrid, threadsPerBlock >> > (d_Vertices, Vert_NUM, d_Object_to_Vertex, index, TranslateX, TranslateY, TranslateZ, rotateX, rotateY, rotateZ, scaleX, scaleY, scaleZ);
    cudaDeviceSynchronize();
    index = 2;
    Transform <<< blocksPerGrid, threadsPerBlock >> > (d_Vertices, Vert_NUM, d_Object_to_Vertex, index, TranslateX, TranslateY, TranslateZ, rotateX, rotateY, rotateZ, scaleX, scaleY, scaleZ);
    cudaDeviceSynchronize();
    index = 3;
    Transform << <blocksPerGrid, threadsPerBlock >> > (d_Vertices, Vert_NUM, d_Object_to_Vertex, index, TranslateX, TranslateY, TranslateZ, rotateX, rotateY, rotateZ, scaleX, scaleY, scaleZ);
    cudaDeviceSynchronize();
    index = 4;
    Transform << <blocksPerGrid, threadsPerBlock >> > (d_Vertices, Vert_NUM, d_Object_to_Vertex, index, TranslateX, TranslateY, TranslateZ, rotateX, rotateY, rotateZ, scaleX, scaleY, scaleZ);
    cudaDeviceSynchronize();
    index = 5;
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
    point3 h_camera_center(120.0, -70.0, -70.0);
    point3  h_camera_focal(120.0 / 2, -70.0 / 2, -70.0 / 2);
    vec3 world_light_dir = vec3(-0.7,1,-1);
    point3* d_camera_center;
    point3* d_camera_focal;
    vec3* d_world_light_dir;
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
    cudaMalloc((void**)&d_world_light_dir, sizeof(vec3));

    cudaMemcpy(d_camera_center, &h_camera_center, sizeof(point3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_camera_focal, &h_camera_focal, sizeof(point3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_world_light_dir, &world_light_dir, sizeof(point3), cudaMemcpyHostToDevice);



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

    Choose_closest <<< gridDim, blockDim >>> (d_distances, Face_NUM, d_colors, d_Planes, d_Materials, d_ray, d_Object_to_Face, d_close_indexes, d_closest_interesections, d_world_light_dir);
    cudaDeviceSynchronize();

    int* Reflected_surface = new int[WIDTH * HEIGHT];

    cudaMemcpy(Colors, d_colors, WIDTH * HEIGHT * 3 * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(Reflected_surface, d_close_indexes, WIDTH * HEIGHT * sizeof(int), cudaMemcpyDeviceToHost);

    saveAsBMP(Colors, WIDTH, HEIGHT, "normals_image.bmp");
    ///////////////////////////////////////////////
    //
    // Dodanie cienia
    //
    ////////////////////////////////////////////////


    Add_shadows <<< gridDim, blockDim >> > (d_closest_interesections, d_shadows, d_normal_index_to_face, d_number_of_vertices_in_one_face, d_Faces, d_Vertices, d_Normals, d_Planes, d_start_face_at_index, Face_NUM, Vert_NUM, Normal_NUM, d_world_light_dir);
    cudaDeviceSynchronize();

    cudaMemcpy(shadows, d_shadows, WIDTH * HEIGHT * 3 * sizeof(float), cudaMemcpyDeviceToHost);


    saveAsBMP(shadows, WIDTH, HEIGHT, "shadowes_image.bmp");




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

    Choose_closest << <gridDim, blockDim >> > (d_distances, Face_NUM, d_colors, d_Planes, d_Materials, d_ray, d_Object_to_Face, d_close_indexes, d_closest_interesections, d_world_light_dir);
    cudaDeviceSynchronize();


    float* Reflections = new float[WIDTH * HEIGHT * 3];

    cudaMemcpy(Reflections, d_colors, WIDTH * HEIGHT * 3 * sizeof(float), cudaMemcpyDeviceToHost);
    saveAsBMP(Reflections, WIDTH, HEIGHT, "reflected_image.bmp");


    cout << "skłądanie obrazów" << endl;

    for (int i = 0; i < WIDTH * HEIGHT; i++)
    {
        if (Reflected_surface[i] != -1) 
        {
            Material Mat;
            Mat = Materials[Object_to_Face[Reflected_surface[i]]];
            Reflections[i * 3] = Mat.Reflectivity * Reflections[i * 3] + Colors[i * 3] * (1 - Mat.Reflectivity);
            Reflections[i * 3 + 1] = Mat.Reflectivity * Reflections[i * 3 + 1] + Colors[i * 3 + 1] * (1 - Mat.Reflectivity);
            Reflections[i * 3 + 2] = Mat.Reflectivity * Reflections[i * 3 + 2] + Colors[i * 3 + 2] * (1 - Mat.Reflectivity);
            if (shadows[i * 3] >= 1) { Reflections[i * 3] -= 0.1f; if (Reflections[i * 3] < 0) { Reflections[i * 3] = 0.0f; } }
            if (shadows[i * 3 + 1] >= 1) { Reflections[i * 3 + 1] -= 0.1f; if (Reflections[i * 3 + 1] < 0) { Reflections[i * 3 + 1] = 0.0f; } }
            if (shadows[i * 3 + 2] >= 1) { Reflections[i * 3 + 2] -= 0.1f; if (Reflections[i * 3 + 2] < 0) { Reflections[i * 3 + 2] = 0.0f; } }
        
        }

    }
    saveAsBMP(Reflections, WIDTH, HEIGHT, "Complete_image.bmp");


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