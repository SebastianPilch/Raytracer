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
#include<cmath>

void saveAsBMP(ray** img, int width, int height, const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);

    if (!file) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    int paddingSize = (4 - (width * 3) % 4) % 4; // Padding required by BMP format

    // BMP header
    int filesize = 54 + (3 * width + paddingSize) * height;
    char fileHeader[54] = { 'B', 'M', 0,0,0,0, 0,0, 0,0, 54,0,0,0, 40,0,0,0, static_cast<char>(width), static_cast<char>(width >> 8), static_cast<char>(width >> 16), static_cast<char>(width >> 24), static_cast<char>(height), static_cast<char>(height >> 8), static_cast<char>(height >> 16), static_cast<char>(height >> 24), 1,0, 24,0, 0,0,0,0, static_cast<char>(filesize), static_cast<char>(filesize >> 8), static_cast<char>(filesize >> 16), static_cast<char>(filesize >> 24), 0,0,0,0, 0,0,0,0 };

    // Write header
    file.write(fileHeader, 54);

    // Write pixel data
    for (int i = height - 1; i >= 0; i--) {
        for (int j = 0; j < width; j++) {
            float focal_len = 5.0f;
            unsigned char color = static_cast<unsigned char>(focal_len / (sqrt(img[i][j].dir.x() * img[i][j].dir.x() + img[i][j].dir.y() * img[i][j].dir.y())) / 4); // Scale value from [0, 1] to [0, 255]
            file.put(color);
            file.put(color);
            file.put(color);
        }

        // Add padding
        for (int k = 0; k < paddingSize; k++) {
            file.put(0);
        }
    }

    file.close();
}


void saveAsBMP2(float* angles, int width, int height, const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);

    if (!file) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    int paddingSize = (4 - (width * 3) % 4) % 4; // Padding required by BMP format

    // BMP header
    int filesize = 54 + (3 * width + paddingSize) * height;
    char fileHeader[54] = {
        'B', 'M',                         // Signature
        static_cast<char>(filesize), static_cast<char>(filesize >> 8), static_cast<char>(filesize >> 16), static_cast<char>(filesize >> 24), // File size
        0,0, 0,0,                         // Reserved
        54,0,0,0,                         // File offset to pixel array
        40,0,0,0,                         // DIB header size
        static_cast<char>(width), static_cast<char>(width >> 8), static_cast<char>(width >> 16), static_cast<char>(width >> 24), // Width
        static_cast<char>(height), static_cast<char>(height >> 8), static_cast<char>(height >> 16), static_cast<char>(height >> 24), // Height
        1,0,                              // Planes
        24,0,                             // Bits per pixel
        0,0,0,0,                          // Compression
        0,0,0,0,                          // Image size (can be 0 for uncompressed)
        0,0,0,0,                          // X pixels per meter (unused)
        0,0,0,0,                          // Y pixels per meter (unused)
        0,0,0,0,                          // Total colors (0 means default)
        0,0,0,0                           // Important colors (0 means all are important)
    };

    // Write header
    file.write(fileHeader, 54);

    // Write pixel data
    for (int i = height - 1; i >= 0; i--) {
        for (int j = 0; j < width; j++) {
            // Retrieve RGB values from angles array
            float red = angles[(i * width + j) * 3 + 0];
            float green = angles[(i * width + j) * 3 + 1];
            float blue = angles[(i * width + j) * 3 + 2];

            unsigned char r = static_cast<unsigned char>(red * 255.99f);
            unsigned char g = static_cast<unsigned char>(green * 255.99f);
            unsigned char b = static_cast<unsigned char>(blue * 255.99f);

            file.put(b); // Blue channel
            file.put(g); // Green channel
            file.put(r); // Red channel
        }

        // Add padding
        for (int k = 0; k < paddingSize; k++) {
            file.put(0);
        }
    }

    file.close();
}


int main() {

    ///////////////////////////////////////////////
    //
    //wczytanie obiektów .obj
    //
    //
    //
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
    Pointer_storage pociety_walec = GetDataFromObj(vert_num3, face_num3, normal_num3, object_couter, "../../../helpers/complete_scene.obj");

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
    //
    //
    //
    //////////////////////////////////////////////////
    int Vert_NUM = vert_num3;
    int Face_NUM = face_num3;
    int Normal_NUM = normal_num3;
    Pointer_storage liczony_objekt = pociety_walec;
    ///////////////////////////////////////////////
    //
    //przepisanie wskaźników
    //
    //
    //
    //
    //////////////////////////////////////////////////
    float** Planes = new float* [Face_NUM];
    Planes[0] = new float[Face_NUM * 4];

    float** Verticies = new float* [Vert_NUM];
    Verticies[0] = new float[Vert_NUM * 3 ];

    float** Normals = new float* [Normal_NUM];
    Normals[0] = new float[Normal_NUM * 3];

    float* Distances = new float[WIDTH * HEIGHT * Face_NUM];
    float* ClosestNormals = new float[WIDTH * HEIGHT * Face_NUM * 3];

    int* number_of_vertices_in_one_face = liczony_objekt.Face_size;
    int* normal_index_to_face = liczony_objekt.Face_to_Normal;

    int* start_face_at_index = new int[Face_NUM];

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
        start_face_at_index[i-1] = current_index;
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
            Planes[i][j] = liczony_objekt.Planes[i][j];
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
    //
    //
    //
    //////////////////////////////////////////////////
    Material* Materials = new Material[3];
    Materials[0] = Material( 0.5f, 0.6f, 0.7f , 0.3f, 0.4f, 0.5f , 0.1f, 0.2f, 0.3f , 1.0f, 0.0f);
    Materials[1] = Material(0.7f, 0.2f, 0.2f, 0.3f, 0.4f, 0.5f, 0.1f, 0.2f, 0.3f, 1.0f, 5.0f);
    Materials[2] = Material(0.5f, 0.6f, 0.7f, 0.3f, 0.4f, 0.5f, 0.1f, 0.2f, 0.3f, 0.7f, 0.0f);
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
    float* d_distances;
    float* d_Vertices;
    float* d_Normals;
    float* d_Planes;
    float* d_closest_normals;
    Material* d_Materials;

    cudaMalloc(&d_Faces, Length_to_Allocate_Faces * sizeof(int));
    cudaMalloc(&d_Planes, 4 * Face_NUM * sizeof(float));
    cudaMalloc(&d_Normals, 3 * Normal_NUM * sizeof(float));
    cudaMalloc(&d_Vertices, 3 * Vert_NUM * sizeof(float));
    cudaMalloc(&d_number_of_vertices_in_one_face, Face_NUM * sizeof(int));
    cudaMalloc(&d_normal_index_to_face, Face_NUM * sizeof(int));
    cudaMalloc(&d_start_face_at_index, Face_NUM * sizeof(int));
    cudaMalloc(&d_distances, WIDTH * HEIGHT * Face_NUM * sizeof(float));
    cudaMalloc(&d_closest_normals, WIDTH * HEIGHT * Face_NUM * 3 * sizeof(float));

    cudaMemcpy(d_Faces, Faces[0], Length_to_Allocate_Faces * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Planes, Planes[0], 4 * Face_NUM * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Normals, Normals[0], 3 * Face_NUM * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Vertices, Verticies[0], 3 * Vert_NUM * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_number_of_vertices_in_one_face, number_of_vertices_in_one_face, Face_NUM * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_normal_index_to_face, normal_index_to_face, Face_NUM * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_start_face_at_index, start_face_at_index, Face_NUM * sizeof(int), cudaMemcpyHostToDevice);

    ///////////////////////////////////////////////
    //
    //Tranformacja obiektu rotacja/skala/przesunięcie
    //
    //
    //
    //
    //////////////////////////////////////////////////

    int threadsPerBlock = 256;
    int blocksPerGrid = (Vert_NUM + threadsPerBlock - 1) / threadsPerBlock;
    float TranslateX = 0.0f;
    float TranslateY = 0.0f;
    float TranslateZ = 0.0f;
    float rotateX = 0.0f;
    float rotateY = 23.0f;
    float rotateZ = 23.0f;
    float scaleX = 1.0f;
    float scaleY = 1.0f;
    float scaleZ = 1.0f;


    Transform<<<blocksPerGrid, threadsPerBlock >>>(d_Vertices, Vert_NUM,  TranslateX,  TranslateY,  TranslateZ,  rotateX,  rotateY,  rotateZ,  scaleX,  scaleY,  scaleZ);
    //cudaMemcpy(Verticies, d_Vertices, Vert_NUM * 3 * sizeof(float), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    threadsPerBlock = 256;
    blocksPerGrid = (Face_NUM + threadsPerBlock - 1) / threadsPerBlock;
    Update_normals_and_Planes << <blocksPerGrid, threadsPerBlock >> > (d_Vertices, d_Faces, d_Normals, d_Planes, d_number_of_vertices_in_one_face, d_normal_index_to_face, d_start_face_at_index,Face_NUM,Normal_NUM);
    cudaDeviceSynchronize();


    ///////////////////////////////////////////////
    //
    //  Wyznaczanie Promieni, uderzenia i dystanse
    //
    //
    //
    //
    //////////////////////////////////////////////////

    double focal_length = 10;
    point3 h_camera_center(45.0, 45.0, 45.0);
    point3 h_camera_focal(1.0,1.0, 1.0);
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

    Generate_rays << <gridDim, blockDim >> > (d_ray, focal_length, d_camera_center, d_camera_focal, d_normal_index_to_face, d_number_of_vertices_in_one_face,
        d_Faces, d_Vertices, d_Normals, d_Planes, d_start_face_at_index, Face_NUM, Vert_NUM, Normal_NUM, d_distances, d_closest_normals);
    cudaDeviceSynchronize();

    cout << endl << "sort start" << endl;

    //quickSortKernel << <gridDim, blockDim >> > (d_distances, Face_NUM);

    cudaMemcpy(h_ray[0], d_ray, WIDTH * HEIGHT * sizeof(ray), cudaMemcpyDeviceToHost);
    cudaMemcpy(Distances, d_distances, WIDTH * HEIGHT * Face_NUM * sizeof(float), cudaMemcpyDeviceToHost);


    ///////////////////////////////////////////////
    //
    // dobór najbliższej ściany z wektora dystansów
    //
    //
    //
    //
    //////////////////////////////////////////////////

    const int BLOCK_SIZE_X = 16;
    const int BLOCK_SIZE_Y = 16;
    dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
    int gridSizeX = (WIDTH + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X;
    int gridSizeY = (HEIGHT + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y;
    dim3 dimGrid(gridSizeX, gridSizeY);

    cudaMalloc(&d_Materials,3 * sizeof(Material));
    cudaMemcpy(d_distances, Distances, WIDTH * HEIGHT * Face_NUM * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Materials, Materials, 3 * sizeof(Material), cudaMemcpyHostToDevice);


    Choose_closest <<< gridDim, blockDim >>> (d_distances, Face_NUM,d_closest_normals,d_Planes,d_Materials,d_ray);
    cudaDeviceSynchronize();

    cudaMemcpy(ClosestNormals, d_closest_normals, WIDTH * HEIGHT * Face_NUM * 3 * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_ray);
    cudaFree(d_distances);
    cudaFree(d_camera_center);
    cudaFree(d_camera_focal);
    cudaFree(d_closest_normals);

    saveAsBMP(h_ray, WIDTH, HEIGHT, "result_image.bmp");
    saveAsBMP2(ClosestNormals, WIDTH, HEIGHT, "normals_image.bmp");

    free(h_ray[0]);
    free(h_ray);
    delete[] ClosestNormals;

    return 0;
}

