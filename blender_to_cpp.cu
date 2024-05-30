
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
//    // Nagłówek BMP
//    int filesize = 54 + (3 * width + paddingSize) * height;
//    char fileHeader[54] = { 'B', 'M', 0,0,0,0, 0,0, 0,0, 54,0,0,0, 40,0,0,0, static_cast<char>(width), static_cast<char>(width >> 8), static_cast<char>(width >> 16), static_cast<char>(width >> 24), static_cast<char>(height), static_cast<char>(height >> 8), static_cast<char>(height >> 16), static_cast<char>(height >> 24), 1,0, 24,0, 0,0,0,0, static_cast<char>(filesize), static_cast<char>(filesize >> 8), static_cast<char>(filesize >> 16), static_cast<char>(filesize >> 24), 0,0,0,0, 0,0,0,0 };
//
//    // Zapisanie nagłówka
//    file.write(fileHeader, 54);
//
//    // Zapisanie danych pikseli
//    for (int i = height - 1; i >= 0; i--) {
//        for (int j = 0; j < width; j++) {
//            unsigned char color = static_cast<unsigned char>(img[i][j] * 255); // Skalowanie wartości z [0, 1] do [0, 255]
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

    return 0;

}


