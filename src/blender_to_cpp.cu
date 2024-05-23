// blender_to_cpp.cpp: definiuje punkt wejścia dla aplikacji.


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "vec3.cuh"
#include "ImportObj.cuh"
#include "Ray.cuh"
#include "Parllel_fun.cuh"

using namespace std;

void saveAsBMP(const std::vector<std::vector<float>>& img, int width, int height, const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);

    if (!file) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    int paddingSize = (4 - (width * 3) % 4) % 4; // Padding wymagany przez format BMP

    // Nagłówek BMP
    int filesize = 54 + (3 * width + paddingSize) * height;
    char fileHeader[54] = { 'B', 'M', 0,0,0,0, 0,0, 0,0, 54,0,0,0, 40,0,0,0, static_cast<char>(width), static_cast<char>(width >> 8), static_cast<char>(width >> 16), static_cast<char>(width >> 24), static_cast<char>(height), static_cast<char>(height >> 8), static_cast<char>(height >> 16), static_cast<char>(height >> 24), 1,0, 24,0, 0,0,0,0, static_cast<char>(filesize), static_cast<char>(filesize >> 8), static_cast<char>(filesize >> 16), static_cast<char>(filesize >> 24), 0,0,0,0, 0,0,0,0 };

    // Zapisanie nagłówka
    file.write(fileHeader, 54);

    // Zapisanie danych pikseli
    for (int i = height - 1; i >= 0; i--) {
        for (int j = 0; j < width; j++) {
            unsigned char color = static_cast<unsigned char>(img[i][j] * 255); // Skalowanie wartości z [0, 1] do [0, 255]
            file.put(color);
            file.put(color);
            file.put(color);
        }
        // Dodanie paddingu
        for (int k = 0; k < paddingSize; k++) {
            file.put(0);
        }
    }

    file.close();
}



int main() {

    int vert_num = 3;
    float** vert = new float* [vert_num];

    for (int i = 0; i < vert_num; ++i) {
        vert[i] = new float[3];  // Allocating 3 floats for each row
    }

    map<int, vector<int>> vertices_to_faces;
    map<int, vec3> faces_normals;
    map<int, Plane> Planes_to_faces;
    GetDataFromObj(vert, vert_num,  vertices_to_faces, faces_normals,Planes_to_faces, "..\\..\\..\\helpers\\Trojkatny_szescian.obj");



    //saveAsBMP(img, width, height, "result_image.bmp");





    //float** x = new float* [4];
    //float** y = new float* [4];
    //float** z= new float* [4];



    //MatrixAddition(4,4,x,y,z);


    //for (int i = 0; i < 4; i++) {
    //    for (int j = 0; j < 4; j++) 
    //    {
    //        cout <<" " << z[i][j];

    //    }
    //    cout << endl;
    //}












    return 0;
}