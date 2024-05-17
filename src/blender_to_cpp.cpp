// blender_to_cpp.cpp: definiuje punkt wejścia dla aplikacji.


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "vec3.h"
#include "ImportObj.h"
#include "Ray.h"

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

    map<int, point3> vertices_coors;
    map<int, vector<int>> vertices_to_faces;
    map<int, vec3> faces_normals;
    map<int, Plane> Planes_to_faces;
    GetDataFromObj(vertices_coors, vertices_to_faces, faces_normals,Planes_to_faces, "..\\..\\..\\helpers\\untitled.obj");

    int img_height = 400;
    int img_width = 400;

    auto focal_length = 1.0;
    auto viewport_height = 2.0;
    auto viewport_width = viewport_height * (double(img_width) / img_height);
    auto camera_center = point3(0, 0, 5);
    auto viewport_u = vec3(viewport_width, 0, 0);
    auto viewport_v = vec3(0, -viewport_height, 0);

    auto pixel_delta_u = viewport_u / img_width;
    auto pixel_delta_v = viewport_v / img_height;
    auto viewport_upper_left = camera_center - vec3(0, 0, focal_length) - viewport_u / 2 - viewport_v / 2;
    auto pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);


    std::vector<std::vector<float>> img(img_height, std::vector<float>(img_width));

    for (int j = 0; j < img_height; j++) {
        for (int i = 0; i < img_width; i++) {
            auto pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
            auto ray_direction = pixel_center - camera_center;
            ray newRay = ray(camera_center, ray_direction);

            float hit_anything = 0.0;
            for (const auto& pair : Planes_to_faces) {
                bool hit = Face_hit(pair.second, newRay, vertices_to_faces[pair.first], vertices_coors);
                if (hit) {
                    hit_anything += 0.2;
                }
                
            }

            img[j][i] = hit_anything;
        }
    }

    int width = img[0].size();
    int height = img.size();

    saveAsBMP(img, width, height, "result_image.bmp");

    //ray mainray = ray(point3(0, 0.1,0), vec3(2, 0, 2));

    //for (const auto& pair : planes_to_faces) {
    //    bool hit = face_hit(pair.second, mainray, vertices_to_faces[pair.first], vertices_coors);
    //    if (hit) { cout << "nastapilo odbicie od sciany "<< pair.first << endl; }
    //    else { cout << "nie nastapilo odbicie od sciany " << pair.first << endl; }    
    //}
 

    return 0;
}