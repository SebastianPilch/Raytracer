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



int main() {

    map<int, point3> vertices_coors;
    map<int, vector<int>> vertices_to_faces;
    map<int, vec3> faces_normals;
    map<int, Plane> Planes_to_faces;
    GetDataFromObj(vertices_coors, vertices_to_faces, faces_normals,Planes_to_faces, "..\\..\\..\\helpers\\untitled.obj");

    ray MainRay = ray(point3(0, 0.1,0), vec3(2, 0, 2));

    for (const auto& pair : Planes_to_faces) {
        bool hit = Face_hit(pair.second, MainRay, vertices_to_faces[pair.first], vertices_coors);
        if (hit) { cout << "Nastapilo odbicie od sciany "<< pair.first << endl; }
        else { cout << "Nie nastapilo odbicie od sciany " << pair.first << endl; }    
    }
    







    return 0;
}