// blender_to_cpp.cpp: definiuje punkt wejścia dla aplikacji.


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "vec3.h"
#include "ImportObj.h"

using namespace std;

int main() {

    map<int, vec3> vertices_coors;
    map<int, vector<int>> vertices_to_faces;
    map<int, vec3> faces_normals;
    vec3 test = vec3(1.1, 2.2, 3.3);
    cout << test << endl;


    void GetDataFromObj(vertices_coors, vertices_to_faces, faces_normals, "..\\..\\..\\helpers\\untitled.obj");


    return 0;
}
