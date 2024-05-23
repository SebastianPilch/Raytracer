#ifndef IMPORTOBJ
#define IMPORTOBJ
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "vec3.cuh"

//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
using namespace std;


size_t split(const string& txt, vector<string>& strs, char ch);

void GetDataFromObj(float*& vertices_coors, int& Vertices_coords_size, map<int, vector<int>> &vertices_to_faces, map<int, vec3> &faces_normals, map<int, Plane>& Planes_to_faces, string file_path);


#endif 