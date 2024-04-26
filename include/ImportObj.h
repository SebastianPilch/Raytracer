#ifndef IMPORTOBJ
#define IMPORTOBJ
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "vec3.h"
#include "ImportObj.h"

using namespace std;


size_t split(const string& txt, vector<string>& strs, char ch);

void GetDataFromObj(map<int, vec3> &vertices_coors, map<int, vector<int>> &vertices_to_faces, map<int, vec3> &faces_normals, string file_path);



#endif //ROWNANIAGAUSS_LINKEDLIST_HPP