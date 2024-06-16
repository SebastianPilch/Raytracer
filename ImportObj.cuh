
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


struct Pointer_storage 
{
	float** Vertices;
	float** Normals;
	int** Faces;
	float** Planes;
	int* Face_size;
	int* Face_to_Normal;

	Pointer_storage(float**  V, float**  N, int** F, float** P,int* F_size, int* F2N)
	{
		this->Vertices = V;
		this->Normals = N;
		this->Faces = F;
		this->Planes = P;
		this->Face_size = F_size;
		this->Face_to_Normal = F2N;
	}	
	Pointer_storage() 
	{
		this->Vertices = nullptr;
		this->Normals = nullptr;
		this->Faces = nullptr;
		this->Planes = nullptr;
		this->Face_size = nullptr;
		this->Face_to_Normal = nullptr;
	}

};



size_t split(const string& txt, vector<string>& strs, char ch);

Pointer_storage GetDataFromObj(int& Vertices_coords_size, int& Face_numer, int& Normals_size, int& object_counter, string file_path);

void Print_Import_data(Pointer_storage object, int ver_size, int nor_size, int Face_size);

#endif 