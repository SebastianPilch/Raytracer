#include <iostream>
#include <fstream>
#include "ImportObj.cuh"
//#include <cuda_runtime.h>
//#include <device_launch_parameters.h>

using namespace std;

size_t split(const std::string& txt, std::vector<std::string>& strs, char ch)
{
    size_t pos = txt.find(ch);
    size_t initialPos = 0;
    strs.clear();

    // Decompose statement
    while (pos != std::string::npos) {
        strs.push_back(txt.substr(initialPos, pos - initialPos));
        initialPos = pos + 1;

        pos = txt.find(ch, initialPos);
    }

    // Add the last one
    strs.push_back(txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1));

    return strs.size();
}


Pointer_storage GetDataFromObj(int& Vertices_coords_size, int& Face_numer,int& Normals_size, string file_path)
{
    ifstream plik(file_path);
    if (!plik.is_open()) 
    {
        cerr << "Nie uda�o si� otworzy� pliku: " << file_path << endl;
        return Pointer_storage();
    }
    
    string linia;
    vector<string> Splited_line;
    vector<string> Backshlash_face_split;
    size_t newSize = 1;
    float** vertices = new float*[Vertices_coords_size];
    for (size_t i = 0; i < Vertices_coords_size; ++i) {
        vertices[i] = new float[3];
        for (size_t j = 0; j < 3; ++j) { vertices[i][j] = 0; }
    }
    float** faces_normals = new float* [Normals_size];
    for (size_t i = 0; i < Normals_size; ++i) {
        faces_normals[i] = new float[3];
    }
    float** Planes_to_faces = new float* [Face_numer];
    for (size_t i = 0; i < Face_numer; ++i) {
        Planes_to_faces[i] = new float[3];
    }
    int** vertices_to_faces = new int* [Face_numer];
    for (size_t i = 0; i < Face_numer; ++i) {
        vertices_to_faces[i] = new int[3];
    }

    int*  noramls_index_to_face = new int[Face_numer];
    int* vertices_in_one_face = new int[Face_numer];


    float** newVertices = nullptr;
    float** newNormals = nullptr;
    float** newFaces = nullptr;
    float** newPlanes = nullptr;

    int* new_face_lengths = nullptr;
    int* new_normals_index_to_face = nullptr;


    int index_v = 0;
    int index_f = 0;
    int index_n = 0;

    while (getline(plik, linia)) 
    {
        split(linia, Splited_line, ' ');
        if (Splited_line[0] == "v") {
            if (index_v >= Vertices_coords_size){
                size_t newSize = Vertices_coords_size * 2;
                float** newVertices = new float* [newSize];
                for (size_t i = 0; i < newSize; ++i) {
                    newVertices[i] = new float[3];
                }
                for (size_t i = 0; i < index_v; ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        newVertices[i][j] = vertices[i][j];
                    }
                    delete[] vertices[i];
                }
                delete[] vertices;
                vertices = newVertices;
                for (size_t i = 0; i < newSize;i++) {
                    vertices[i] = newVertices[i];
                }
                Vertices_coords_size = newSize;
            }
            vertices[index_v] = new float[3];
            vertices[index_v][0] = (float)stof(Splited_line[1]);
            vertices[index_v][1] = (float)stof(Splited_line[2]);
            vertices[index_v][2] = (float)stof(Splited_line[3]);
            index_v++;
        }


        //cout << "v - posz�o" << endl;

        if (Splited_line[0] == "vn")
        {
            if (index_n >= Normals_size) 
            {
                size_t newSize = Normals_size * 2;
                float** newNormals = new float* [newSize];
                for (size_t i = 0; i < newSize; ++i) {
                    newNormals[i] = new float[3];
                }
                // Przekopiowanie danych
                for (size_t i = 0; i < index_n; ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        newNormals[i][j] = faces_normals[i][j];
                    }
                    delete[] faces_normals[i];
                }
                delete[] faces_normals;
                faces_normals = newNormals;
                for (size_t i = 0; i < newSize; ++i) {
                    faces_normals[i] = newNormals[i];
                }
                Normals_size = newSize;

            }
            faces_normals[index_n] = new float[3];
            faces_normals[index_n][0] = stof(Splited_line[1]);
            faces_normals[index_n][1] = stof(Splited_line[2]);
            faces_normals[index_n][2] = stof(Splited_line[3]);
            index_n++;
        }




        if (Splited_line[0] == "f")
        {
            if (index_f >= Face_numer) {
                size_t newSize = Face_numer * 2;
                
                // alokacja nowego wektora przypisuj�cego ile wierzcho�k�w tworzy i-t� �cian�
                int* new_face_lengths = new int[newSize];

                for (size_t i = 0; i < index_f; i++) {
                       new_face_lengths[i] = vertices_in_one_face[i];
                }
                delete[] vertices_in_one_face;
                vertices_in_one_face = new_face_lengths;
                // koniec alokacji



                //alokacja nowego wektora przypisuj�cego do indeksu �ciany indeks jej nomalnej

                int* new_normals_index_to_face = new int[newSize];
                for (size_t i = 0; i < index_f; i++) {
                    new_normals_index_to_face[i] = noramls_index_to_face[i];
                }
                delete[] noramls_index_to_face;
                noramls_index_to_face = new_normals_index_to_face;

                //koniec alokacji


                //alokacja nowego wektora p�aszczyzn zawieraj�cych �ciany
                float** NewPlanes = new float* [newSize];
                for (size_t i = 0; i < index_f; i++) {
                    NewPlanes[i] = new float[4];
                }
                for (size_t i = 0; i < index_f; i++) {
                    for (size_t j = 0; j < 4; ++j) {
                        NewPlanes[i][j] = Planes_to_faces[i][j];
                    }
                    delete[] Planes_to_faces[i];
                }
                delete[] Planes_to_faces;
                Planes_to_faces = NewPlanes;
                for (size_t i = 0; i < newSize; i++) {
                    Planes_to_faces[i] = NewPlanes[i];
                }
                //koniec alokacji



                // alokacja nowego wektora przypisuj�ca do indeksu �ciany indeksy wierzcho�k�w kt�re j� twrz�
                int** newFaces = new int* [newSize];
                for (size_t i = 0; i < index_f; i++) {
                    newFaces[i] = new int[vertices_in_one_face[i]];
                }
                for (size_t i = 0; i < index_f; i++) {
                    for (size_t j = 0; j < vertices_in_one_face[i]; ++j) {
                        newFaces[i][j] = vertices_to_faces[i][j];
                    }
                    delete[] vertices_to_faces[i];
                }
                delete[] vertices_to_faces;
                vertices_to_faces = newFaces;
                for (size_t i = 0; i < newSize; i++) {
                    vertices_to_faces[i] = newFaces[i];
                }
                // koniec alokacji
                Face_numer = newSize;
            }
            // przypisanie ilo�ci wierzcho�k�w do �ciany i stworzenie odpowiadaj�cego tej d�ugo�ci wektora w li�cie �cian:
            vertices_in_one_face[index_f] = Splited_line.size() - 1;
            vertices_to_faces[index_f] = new int[vertices_in_one_face[index_f]];
            //split po  "/" w obr�bie danych przypisanych do �ciany  f wierzcho�ek/UV wierzcho�ka/normalna �ciany
            for(int i = 0;i < vertices_in_one_face[index_f];i++)
            {
                vertices_to_faces[index_f][i] = (int)stof(Splited_line[i+1]);
            }
            split(Splited_line[1], Backshlash_face_split, '/');
            noramls_index_to_face[index_f] = (int)stof(Backshlash_face_split[2]);
            
            int current_face_normal = noramls_index_to_face[index_f];
            float x = vertices[vertices_to_faces[index_f][0]][0];
            float y = vertices[vertices_to_faces[index_f][0]][1];
            float z = vertices[vertices_to_faces[index_f][0]][2];
            float nor_x = faces_normals[current_face_normal-1][0];
            float nor_y = faces_normals[current_face_normal-1][1];
            float nor_z = faces_normals[current_face_normal-1][2];
            Planes_to_faces[index_f] = new float[4];
            Planes_to_faces[index_f][0] = nor_x;
            Planes_to_faces[index_f][1] = nor_y;
            Planes_to_faces[index_f][2] = nor_z;
            Planes_to_faces[index_f][3] = -(nor_x * x + nor_y * y + nor_z * z);
            index_f++;
        }
    }

    Vertices_coords_size = index_v;
    Face_numer = index_f;
    Normals_size = index_n;


    plik.close();

    delete[] newVertices  ;
    delete[] newNormals ;
    delete[] newFaces;
    delete[] newPlanes;
    delete[] new_face_lengths;
    delete[] new_normals_index_to_face;
    
    return Pointer_storage(vertices, faces_normals, vertices_to_faces, Planes_to_faces, vertices_in_one_face, noramls_index_to_face);
}


void Print_Import_data(Pointer_storage object,int ver_size,int nor_size,int Face_size) 
{
    for (int i = 0; i < ver_size; i++)
    {
        cout << "Vertex " << i + 1 << " :  " << object.Vertices[i][0] << " ,  " << object.Vertices[i][1] << " ,  "
            << object.Vertices[i][2] << endl;
    }

    cout << endl << endl;


    for (int i = 0; i < nor_size; i++)
    {
        cout << "Normal " << i + 1 << " :  " << object.Normals[i][0] << " ,  " << object.Normals[i][1] << " ,  "
            << object.Normals[i][2] << endl;
    }
    cout << endl << endl;


    for (int i = 0; i < Face_size; i++)
    {
        cout << "Vertices creating face " << i + 1 << " :  " << object.Face_size[i] << endl;
    }

    cout << endl << endl;

    for (int i = 0; i < Face_size; i++)
    {
        cout << "Normal index to face " << i + 1 << " :  " << object.Face_to_Normal[i]  << endl;
    }
    cout << endl << endl;


    for (int i = 0; i < Face_size; i++)
    {
        cout << "Plane " << i + 1 << " :  ";
        for (int j = 0; j < 4; j++)
        {
            cout << object.Planes[i][j] << " ,  ";
        }
        cout << endl;
    }


    cout << endl << endl;

    for (int i = 0; i < Face_size; i++)
    {
        cout << "Face " << i + 1 << " :  ";
        for (int j = 0; j < object.Face_size[i]; j++) 
        {
            cout << object.Faces[i][j] << " ,  ";
        }
        cout  << endl;
    }

}