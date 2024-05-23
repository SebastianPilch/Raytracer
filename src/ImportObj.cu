#include <iostream>
#include <fstream>
#include "ImportObj.cuh"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

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


void GetDataFromObj(float*& vertices_coors, int& Vertices_coords_size, map<int, vector<int>>& vertices_to_faces, map<int, vec3>& faces_normals, map<int, Plane>& Planes_to_faces, string file_path)
{
    ifstream plik;
    plik.open(file_path);
    if (!plik.is_open()) {
        cerr << "Nie uda³o siê otworzyæ pliku: " << file_path << endl;
        return;
    }
    string linia;

    vector<string> Splited_line;

    int index_v = 0;
    int index_f = 0;
    int index_normal = 0;
    bool UV_map = 0;


    while (getline(plik, linia)) {
        split(linia, Splited_line, ' ');
        if (Splited_line[0] == "v") 
        {
            if (index_v >= Vertices_coords_size) 
            {
                size_t newSize = Vertices_coords_size * 2;
                float** newVertices = new float* [newSize];
                for (size_t i = 0; i < newSize; ++i) 
                {
                    newVertices[i] = new float[3];
                }
                for (size_t i = 0; i < index_v; ++i) 
                {
                    newVertices[i][0] = vertices_coors[i][0];
                    newVertices[i][1] = vertices_coors[i][1];
                    newVertices[i][2] = vertices_coors[i][2];
                    delete[] vertices_coors[i]; // Usuniêcie starej pamiêci
                }
                delete[] vertices_coors; // Usuniêcie starej tablicy wskaŸników
                vertices_coors = newVertices;
                Vertices_coords_size = newSize;
            }

            vertices_coors[index_v * 3] = (float)stof(Splited_line[1]);
            vertices_coors[index_v * 3 + 1] = (float)stof(Splited_line[2]);
            vertices_coors[index_v * 3 + 2] = (float)stof(Splited_line[3]);
            index_v++;
        }
    




            if (Splited_line[0] == "vn")
            {

            }
            if (Splited_line[0] == "f")
            {

            }
            //for (const string& element : v)

            //        if (index_f == 0 && index_v != 0 && index_normal == 0)
            //        {
            //            axis.push_back(stof(element));
            //            vertices_coors[index_v] = point3((double)axis[0], (double)axis[1], (double)axis[2]);
            //        }
            //        if (index_normal != 0 && UV_map == 0)
            //        {
            //            normal_indexes.push_back(stof(element));
            //            faces_normals[index_normal] = vec3((double)normal_indexes[0], (double)normal_indexes[1], (double)normal_indexes[2]);
            //        }
            //        if (index_f != 0)
            //        {
            //            vertex_indexes.push_back(stof(element));
            //            vertices_to_faces[index_f] = vertex_indexes;
            //        }
            //    }

            //}

        Vertices_coords_size = index_v;
    }
    plik.close();

    for (int i = 0; i < Vertices_coords_size; i++) {
        cout << "Vertex: " << i + 1 << ", Coordinates: ";
        cout << vertices_coors[Vertices_coords_size][0] << ",  " << vertices_coors[Vertices_coords_size][1] << ",  " << vertices_coors[Vertices_coords_size][2];
        cout << endl;
    }

    cout << endl;
    cout << endl;
    cout << endl;

    //for (const auto& pair : faces_normals) {
    //    cout << "Face: " << pair.first << ", Normal: ";
    //    cout << pair.second << " ";
    //    Planes_to_faces[pair.first] = Plane(pair.second, vertices_coors[vertices_to_faces[pair.first][0]]);
    //    cout << endl;
    //}

    //cout << endl;
    //cout << endl;
    //cout << endl;

    //for (const auto& pair : vertices_to_faces) {
    //    cout << "Face: " << pair.first << ", Vertex list: ";
    //    for (int value : pair.second) {
    //        cout << value << " ";
    //    }
    //    cout << endl;
    //}
    //cout << endl;
    //cout << endl;
    //cout << endl;

    //for (const auto& pair : Planes_to_faces) {
    //    cout << "Face: " << pair.first << ", Plane: " << pair.second << " " << endl;
}
