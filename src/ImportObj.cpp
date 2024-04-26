#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "vec3.h"
#include "ImportObj.h"
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


void GetDataFromObj(map<int, vec3> &vertices_coors, map<int, vector<int>> &vertices_to_faces, map<int, vec3> &faces_normals, string file_path)
{
    ifstream plik;
    plik.open(file_path);
    if (!plik.is_open()) {
        cerr << "Nie uda³o siê otworzyæ pliku: " << file_path << endl;
        return;
    }
    string linia;

    vector<string> v;

    int index_v = 0;
    int index_f = 0;
    int index_normal = 0;
    bool UV_map = 0;

    while (getline(plik, linia))
    {
        vector<float> axis;
        vector<float> normal_indexes;
        vector<int> vertex_indexes;
        if (linia[0] == 'v' || linia[0] == 'f') {
            split(linia, v, ' ');
            for (const auto& element : v)
            {
                if (element == "v") { index_v++; }
                else if (element == "vn") { index_normal++; }
                else if (element == "f") { index_f++; }
                else if (element == "vt") { UV_map = 1; }
                else
                {

                    if (index_f == 0 && index_v != 0 && index_normal == 0)
                    {
                        axis.push_back(stof(element));
                        vertices_coors[index_v] = vec3((double)axis[0], (double)axis[1], (double)axis[2]);
                    }
                    if (index_normal != 0 && UV_map == 0)
                    {
                        normal_indexes.push_back(stof(element));
                        faces_normals[index_normal] = vec3((double)normal_indexes[0], (double)normal_indexes[1], (double)normal_indexes[2]);
                    }
                    if (index_f != 0)
                    {
                        vertex_indexes.push_back(stof(element));
                        vertices_to_faces[index_f] = vertex_indexes;
                    }
                }

            }

        }
    }
    plik.close();

    for (const auto& pair : vertices_coors) {
        cout << "Vertex: " << pair.first << ", Koordynaty: ";
        cout << pair.second << " ";
        cout << endl;
    }

    cout << endl;
    cout << endl;
    cout << endl;

    for (const auto& pair : faces_normals) {
        cout << "Vertex: " << pair.first << ", Normalna: ";
        cout << pair.second << " ";
        cout << endl;

    }

    cout << endl;
    cout << endl;
    cout << endl;

    for (const auto& pair : vertices_to_faces) {
        cout << "Face: " << pair.first << ", Lista vertex: ";
        for (int value : pair.second) {
            cout << value << " ";
        }
        cout << endl;
    }

}