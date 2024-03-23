// blender_to_cpp.cpp: definiuje punkt wejścia dla aplikacji.


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

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


int main() {
    map<int, vector<float>> vertices_coors;
    map<int, vector<int>> vertices_to_faces;
    map<int, vector<float>> faces_normals;
//    std::string nazwaPliku = "teks.txt";
    ifstream plik;
    plik.open("C:\\Users\\Sebastian\\Desktop\\blender\\blender_to_cpp\\blender_to_cpp\\untitled.obj");
    //plik.open("C:\\Users\\Sebastian\\Desktop\\blender\\blender_to_cpp\\blender_to_cpp\\twarz_bez_refki.obj");
    if (!plik.is_open()) {
        cerr << "Nie udało się otworzyć pliku: " << "teks.txt" << endl;
        return 1;
    }
    string linia;

    std::vector<std::string> v;

    int index_v = 0;
    int index_f = 0;
    int index_normal = 0;
    bool UV_map = 0;

    while (getline(plik, linia))
    {
        vector<float> axis;
        vector<float> normal_indexes;
        vector<int> vertex_indexes;
        if (linia[0] == 'v' or linia[0] == 'f') {
            split(linia, v, ' ');
            for (const auto& element : v)
            {
                if (element == "v") { index_v++; }
                else if (element == "vn") { index_normal++; }
                else if (element == "f") { index_f++; }
                else if (element == "vt") { UV_map = 1; }
                else
                {

                    if (index_f == 0 and index_v != 0 and index_normal == 0)
                    {
                        axis.push_back(stof(element));
                        vertices_coors[index_v] = axis;
                    }
                    if (index_normal != 0 and UV_map == 0)
                    {
                        normal_indexes.push_back(stof(element));
                        faces_normals[index_normal] = normal_indexes;
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
        std::cout << "Vertex: " << pair.first << ", Koordynaty: ";
        for (float value : pair.second) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    for (const auto& pair : faces_normals) {
        std::cout << "Vertex: " << pair.first << ", Normalna: ";
        for (float value : pair.second) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    for (const auto& pair : vertices_to_faces) {
        std::cout << "Face: " << pair.first << ", Lista vertex: ";

        
        for (int value : pair.second) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
