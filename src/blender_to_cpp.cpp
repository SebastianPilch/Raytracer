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
    GetDataFromObj(vertices_coors, vertices_to_faces, faces_normals, "..\\..\\..\\helpers\\untitled.obj");


    return 0;
}


 //sprawdzanie płaszczyzny i punktu przecięcia

//#include <iostream>
//#include <cmath>
//#include <vector>
//struct Vector3 {
//    double x, y, z;
//
//    Vector3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
//    Vector3 operator-(const Vector3& other) const {
//        return { x - other.x, y - other.y, z - other.z };
//    }
//};

//// Struktura reprezentująca płaszczyznę w postaci ogólnej
//struct Plane {
//    double A, B, C, D;
//
//    Plane(double _A, double _B, double _C, double _D) : A(_A), B(_B), C(_C), D(_D) {}
//};
//
//// Struktura reprezentująca prostą w postaci parametrycznej
//struct Line {
//    Vector3 p; // Punkt na prostej
//    Vector3 d; // Wektor kierunkowy prostej
//
//    Line(Vector3 _p, Vector3 _d) : p(_p), d(_d) {}
//};
//
//// Funkcja obliczająca punkt wspólny prostej i płaszczyzny
//Vector3 findIntersection(const Line& line, const Plane& plane) {
//    double t = -(plane.A * line.p.x + plane.B * line.p.y + plane.C * line.p.z + plane.D) /
//        (plane.A * line.d.x + plane.B * line.d.y + plane.C * line.d.z);
//
//    Vector3 intersection(line.p.x + t * line.d.x,
//        line.p.y + t * line.d.y,
//        line.p.z + t * line.d.z);
//
//    return intersection;
//}
//
//int main() {
//    // Przykładowe dane - płaszczyzna i prosta
//    Plane plane(0, -1, 0, 1); 
//    Plane plane2(0, 1, 0, 1);
//
//    Line line(Vector3(0, 0, 0), Vector3(1, 1, 1)); 
//
//    Vector3 intersection = findIntersection(line, plane);
//
//    // Wyświetlenie współrzędnych punktu wspólnego
//    std::cout << "Punkt wspolny prostej i plaszczyzny: (" << intersection.x << ", "
//        << intersection.y << ", " << intersection.z << ")" << std::endl;
//    intersection = findIntersection(line, plane2);
//
//    // Wyświetlenie współrzędnych punktu wspólnego
//    std::cout << "Punkt wspolny prostej i plaszczyzny: (" << intersection.x << ", "
//        << intersection.y << ", " << intersection.z << ")" << std::endl;
//    return 0;
//}



// Punkt wewnątrz Poligonu


//Vector3 crossProduct(const Vector3& a, const Vector3& b) {
//    Vector3 result(0,0,0);
//    result.x = a.y * b.z - a.z * b.y;
//    result.y = a.z * b.x - a.x * b.z;
//    result.z = a.x * b.y - a.y * b.x;
//    return result;
//}
//
//// Function to calculate dot product of two vectors
//double dotProduct(const Vector3& a, const Vector3& b) {
//    return a.x * b.x + a.y * b.y + a.z * b.z;
//}
//
//bool isPointInsidePolygon(const Vector3& point, const std::vector<Vector3>& polygon) {
//    Vector3 normal = crossProduct(polygon[1] - polygon[0], polygon[2] - polygon[0]);
//
//    // Check if the point is inside the polygon using ray casting algorithm
//    Vector3 edge(0,0,0);
//    for (size_t i = 0; i < polygon.size(); ++i) {
//        edge = polygon[(i + 1) % polygon.size()] - polygon[i];
//        Vector3 vp = point - polygon[i];
//        Vector3 n = crossProduct(edge, vp);
//        if (dotProduct(n, normal) < 0) {
//            return false;
//        }
//    }
//    return true;
//}
//
//int main() 
//{
//    std::vector<Vector3> Poligon;
//    Poligon = { Vector3(2,1,0), Vector3(0,1,0),Vector3(0,1,2),Vector3(2,1,2) };
//    Vector3 point(1,1,1);
//
//    std::cout << isPointInsidePolygon(point, Poligon) << std::endl;
//}