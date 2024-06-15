#include "vec3.cuh"
#include "Ray.cuh"

#ifndef PAR_FUN_H
#define PAR_FUN_H


#define WIDTH 400
#define HEIGHT 400
#define VIEWPORT_HEIGHT  2.0
#define VIEWPORT_WIDTH VIEWPORT_HEIGHT * (double)(WIDTH / HEIGHT)
#define VIEWPORT_U vec3(VIEWPORT_WIDTH, 0, 0)
#define VIEWPORT_V vec3(0, -VIEWPORT_HEIGHT, 0)
#define DELTA_U VIEWPORT_U / WIDTH
#define DELTA_V VIEWPORT_V / HEIGHT


__global__ void Generate_rays(ray* viewport_rays, double focal_length, point3* camera_center,
    point3* camera_focal, int* d_normal_index_to_face, int* d_number_of_vertices_in_one_face,
    int* d_Faces, float* d_Vertices, float* d_Normals, float* d_Planes, int* start_face_at_index,
    int Face_NUM, int Vertex_NUM, int Normal_NUM, float* d_distances, float* d_closest_normals);


__global__ void Choose_closest(float* d_distances, int Face_NUM, float* d_closest_normals, float* d_Planes);


__device__ void matrixVectorMul(float* matrix, float* vector, float* result);

__global__ void Transform(float* Vertexes, int numVertices, float TranslateX, float TranslateY, float TranslateZ, float rotateX, float rotateY, float rotateZ, float scaleX,float scaleY,float scaleZ);


__global__ void Update_normals_and_Planes(float* d_Vertices, int* d_Faces, float* d_Normals, float* d_Planes, int* d_number_of_vertices_in_one_face, int* d_normal_index_to_face, int* d_start_face_at_index, int Face_NUM, int Normals_NUM);





#endif // PAR_FUN_H