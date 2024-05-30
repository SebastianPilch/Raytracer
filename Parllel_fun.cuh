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


__global__ void Generate_rays(ray* viewport_rays, double focal_length, point3 *camera_center,
	point3 *camera_focal, int* d_normal_index_to_face, int* d_number_of_vertices_in_one_face,int* d_Faces,
	float* d_Vertices,float* d_Normals,float* d_Planes);

#endif // PAR_FUN_H