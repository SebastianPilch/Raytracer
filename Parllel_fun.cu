#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include "vec3.cuh"
#include "Parllel_fun.cuh"
#include "Ray.cuh"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>

using namespace std;
#ifndef __CUDACC__ 
#define __CUDACC__
#endif

__global__ void Generate_rays(ray* viewport_rays, double focal_length, point3* camera_center,
    point3* camera_focal, int* d_normal_index_to_face, int* d_number_of_vertices_in_one_face,
    int* d_Faces, float* d_Vertices, float* d_Normals, float* d_Planes, int* start_face_at_index,
    int Face_NUM, int Vertex_NUM, int Normal_NUM, float* d_distances, float* d_closest_angles) {



    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int f = blockIdx.z * blockDim.z + threadIdx.z;

    if (f >= Face_NUM || i >= WIDTH || j >= HEIGHT) return;


    float aspect_ratio = float(WIDTH) / float(HEIGHT);
    float viewport_height = 2.0;
    float viewport_width = aspect_ratio * viewport_height;
    float focal_dist = focal_length;

    vec3 w = (*camera_center - *camera_focal) / ((*camera_center - *camera_focal).length());
    vec3 u = (crossProduct_(vec3(0, 1, 0), w)) / ((crossProduct_(vec3(0, 1, 0), w)).length());
    vec3 v = crossProduct_(w,u);

    vec3 horizontal = u * viewport_width;
    vec3 vertical = v * viewport_height;
    vec3 lower_left_corner = *camera_center - horizontal / 2 - vertical / 2 - w * focal_dist;

    float u_offset = float(i) / (WIDTH - 1);
    float v_offset = float(j) / (HEIGHT - 1);
    vec3 pixel_center = lower_left_corner + horizontal * u_offset + vertical * v_offset;

/*    vec3 viewport_upper_left = *camera_center - vec3(focal_length, focal_length, focal_length) - VIEWPORT_U / 2 - VIEWPORT_V / 2;
    vec3 pixel00_loc = viewport_upper_left + 0.5 * (DELTA_U + DELTA_V);

    vec3 pixel_center = pixel00_loc + (i * DELTA_U) + (j * DELTA_V) + *camera_focal;*/


    vec3 ray_direction = pixel_center - *camera_center;
    ray UV_ray = ray(*camera_center, ray_direction);
    viewport_rays[j * HEIGHT + i] = ray(*camera_center, ray_direction);

    __syncthreads();

    Plane plane = Plane((double)d_Planes[f * 4], (double)d_Planes[f * 4 + 1], (double)d_Planes[f * 4 + 2], (double)d_Planes[f * 4 + 3]);
    point3 intersection = UV_ray.findIntersection(plane);
    bool Is_Hitten_correct = true;
    if (intersection[0] == INFINITY || intersection[0] == -INFINITY || intersection[1] == INFINITY ||
        intersection[1] == -INFINITY || intersection[2] == INFINITY || intersection[2] == -INFINITY)
    {
        Is_Hitten_correct = false;
    }

    vec3 edge;
    size_t next_index;
    for (size_t idx = start_face_at_index[f]; idx < start_face_at_index[f] + d_number_of_vertices_in_one_face[f]; ++idx)
    {
        size_t next_index;
        if (idx + 1 >= start_face_at_index[f] + d_number_of_vertices_in_one_face[f])
        {
            next_index = start_face_at_index[f];
        }
        else
        {
            next_index = idx + 1;
        }
        int vertex_index = d_Faces[idx] - 1;
        int vertex_next = d_Faces[next_index] - 1;

        point3 current_vertex = point3(
            (double)d_Vertices[3 * vertex_index],
            (double)d_Vertices[3 * vertex_index + 1],
            (double)d_Vertices[3 * vertex_index + 2]
        );

        point3 next_vertex = point3(
            (double)d_Vertices[3 * vertex_next],
            (double)d_Vertices[3 * vertex_next + 1],
            (double)d_Vertices[3 * vertex_next + 2]
        );

        __syncthreads();
        edge = next_vertex - current_vertex;
        vec3 vp = intersection - current_vertex;
        vec3 n = crossProduct_(edge, vp);
        vec3 normal = vec3(plane.A, plane.B, plane.C);
        if (dotProduct_(n, normal) < 0) 
        {
            Is_Hitten_correct = false;
        }
        if (current_vertex[0] == 0 && current_vertex[1] == 0 && current_vertex[2] == 0)
        {
            //printf("%d", vertex_index);
        }
    }
    vec3 dis = *camera_center - intersection;

    if (!Is_Hitten_correct)
    {
        d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] = INFINITY;
        //printf("%f  %f  %f \n", intersection[0], intersection[1], intersection[2]);

    }
    else
    {
        float distance = sqrt(dis[0] * dis[0] + dis[1] * dis[1] + dis[2] * dis[2]);
        d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] = distance;

    }


    __syncthreads();
}

__global__ void Choose_closest(float* d_distances,int Face_NUM, float* d_closest_normals, float* d_Planes, Material* Mats, ray* rays, int* Mats_to_face)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;


    if (i < WIDTH && j < HEIGHT) 
    {



        color light_ambient = vec3(0.2f, 0.2f, 0.2f);
        color light_diffuse = vec3(0.7f, 0.7f, 0.7f);
        color light_specular = vec3(1.0f, 1.0f, 1.0f);
        vec3 camera_vector = rays[j * WIDTH + i].dir;

        vec3 L = vec3(-15, -15, -15); // tu trzeba jakis wektor swiatla globalnego
        //L = L / L.length();


        int closest = -1;
        float lowest_distance = INFINITY;
        float far = -1.0f;

        for (int f = 0; f < Face_NUM; f++)
        {
            if (d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] < lowest_distance)
            {
                closest = f;
                lowest_distance = d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f];
            }
            if (d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] != INFINITY && d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] > far)
            {
                far = d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f];
            }

        }

        int Mat_idx = Mats_to_face[closest] ;

        color diffuse = vec3(Mats[Mat_idx].Diffuse[0], Mats[Mat_idx].Diffuse[1], Mats[Mat_idx].Diffuse[2]);
        color specular = vec3(Mats[Mat_idx].Specular[0], Mats[Mat_idx].Specular[1], Mats[Mat_idx].Specular[2]);
        color ambient = vec3(Mats[Mat_idx].Ambient[0], Mats[Mat_idx].Ambient[1], Mats[Mat_idx].Ambient[2]);
        float shininess = Mats[Mat_idx].Shininess;
        float alpha = Mats[Mat_idx].Alpha;

        Plane plane = Plane((double)d_Planes[closest * 4], (double)d_Planes[closest * 4 + 1], (double)d_Planes[closest * 4 + 2], (double)d_Planes[closest * 4 + 3]);

        if (d_distances[j * Face_NUM * WIDTH + i * Face_NUM + closest] == INFINITY || d_distances[j * Face_NUM * WIDTH + i * Face_NUM + closest] < 0.0f)
        {
            d_closest_normals[(j * WIDTH + i) * 3 + 0] = 0.0f;
            d_closest_normals[(j * WIDTH + i) * 3 + 1] = 0.0f;
            d_closest_normals[(j * WIDTH + i) * 3 + 2] = 0.0f;
        }
        else
        {



            vec3 normal = vec3((float)plane.A, (float)plane.B, (float)plane.C);
            //color face_color;
            float norm_length = sqrt(normal.x() * normal.x() + normal.y() * normal.y() + normal.z() * normal.z());
            normal = normal / norm_length;
            vec3 reflection_vector = 2  * dotProduct_(normal, L) * normal - L;
            color face_color = {
                ambient[0] * light_ambient[0] +
                alpha * (
                    diffuse[0] * light_diffuse[0] * max(0.0f, dotProduct_(normal, L) / L.length()) +
                    specular[0] * light_specular[0] * pow(max(0.0f, -dotProduct_(reflection_vector, camera_vector) / (reflection_vector.length() * camera_vector.length())), shininess)
                ),

                ambient[1] * light_ambient[1] +
                alpha * (
                    diffuse[1] * light_diffuse[1] * max(0.0f, dotProduct_(normal, L) / L.length()) +
                    specular[1] * light_specular[1] * pow(max(0.0f, -dotProduct_(reflection_vector, camera_vector) / (reflection_vector.length() * camera_vector.length())), shininess)
                ),

                ambient[2] * light_ambient[2] +
                alpha * (
                    diffuse[2] * light_diffuse[2] * max(0.0f, dotProduct_(normal, L) / L.length()) +
                    specular[2] * light_specular[2] * pow(max(0.0f, -dotProduct_(reflection_vector, camera_vector) / (reflection_vector.length() * camera_vector.length())), shininess)
                )
            };


            d_closest_normals[(j * WIDTH + i) * 3 + 0] = face_color[0];
            d_closest_normals[(j * WIDTH + i) * 3 + 1] = face_color[1];
            d_closest_normals[(j * WIDTH + i) * 3 + 2] = face_color[2];

        }
    }
}

__device__ void matrixVectorMul(float* matrix, float* vector, float* result) 
{
    float x = vector[0];
    float y = vector[1];
    float z = vector[2];
    result[0] = matrix[0] * x + matrix[4] * y + matrix[8] * z + matrix[12];
    result[1] = matrix[1] * x + matrix[5] * y + matrix[9] * z + matrix[13];
    result[2] = matrix[2] * x + matrix[6] * y + matrix[10] * z + matrix[14];
}


__global__ void Transform(float* Vertexes, int numVertices, int* Object_to_Vertex, int index, float TranslateX, float TranslateY, float TranslateZ, float rotateX, float rotateY, float rotateZ, float scaleX, float scaleY, float scaleZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < numVertices) 
    {
        if (Object_to_Vertex[i] == index) 
        {
            float vertex[4] = { Vertexes[i * 3], Vertexes[i * 3 + 1], Vertexes[i * 3 + 2], 1.0f };
            float result[4];
            float transformationMatrix[16];

            transformationMatrix[0] = scaleX * (cos(rotateY) * cos(rotateZ));
            transformationMatrix[1] = scaleY * (cos(rotateZ) * sin(rotateX) * sin(rotateY) - cos(rotateX) * sin(rotateZ));
            transformationMatrix[2] = scaleZ * (cos(rotateX) * cos(rotateZ) * sin(rotateY) + sin(rotateX) * sin(rotateZ));
            transformationMatrix[3] = TranslateX;
            transformationMatrix[4] = scaleX * (cos(rotateY) * sin(rotateZ));
            transformationMatrix[5] = scaleY * (cos(rotateX) * cos(rotateZ) + sin(rotateX) * sin(rotateY) * sin(rotateZ));
            transformationMatrix[6] = scaleZ * (cos(rotateX) * sin(rotateY) * sin(rotateZ) - cos(rotateZ) * sin(rotateX));
            transformationMatrix[7] = TranslateY;
            transformationMatrix[8] = scaleX * (-sin(rotateY));
            transformationMatrix[9] = scaleY * (cos(rotateY) * sin(rotateX));
            transformationMatrix[10] = scaleZ * (cos(rotateX) * cos(rotateY));
            transformationMatrix[11] = TranslateZ;
            transformationMatrix[12] = 0.0f;
            transformationMatrix[13] = 0.0f;
            transformationMatrix[14] = 0.0f;
            transformationMatrix[15] = 1.0f;

            matrixVectorMul(transformationMatrix, vertex, result);

            Vertexes[i * 3] = result[0] + TranslateX;
            Vertexes[i * 3 + 1] = result[1] + TranslateY;
            Vertexes[i * 3 + 2] = result[2] + TranslateZ;
        }
    }
}

__global__ void Update_normals_and_Planes(float* d_Vertices, int* d_Faces, int* Object_to_Face, int index, float* d_Normals, float* d_Planes, int* d_number_of_vertices_in_one_face, int* d_normal_index_to_face, int* d_start_face_at_index, int Face_NUM, int Normals_NUM)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < Face_NUM) 
    {
            int modify_normal_idx = d_normal_index_to_face[i] - 1;

            int Vertex_1_idx = d_Faces[d_start_face_at_index[i]] - 1;
            int Vertex_2_idx = d_Faces[d_start_face_at_index[i] + 1] - 1;
            int Vertex_3_idx = d_Faces[d_start_face_at_index[i] + 2] - 1;

            point3 Vertex_1 = point3(d_Vertices[3 * Vertex_1_idx], d_Vertices[3 * Vertex_1_idx + 1], d_Vertices[3 * Vertex_1_idx + 2]);
            point3 Vertex_2 = point3(d_Vertices[3 * Vertex_2_idx], d_Vertices[3 * Vertex_2_idx + 1], d_Vertices[3 * Vertex_2_idx + 2]);
            point3 Vertex_3 = point3(d_Vertices[3 * Vertex_3_idx], d_Vertices[3 * Vertex_3_idx + 1], d_Vertices[3 * Vertex_3_idx + 2]);

            vec3 Vec_12 = vec3(Vertex_2[0] - Vertex_1[0], Vertex_2[1] - Vertex_1[1], Vertex_2[2] - Vertex_1[2]);
            vec3 Vec_13 = vec3(Vertex_3[0] - Vertex_1[0], Vertex_3[1] - Vertex_1[1], Vertex_3[2] - Vertex_1[2]);

            vec3 new_Norm = crossProduct_(Vec_12, Vec_13);
            Plane new_Plane = Plane(new_Norm, Vertex_1);

            d_Planes[i * 4] = new_Plane.A;
            d_Planes[i * 4 + 1] = new_Plane.B;
            d_Planes[i * 4 + 2] = new_Plane.C;
            d_Planes[i * 4 + 3] = new_Plane.D;

            d_Normals[modify_normal_idx * 3] = new_Norm[0];
            d_Normals[modify_normal_idx * 3 + 1] = new_Norm[1];
            d_Normals[modify_normal_idx * 3 + 2] = new_Norm[2];
    }
}



__device__ int partition(float* data, int left, int right) {
    float pivot = data[right];
    int i = left - 1;
    for (int j = left; j < right; ++j) {
        if (data[j] < pivot) {
            ++i;
            float temp = data[i];
            data[i] = data[j];
            data[j] = temp;
        }
    }
    float temp = data[i + 1];
    data[i + 1] = data[right];
    data[right] = temp;
    return i + 1;
}

__device__ void quickSort(float* data, int left, int right) {
    if (left < right) {
        int pi = partition(data, left, right);

        quickSort(data, left, pi - 1);
        quickSort(data, pi + 1, right);
    }
}

__global__ void quickSortKernel(float* d_distances, int Face_NUM) {
    int i_ind = blockIdx.x * blockDim.x + threadIdx.x;
    int j_ind = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i_ind < WIDTH && j_ind < HEIGHT) 
    {
        int start = j_ind * Face_NUM * WIDTH + i_ind * Face_NUM;
        int end = start + Face_NUM - 1;
        quickSort(d_distances, start, end);
    }
}

