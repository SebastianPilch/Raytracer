#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include "vec3.cuh"
#include "Parllel_fun.cuh"
#include "Ray.cuh"
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

__global__ void Choose_closest(float* d_distances,int Face_NUM, float* d_closest_normals, float* d_Planes, Material* Mats, ray* rays)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < WIDTH && j < HEIGHT) {
        int closest = -1;
        float lowest_distance = INFINITY;

        for (int f = 0; f < Face_NUM; f++)
        {
            if (d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] < lowest_distance)
            {
                closest = f;
                lowest_distance = d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f];
            }
        }

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
            float norm_length = sqrt(normal.x() * normal.x() + normal.y() * normal.y() + normal.z() * normal.z());
            normal = normal / norm_length;


            vec3 L = vec3(-5, -5, -5); // tu trzeba jakis wektor swiatla globalnego

            vec3 reflection_vector = 2 * (normal * L) * normal - L;

            color diffuse = vec3(Mats[1].Diffuse[0], Mats[1].Diffuse[1], Mats[1].Diffuse[2]);
            color specular = vec3(Mats[1].Specular[0], Mats[1].Specular[1], Mats[1].Specular[2]);
            color ambient = vec3(Mats[1].Ambient[0], Mats[1].Ambient[1], Mats[1].Ambient[2]);
            float shininess = Mats[1].Shininess;
            float alpha = Mats[1].Alpha;

            color light_ambient = vec3(0.1f, 0.1f, 0.1f);
            color light_diffuse = vec3(0.1f, 0.1f, 0.1f);
            color light_specular = vec3(0.2f, 0.2f, 0.2f);
            vec3 camera_vector = rays[j * WIDTH + i].dir;

            color face_color = ambient * light_ambient + diffuse * light_diffuse * max(dotProduct_(normal, L), 0.0f) + specular * light_specular * pow(max(dotProduct_(reflection_vector, camera_vector), 0.0f), shininess) * alpha;
            //face_color = diffuse * light_diffuse * max(dotProduct_(normal, L), 0.0f);
            //printf("Color: %f, %f, %f   \n", Mats[0].Diffuse[0], Mats[0].Diffuse[1], Mats[0].Diffuse[2]);

            d_closest_normals[(j * WIDTH + i) * 3 + 0] = face_color[0];
            d_closest_normals[(j * WIDTH + i) * 3 + 1] = face_color[1];
            d_closest_normals[(j * WIDTH + i) * 3 + 2] = face_color[2];
        }
        //printf("%f \n", d_distances[j * Face_NUM * WIDTH + i * Face_NUM + closest]);
    }
}

__device__ void matrixVectorMul(float* matrix, float* vector, float* result) {
    float x = vector[0];
    float y = vector[1];
    float z = vector[2];
    result[0] = matrix[0] * x + matrix[4] * y + matrix[8] * z + matrix[12];
    result[1] = matrix[1] * x + matrix[5] * y + matrix[9] * z + matrix[13];
    result[2] = matrix[2] * x + matrix[6] * y + matrix[10] * z + matrix[14];
}


__global__ void Transform(float* Vertexes, int numVertices, float TranslateX, float TranslateY, float TranslateZ, float rotateX, float rotateY, float rotateZ, float scaleX, float scaleY, float scaleZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < numVertices) 
    {
        float vertex[4] = { Vertexes[i * 3], Vertexes[i * 3 + 1], Vertexes[i * 3 + 2], 1.0f };
        float result[4];
        float transformationMatrix[16];

        transformationMatrix[0] = scaleX * (cos(rotateY)* cos(rotateZ));
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

__global__ void Update_normals_and_Planes(float* d_Vertices, int* d_Faces, float* d_Normals, float* d_Planes, int* d_number_of_vertices_in_one_face, int* d_normal_index_to_face, int* d_start_face_at_index, int Face_NUM, int Normals_NUM)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < Face_NUM) 
    {
        int modify_normal_idx = d_normal_index_to_face[i] - 1;

        int Vertex_1_idx = d_Faces[d_start_face_at_index[i]] - 1;
        int Vertex_2_idx = d_Faces[d_start_face_at_index[i]+1] - 1;
        int Vertex_3_idx = d_Faces[d_start_face_at_index[i]+2] - 1;

        point3 Vertex_1 = point3(d_Vertices[3* Vertex_1_idx], d_Vertices[3 * Vertex_1_idx + 1], d_Vertices[3 * Vertex_1_idx + 2]);
        point3 Vertex_2 = point3(d_Vertices[3 * Vertex_2_idx], d_Vertices[3 * Vertex_2_idx + 1], d_Vertices[3 * Vertex_2_idx + 2]);
        point3 Vertex_3 = point3(d_Vertices[3 * Vertex_3_idx], d_Vertices[3 * Vertex_3_idx + 1], d_Vertices[3 * Vertex_3_idx + 2]);

        vec3 Vec_12 = vec3(Vertex_2[0] - Vertex_1[0], Vertex_2[1] - Vertex_1[1], Vertex_2[2] - Vertex_1[2]);
        vec3 Vec_13 = vec3(Vertex_3[0] - Vertex_1[0], Vertex_3[1] - Vertex_1[1], Vertex_3[2] - Vertex_1[2]);

        vec3 new_Norm =  crossProduct_(Vec_12, Vec_13);
        Plane new_Plane = Plane(new_Norm,Vertex_1);

        d_Planes[i * 4] = new_Plane.A;
        d_Planes[i * 4+1] = new_Plane.B;
        d_Planes[i * 4+2] = new_Plane.C;
        d_Planes[i * 4+3] = new_Plane.D;

        d_Normals[modify_normal_idx * 3] = new_Norm[0];
        d_Normals[modify_normal_idx * 3 + 1] = new_Norm[1];
        d_Normals[modify_normal_idx * 3 + 2] = new_Norm[2];


    }
}