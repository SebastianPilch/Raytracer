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


__global__ void Generate_rays(ray* viewport_rays, double focal_length, point3* camera_center, point3* camera_focal)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= WIDTH || j >= HEIGHT || z != 0) {
        return;
    }

    float aspect_ratio = float(WIDTH) / float(HEIGHT);
    float viewport_height = 2.0;
    float viewport_width = aspect_ratio * viewport_height;
    float focal_dist = focal_length;

    vec3 w = (*camera_center - *camera_focal) / ((*camera_center - *camera_focal).length());
    vec3 u = (crossProduct_(vec3(0, 1, 0), w)) / ((crossProduct_(vec3(0, 1, 0), w)).length());
    vec3 v = crossProduct_(w, u);

    vec3 horizontal = u * viewport_width;
    vec3 vertical = v * viewport_height;
    vec3 lower_left_corner = *camera_center - horizontal / 2 - vertical / 2 - w * focal_dist;

    float u_offset = float(i) / (WIDTH - 1);
    float v_offset = float(j) / (HEIGHT - 1);
    vec3 pixel_center = lower_left_corner + horizontal * u_offset + vertical * v_offset;

    vec3 ray_direction = pixel_center - *camera_center;

    viewport_rays[j * WIDTH + i] = ray(*camera_center, ray_direction);

    __syncthreads();

}

__global__ void Update_rays(ray* viewport_rays, float* First_Intersections, int* Intersected_face_idx, int* d_normal_index_to_face, float* d_Normals)
{    
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        int j = blockIdx.y * blockDim.y + threadIdx.y;
        int z = blockIdx.z * blockDim.z + threadIdx.z;

        if (i >= WIDTH || j >= HEIGHT || z != 0) 
        {
            return;
        } 


        vec3 UV_ray_dir = viewport_rays[j * WIDTH + i].dir;
        int f = Intersected_face_idx[j * WIDTH + i];
        if (f == -1) 
        {
            return;
        }
        int normal_idx = d_normal_index_to_face[f];
        vec3 normal = vec3(d_Normals[3 * normal_idx], d_Normals[3 * normal_idx + 1], d_Normals[3 * normal_idx + 2]);
        normal = normal / normal.length();
        UV_ray_dir = UV_ray_dir / UV_ray_dir.length();
        vec3 ray_direction = UV_ray_dir - 2 * dotProduct_(normal, UV_ray_dir) * normal ;

        point3 pixel_center = vec3(First_Intersections[(j * WIDTH + i) * 3], First_Intersections[(j * WIDTH + i) * 3 + 1], First_Intersections[(j * WIDTH + i) * 3 + 2]);

        viewport_rays[j * WIDTH + i] = ray(pixel_center, ray_direction);

       //printf("%f , %f , %f \n", First_Intersections[j * WIDTH + i], First_Intersections[j * WIDTH + i+1], First_Intersections[j * WIDTH + i+2]);


        __syncthreads();

}



__global__ void Generate_distances(ray* viewport_rays,point3* camera_center, float* d_intersections, int* d_normal_index_to_face, int* d_number_of_vertices_in_one_face,
    int* d_Faces, float* d_Vertices, float* d_Normals, float* d_Planes, int* start_face_at_index,
    int Face_NUM, int Vertex_NUM, int Normal_NUM, float* d_distances, int reflections)
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int f = blockIdx.z * blockDim.z + threadIdx.z;

    if (f >= Face_NUM || i >= WIDTH || j >= HEIGHT) 
    {
        return;
    }
    bool Is_Hitten_correct = true;

    ray UV_ray = viewport_rays[j * HEIGHT + i];

    Plane plane = Plane((double)d_Planes[f * 4], (double)d_Planes[f * 4 + 1], (double)d_Planes[f * 4 + 2], (double)d_Planes[f * 4 + 3]);
    inter_ inter_data = UV_ray.findIntersection(plane);
    point3 intersection = inter_data.intersection;
    double t = inter_data.t;


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
    vec3 dis;
    if (reflections == 0) 
    {
        dis = *camera_center - intersection;
    }
    else 
    {
        point3 last_inter = vec3(d_intersections[(j * WIDTH + i) * 3], d_intersections[(j * WIDTH + i) * 3 + 1], d_intersections[(j * WIDTH + i)*3 + 2]);
        dis = last_inter - intersection;

    }

    if (intersection[0] == INFINITY || intersection[0] == -INFINITY || intersection[1] == INFINITY ||
        intersection[1] == -INFINITY || intersection[2] == INFINITY || intersection[2] == -INFINITY || t < 0)
    {
        Is_Hitten_correct = false;
    }



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

__global__ void Choose_closest(float* d_distances,int Face_NUM, float* d_colors, float* d_Planes, Material* Mats, ray* rays, int* Mats_to_face, int* Intersected_face_idx, float* d_intersections)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    

    if (i < WIDTH && j < HEIGHT && z == 1)
    {
        ray UV_ray = rays[j * WIDTH + i];


        color light_ambient = vec3(0.2f, 0.2f, 0.2f);
        color light_diffuse = vec3(0.7f, 0.7f, 0.7f);
        color light_specular = vec3(1.0f, 1.0f, 1.0f);
        vec3 camera_vector = rays[j * WIDTH + i].dir;

        vec3 L = vec3(-10, 12, 15); // tu trzeba jakis wektor swiatla globalnego
        //L = L / L.length();


        int closest = -1;
        float lowest_distance = INFINITY;
        float far = -1.0f;

        for (int f = 0; f < Face_NUM; f++)
        {
            if (d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f] < lowest_distance && f != Intersected_face_idx[(j * WIDTH + i)])
            {
                closest = f;
                lowest_distance = d_distances[j * Face_NUM * WIDTH + i * Face_NUM + f];
            }


        }
        if (lowest_distance == INFINITY || closest == -1)
        {
            d_colors[(j * WIDTH + i) * 3 + 0] = 0.5f;
            d_colors[(j * WIDTH + i) * 3 + 1] = 0.5f;
            d_colors[(j * WIDTH + i) * 3 + 2] = 0.5f;

            //printf("Indeks: %d, a f: %d\n", j * WIDTH + i, z);

            d_intersections[(j * WIDTH + i)*3] = INFINITY;
            d_intersections[(j * WIDTH + i) * 3 + 1] = INFINITY;
            d_intersections[(j * WIDTH + i) * 3 + 2] = INFINITY;

            return;
        }
        Plane plane = Plane((double)d_Planes[closest * 4], (double)d_Planes[closest * 4 + 1], (double)d_Planes[closest * 4 + 2], (double)d_Planes[closest * 4 + 3]);
        inter_ inter_data  = UV_ray.findIntersection(plane);
        point3 intersection = inter_data.intersection;
        double t = inter_data.t;



        d_intersections[(j * WIDTH + i) * 3] = intersection[0];
        d_intersections[(j * WIDTH + i) * 3 + 1] = intersection[1];
        d_intersections[(j * WIDTH + i) * 3 + 2] = intersection[2];
        //printf("%f , %f , %f \n", d_intersections[j * WIDTH + i], d_intersections[j * WIDTH + i+1], d_intersections[j * WIDTH + i+2]);


        Intersected_face_idx[j * WIDTH + i] = closest;

        int Mat_idx = Mats_to_face[closest];

        color diffuse = vec3(Mats[Mat_idx].Diffuse[0], Mats[Mat_idx].Diffuse[1], Mats[Mat_idx].Diffuse[2]);
        color specular = vec3(Mats[Mat_idx].Specular[0], Mats[Mat_idx].Specular[1], Mats[Mat_idx].Specular[2]);
        color ambient = vec3(Mats[Mat_idx].Ambient[0], Mats[Mat_idx].Ambient[1], Mats[Mat_idx].Ambient[2]);
        float shininess = Mats[Mat_idx].Shininess;
        float alpha = Mats[Mat_idx].Alpha;



            vec3 normal = vec3((float)plane.A, (float)plane.B, (float)plane.C);
            normal = - normal / normal.length();

            L = L / L.length();
            vec3 reflection_vector = 2 * dotProduct_(normal, L) * normal - L;
            reflection_vector = reflection_vector / reflection_vector.length();
            camera_vector = camera_vector / camera_vector.length();

            color face_color = {
                alpha * (
                    ambient[0] * light_ambient[0] +
                    diffuse[0] * light_diffuse[0] * max(0.0f, dotProduct_(normal, L)) +
                    specular[0] * light_specular[0] * pow(max(0.0f, dotProduct_(reflection_vector, camera_vector)), shininess)
                ),

                alpha * (
                    ambient[1] * light_ambient[1] +
                    diffuse[1] * light_diffuse[1] * max(0.0f, dotProduct_(normal, L)) +
                    specular[1] * light_specular[1] * pow(max(0.0f, dotProduct_(reflection_vector, camera_vector)), shininess)
                ),

                alpha * (
                    ambient[2] * light_ambient[2] +
                    diffuse[2] * light_diffuse[2] * max(0.0f, dotProduct_(normal, L)) +
                    specular[2] * light_specular[2] * pow(max(0.0f, dotProduct_(reflection_vector, camera_vector)), shininess)
                )
            };

            face_color[0] = min(1.0f, max(0.0f, face_color[0]));
            face_color[1] = min(1.0f, max(0.0f, face_color[1]));
            face_color[2] = min(1.0f, max(0.0f, face_color[2]));

                //printf("%f\n", face_color[0]);
            d_colors[(j * WIDTH + i) * 3 + 0] = face_color[0];
            d_colors[(j * WIDTH + i) * 3 + 1] = face_color[1];
            d_colors[(j * WIDTH + i) * 3 + 2] = face_color[2];



    }
}




__global__ void Add_shadows(float* d_intersections, float* d_shadows, int* d_normal_index_to_face, int* d_number_of_vertices_in_one_face, int* d_Faces, float* d_Vertices, float* d_Normals, float* d_Planes, int* start_face_at_index, int Face_NUM, int Vertex_NUM, int Normal_NUM)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int f = blockIdx.z * blockDim.z + threadIdx.z;



    if (i < WIDTH && j < HEIGHT && f < Face_NUM)
    {


        //if (d_intersections[(j * WIDTH + i) * 3] == INFINITY || d_intersections[(j * WIDTH + i) * 3 + 1] == INFINITY || d_intersections[(j * WIDTH + i) * 3 + 2] == INFINITY)
        //{
        //    return;
        //}
        ray Reverse_Light = ray(point3(d_intersections[(j * WIDTH + i) * 3], d_intersections[(j * WIDTH + i) * 3 + 1], d_intersections[(j * WIDTH + i) * 3 + 2]), vec3(-1,-1, -1));

        bool Is_Hitten_correct = true;
        Plane plane = Plane((double)d_Planes[f * 4], (double)d_Planes[f * 4 + 1], (double)d_Planes[f * 4 + 2], (double)d_Planes[f * 4 + 3]);
        inter_ inter_data = Reverse_Light.findIntersection(plane);
        point3 intersection = inter_data.intersection;
        double t = inter_data.t;

        if (intersection[0] == INFINITY || intersection[0] == -INFINITY || intersection[1] == INFINITY ||
            intersection[1] == -INFINITY || intersection[2] == INFINITY || intersection[2] == -INFINITY)
        {
            Is_Hitten_correct = false;
        }
        //printf("%f , %f , %f \n", d_intersections[(j * WIDTH + i) * 3], d_intersections[(j * WIDTH + i) * 3 + 1], d_intersections[(j * WIDTH + i) * 3 + 2]);

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
        }


        if(Is_Hitten_correct && t > 0 )
        {
            d_shadows[(j * WIDTH + i)*3] = 1;
            d_shadows[(j * WIDTH + i)*3 + 1] = 1;
            d_shadows[(j * WIDTH + i)*3 + 2] = 1;

        }


        __syncthreads();
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

