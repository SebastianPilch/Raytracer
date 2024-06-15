#include "material.cuh"

__host__ __device__ Material::Material() {
    Diffuse[0] = 0.0f; Diffuse[1] = 0.0f; Diffuse[2] = 0.0f;
    Specular[0] = 0.0f; Specular[1] = 0.0f; Specular[2] = 0.0f;
    Ambient[0] = 0.0f; Ambient[1] = 0.0f; Ambient[2] = 0.0f;
    Alpha = 1.0f;
    Shininess = 0.0f;
}

__host__ __device__ Material::Material(float d0, float d1, float d2,
                                       float s0, float s1, float s2,
                                       float a0, float a1, float a2,
                                       float alpha, float shininess) {
    Diffuse[0] = clamp(d0, 0.0f, 1.0f);
    Diffuse[1] = clamp(d1, 0.0f, 1.0f);
    Diffuse[2] = clamp(d2, 0.0f, 1.0f);

    Specular[0] = clamp(s0, 0.0f, 1.0f);
    Specular[1] = clamp(s1, 0.0f, 1.0f);
    Specular[2] = clamp(s2, 0.0f, 1.0f);

    Ambient[0] = clamp(a0, 0.0f, 1.0f);
    Ambient[1] = clamp(a1, 0.0f, 1.0f);
    Ambient[2] = clamp(a2, 0.0f, 1.0f);

    Alpha = clamp(alpha, 0.0f, 1.0f);
    Shininess = clamp(shininess, 0.0f, 1.0f);
}

__host__ __device__ float Material::clamp(float value, float min, float max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}
