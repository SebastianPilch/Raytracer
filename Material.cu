#include "material.cuh"

__host__ __device__ Material::Material() {
    Diffuse[0] = 0.0f; Diffuse[1] = 0.0f; Diffuse[2] = 0.0f;
    Specular[0] = 0.0f; Specular[1] = 0.0f; Specular[2] = 0.0f;
    Ambient[0] = 0.0f; Ambient[1] = 0.0f; Ambient[2] = 0.0f;
    Alpha = 1.0f;
    Shininess = 0.0f;
}

__host__ __device__ Material::Material(const float diffuse[3], const float specular[3], const float ambient[3], float alpha, float shininess) {
    for (int i = 0; i < 3; i++) {
        Diffuse[i] = clamp(diffuse[i], 0.0f, 1.0f);
        Specular[i] = clamp(specular[i], 0.0f, 1.0f);
        Ambient[i] = clamp(ambient[i], 0.0f, 1.0f);
    }
    Alpha = clamp(alpha, 0.0f, 1.0f);
    Shininess = clamp(shininess, 0.0f, 1.0f);
}

__host__ __device__ float Material::clamp(float value, float min, float max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}
