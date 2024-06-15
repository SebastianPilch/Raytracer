#ifndef MATERIAL_CUH
#define MATERIAL_CUH

class Material {
public:
    float Diffuse[3];
    float Specular[3];
    float Ambient[3];
    float Alpha;
    float Shininess;

    __host__ __device__ Material();

    __host__ __device__ Material(const float diffuse[3], const float specular[3], const float ambient[3], float alpha, float shininess);

private:
    __host__ __device__ float clamp(float value, float min, float max);
};

#endif 
