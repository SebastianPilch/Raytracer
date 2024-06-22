#ifndef MATERIAL_CUH
#define MATERIAL_CUH

class Material {
public:
    float Diffuse[3];
    float Specular[3];
    float Ambient[3];
    float Alpha;
    float Shininess;
    float Reflectivity;

    __host__ __device__ Material();

    __host__ __device__ Material(float d0, float d1, float d2,
                                 float s0, float s1, float s2,
                                 float a0, float a1, float a2,
                                 float alpha, float shininess, 
                                 float reflectivity);

private:
    __host__ __device__ float clamp(float value, float min, float max);
};

#endif 
