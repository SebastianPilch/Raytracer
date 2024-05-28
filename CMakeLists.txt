cmake_minimum_required(VERSION 3.18)
project(ProjektCUDA LANGUAGES CXX CUDA)

# Ustawienie wersji standardów
set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

# Dodanie ścieżki do plików nagłówkowych
include_directories(${PROJECT_SOURCE_DIR})

# Lista plików źródłowych CUDA
file(GLOB CUDA_SOURCES ${PROJECT_SOURCE_DIR}/*.cu)

# Dodanie głównego pliku wykonywalnego
add_executable(projekt ${CUDA_SOURCES})

# Znalezienie i linkowanie biblioteki CUDA
find_package(CUDAToolkit REQUIRED)
target_link_libraries(projekt CUDA::cudart)

# Ustawienia specyficzne dla CUDA
set_target_properties(projekt PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "all"
)