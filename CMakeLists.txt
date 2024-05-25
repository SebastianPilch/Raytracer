cmake_minimum_required(VERSION 3.18)
project(ProjektCUDA LANGUAGES CXX CUDA)

# Ustawienie wersji CUDA, jeżeli jest wymagane
set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

# Dodanie ścieżki do plików nagłówkowych
include_directories(${PROJECT_SOURCE_DIR}/include)

# Lista plików źródłowych CUDA
file(GLOB_RECURSE CUDA_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cu ${PROJECT_SOURCE_DIR}/helpers/*.cu)

# Dodanie głównego pliku wykonywalnego
add_executable(projekt ${CUDA_SOURCES})

# Ustawienia specyficzne dla CUDA
set_target_properties(projekt PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "all"
)