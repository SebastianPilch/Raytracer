cmake_minimum_required(VERSION 3.18)
project(ProjektCUDA LANGUAGES CXX CUDA)

# Set the C++ and CUDA standard versions
set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

# Add include directories
include_directories(${PROJECT_SOURCE_DIR})

# List of CUDA source files
file(GLOB CUDA_SOURCES ${PROJECT_SOURCE_DIR}/*.cu)

# Add the main executable
add_executable(projekt ${CUDA_SOURCES})

# Find CUDA Toolkit and specify static linking
find_package(CUDAToolkit REQUIRED)

# Determine the path to CUDA libraries
if(CUDAToolkit_FOUND)
    set(CUDA_LIB_DIR ${CUDAToolkit_LIBRARY_DIR})
else()
    message(FATAL_ERROR "CUDA Toolkit not found")
endif()

# Specify the CUDA runtime library to link statically
target_link_libraries(projekt
    ${CUDA_LIB_DIR}/cudart_static.lib
    ${CUDA_LIB_DIR}/cuda.lib
    ${CUDA_LIB_DIR}/cudadevrt.lib
)

# Specific settings for CUDA
set_target_properties(projekt PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "all"
)
