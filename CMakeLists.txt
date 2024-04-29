cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Ustawianie standardu C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Dodawanie folderu include
include_directories(include)

# Dodawanie ścieżki do plików pomocniczych
include_directories(helpers)
# Znajdowanie wszystkich plików źródłowych w folderze src
file(GLOB SOURCES "src/*.cpp")

# Tworzenie wykonywalnego pliku z plikami źródłowymi
add_executable(MyExecutable ${SOURCES} "include/ImportObj.h" "src/ImportObj.cpp" "include/vec3.h"  "include/Ray.h" "src/Ray.cpp")