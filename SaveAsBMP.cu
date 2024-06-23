#include "SaveAsBMP.cuh"
#include <iostream>
#include <fstream>
#include <string>

void saveAsBMP(float* angles, int width, int height, const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);

    if (!file) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    int paddingSize = (4 - (width * 3) % 4) % 4; // Padding required by BMP format

    // BMP header
    int filesize = 54 + (3 * width + paddingSize) * height;
    char fileHeader[54] = {
        'B', 'M',                         // Signature
        static_cast<char>(filesize), static_cast<char>(filesize >> 8), static_cast<char>(filesize >> 16), static_cast<char>(filesize >> 24), // File size
        0,0, 0,0,                         // Reserved
        54,0,0,0,                         // File offset to pixel array
        40,0,0,0,                         // DIB header size
        static_cast<char>(width), static_cast<char>(width >> 8), static_cast<char>(width >> 16), static_cast<char>(width >> 24), // Width
        static_cast<char>(height), static_cast<char>(height >> 8), static_cast<char>(height >> 16), static_cast<char>(height >> 24), // Height
        1,0,                              // Planes
        24,0,                             // Bits per pixel
        0,0,0,0,                          // Compression
        0,0,0,0,                          // Image size (can be 0 for uncompressed)
        0,0,0,0,                          // X pixels per meter (unused)
        0,0,0,0,                          // Y pixels per meter (unused)
        0,0,0,0,                          // Total colors (0 means default)
        0,0,0,0                           // Important colors (0 means all are important)
    };

    // Write header
    file.write(fileHeader, 54);

    // Write pixel data
    for (int i = height - 1; i >= 0; i--) {
        for (int j = 0; j < width; j++) {
            // Retrieve RGB values from angles array
            float red = angles[(i * width + j) * 3 + 0];
            float green = angles[(i * width + j) * 3 + 1];
            float blue = angles[(i * width + j) * 3 + 2];

            unsigned char r = static_cast<unsigned char>(red * 255.99f);
            unsigned char g = static_cast<unsigned char>(green * 255.99f);
            unsigned char b = static_cast<unsigned char>(blue * 255.99f);

            file.put(b); // Blue channel
            file.put(g); // Green channel
            file.put(r); // Red channel
        }

        // Add padding
        for (int k = 0; k < paddingSize; k++) {
            file.put(0);
        }
    }

    file.close();
}