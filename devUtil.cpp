#include <iostream>
#include <vector>
#include <array>
#include "math.hpp"

void printMatrix44f(const Matrix44f& m)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << m[i][j] << ' ';
        }

        std::cout << std::endl;
    }
}

// ******************* Code for storing prototype asteroid mesh in memory ******************
constexpr int numVerts = 8;
constexpr int numPoly = 12;

void retrievePrototypeMeshData(std::vector<Vec3f>& verts, std::vector<std::array<int, 3>>& polyIdxs, float scale)
{   
    verts.reserve(numVerts);
    polyIdxs.reserve(numPoly);
    
    // Vertices for a cube of size 1 centered around the origin
    verts =
    {
        {0.5, 0.5, 0.5}, // 0
        {-0.5, 0.5, 0.5}, // 1
        {-0.5, -0.5, 0.5}, // 2
        {0.5, -0.5, 0.5}, // 3
        {0.5, -0.5, -0.5}, // 4
        {-0.5, -0.5, -0.5}, // 5
        {-0.5, 0.5, -0.5}, // 6
        {0.5, 0.5, -0.5} // 7
    };

    for (Vec3f& v : verts)
        v = v * scale;

    polyIdxs =
    {
        // Front
        {0, 1, 3},
        {2, 3, 1},
        // Left
        {2, 5, 1},
        {6, 1, 5},
        // Back
        {7, 6, 4},
        {5, 4, 6},
        // Right
        {3, 4, 0},
        {7, 0, 4},
        // Top
        {1, 0, 6},
        {7, 6, 0},
        // Bottom
        {2, 3, 5},
        {4, 5, 3}
    };
}
// ****************************************************************************************