#include <iostream>
#include "3DRenderer.hpp"

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

void printVec3f(const Vec3f& v)
{
    std::cout << "(" << v.x << ", " << v.y << ", " << v.z << ")" << std::endl;
}