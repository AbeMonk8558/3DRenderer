#pragma once

#include <vector>
#include <array>
#include "math.hpp"

void printMatrix44f(const Matrix44f& m);

void retrievePrototypeMeshData(std::vector<Vec3f>& verts, std::vector<std::array<int, 3>>& polyIdxs, float scale = 1.0f);

static inline void printVec3f(const Vec3f& v)
{
    std::cout << "(" << v.x << ", " << v.y << ", " << v.z << ")" << std::endl;
}

static inline void printVec2f(const Vec2f& v)
{
    std::cout << "(" << v.x << ", " << v.y << ")" << std::endl;
}

template <typename T>
static inline void printPrim(T prim)
{
    std::cout << prim << std::endl;
}