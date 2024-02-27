#pragma once

#include <vector>
#include <array>
#include "math.hpp"

void printMatrix44f(const Matrix44f& m);
void printVec3f(const Vec3f& v);

void retrievePrototypeMeshData(std::vector<Vec3f>& verts, std::vector<std::array<int, 3>>& polyIdxs, float scale = 1.0f);