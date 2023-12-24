#pragma once

#include <vector>
#include "3DRenderer.hpp"

void printMatrix44f(const Matrix44f& m);
void printVec3f(const Vec3f& v);

void retrievePrototypeMeshData(std::vector<Vec3f>& verts, std::vector<std::vector<int>>& polyIdxs, float scale = 1.0f);