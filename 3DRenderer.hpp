#pragma once

#include <initializer_list>
#include <immintrin.h>
#include <raylib.h>
#include "linearMath.hpp"
#include "SIMD.hpp"

using Vec2f = Vec2<float>;
using Vec3f = Vec3<float>;
using Matrix44f = Matrix44<float>;

using Vec2f_m256 = Vec2<float_m256>;
using Vec3f_m256 = Vec3<float_m256>;
using Matrix44f_m256 = Matrix44<float_m256>;