#pragma once

#include <initializer_list>
#include <raylib.h>

template <typename T>
class Vec2;
template <typename T>
class Vec3;
template <typename T>
class Matrix44;

using Vec2f = Vec2<float>;
using Vec3f = Vec3<float>;
using Matrix44f = Matrix44<float>;

namespace simd
{
    class float_m256;

    using Vec2f_m256 = Vec2<simd::float_m256>;
    using Vec3f_m256 = Vec3<simd::float_m256>;
    using Matrix44f_m256 = Matrix44<simd::float_m256>;

    float_m256 max_m256(const float_m256& left, const float_m256& right);

    float_m256 min_m256(const float_m256& left, const float_m256& right);

    float_m256 abs_m256(const float_m256& f);

    template <typename T>
    T absf(const T& val);

    template <>
    float_m256 absf<float_m256>(const float_m256& fv);

    Vec3f_m256 vec3fToVec3f_m256(const Vec3f* vecs);

    Vec3f_m256 vec3fToVec3f_m256(const Vec3f* vecs, int num);
}