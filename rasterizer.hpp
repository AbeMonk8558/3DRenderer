#pragma once

#include "linearMath.hpp"
#include "SIMD.hpp"

class Rasterizer
{
public:
    static void start();

private:
    constexpr static float _cubeSize = 850;
    constexpr static int _screenSize = 900;
    constexpr static float _nearZ = (float)_screenSize / 2;

    static simd::Vec2f_m256* _proj;
    static float* _invZ;
    static float* _zBuffer;

    template <typename TVec>
    static TVec pinedaEdge(const Vec2<TVec>& v1, const Vec2<TVec>& v2, const Vec2<TVec>& p);

    static __m256 pointOnLine(const simd::Vec2f_m256& v1, const simd::Vec2f_m256& v2, const simd::Vec2f_m256& p);

    static Vec2f getSerialVec2(int idx);
};