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
    static void perspectiveProject(const int& i, const Vec3<TVec>& p, const Matrix44<TVec>& worldToCamera, const float& canvasSize);

    template <typename TVec>
    static inline TVec pinedaEdge(const Vec2<TVec>& v1, const Vec2<TVec>& v2, const Vec2<TVec>& p);

    template <typename TVec>
    static inline TVec pinedaEdgeGetInitial(const Vec2<TVec>& v1, const Vec2<TVec>& v2, const int& minx, const int& miny);

    template <typename TVec>
    static inline void pinedaEdgeIncrX(TVec& a, const Vec2<TVec>& v1, const Vec2<TVec>& v2);

    template <typename TVec>
    static inline void pinedaEdgeSetY(TVec& a, const TVec& aInitial, const Vec2<TVec>& v1, const Vec2<TVec>& v2, const int& y, const int& miny);

    static bool pointOnLine(const Vec2f& v1, const Vec2f& v2, const Vec2f& p);

    static Vec2f getSerialVec2(int idx);
};