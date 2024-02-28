#pragma once

#include <vector>
#include <array>
#include "math.hpp"
#include "renderTile.hpp"
#include "SIMD.hpp"

class Rasterizer
{
public:
    static void start();

private:
    constexpr static float _cubeSize = 400;
    constexpr static int _screenHeight = 560;
    constexpr static int _screenWidth = _screenHeight * 2;
    constexpr static float _nearZ = (float)_screenWidth / 2;

    constexpr static int _tileWidth = 32;
    constexpr static int _tileHeight = 16;

    static simd::float_m256 increments[4];

    static std::vector<Vec2f> _proj;
    static std::vector<float> _invZ;
    static std::vector<std::array<float, _screenWidth>> _zBuffer;
    static std::vector<RenderTileList> _tileBins;
    static std::vector<RenderTriangle> _triangles;

    static void setupStaticMembers();

    static float yToRaster(float y) noexcept;

    static simd::float_m256 getEdgeFunctionIncrements(const Vec2f& v1, const Vec2f& v2, float initial, const Vec2f& off);

    template <typename TVec>
    static TVec pinedaEdge(const Vec2<TVec>& v1, const Vec2<TVec>& v2, const Vec2<TVec>& p);

    static void perspectiveProject(const int& i, const simd::Vec3f_m256& p, const simd::Matrix44f_m256& worldToCamera, 
        const float& canvasSizeX, const float& canvasSizeY);

    static Vec2f getTrivialRejectOffset(const Vec2f& v1, const Vec2f& v2, int tileWidth, int tileHeight);

    static Vec2f getTrivialAcceptOffset(const Vec2f& trivialReject, int tileWidth, int tileHeight);
};