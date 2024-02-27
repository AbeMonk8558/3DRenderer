#include <cmath>
#include <vector>
#include <raylib.h>
#include "rasterizer.hpp"
#include "SIMD.hpp"
#include "math.hpp"

#include <fstream>
#include <iostream>
#include "devUtil.hpp"


Vec2f* Rasterizer::_proj = nullptr;
float* Rasterizer::_invZ = nullptr;
float** Rasterizer::_zBuffer = nullptr;
RenderTileList* Rasterizer::_tileBins = nullptr;
RenderTriangle* Rasterizer::_triangles = nullptr;


simd::float_m256 Rasterizer::getEdgeFunctionIncrements(const Vec2f& v1, const Vec2f& v2, float initial, const Vec2f& off)
{
    return simd::float_m256 { 4, 5 };
}

template <typename TVec>
TVec Rasterizer::pinedaEdge(const Vec2<TVec>& v1, const Vec2<TVec>& v2, const Vec2<TVec>& p)
{
    // Using determinants, we can calculate whether a point is the the right of left of an edge
    // (positive means right and vice versa). The 2D cross product is effectively the determinant
    // of a 2x2 matrix.
    return (v2 - v1).cross(p - v1);
}

void Rasterizer::perspectiveProject(const int& i, const simd::Vec3f_m256& p, const simd::Matrix44f_m256& worldToCamera, 
    const float& canvasSizeX, const float& canvasSizeY)
{
    simd::Vec3f_m256 v = worldToCamera * p;

    simd::float_m256 invZCurr = simd::float_m256(1) / -v.z;
    
    // Perform perspective divide, considering near clipping plane
    simd::Vec2f_m256 screen;
    screen.x = v.x * _nearZ * invZCurr;
    screen.y = v.y * _nearZ * invZCurr;

    // Normalized Device Coordinate in range [-1, 1]
    simd::Vec2f_m256 NDC;
    NDC.x = simd::float_m256(2) * screen.x / canvasSizeX;
    NDC.y = simd::float_m256(2) * screen.y / canvasSizeY;

    simd::Vec2f_m256 raster;
    raster.x = (simd::float_m256(1) + NDC.x) / 2 * (float)_screenWidth;
    raster.y = (simd::float_m256(1) + NDC.y) / 2 * (float)_screenHeight;

    invZCurr.storeData(&_invZ[i]);

    for (int j = 0; j < 8; j++)
    {
        _proj[i + j] = Vec2f(raster.x[j], raster.y[j]);
    }
}

Vec2f Rasterizer::getTrivialRejectOffset(const Vec2f& v1, const Vec2f& v2, int tileWidth, int tileHeight)
{
    // Uses the direction of the surface normal of the triangle's edge to determine
    // which corner of the tile should be used in testing.
    Vec2f n = (v2 - v1).surfaceNorm();

    return Vec2f(n.x < 0 ? 0 : tileWidth, n.y < 0 ? 0 : tileHeight);
}

Vec2f Rasterizer::getTrivialAcceptOffset(const Vec2f& r, int tileWidth, int tileHeight)
{
    // r is trivial reject value for the tile
    return Vec2f(((int)r.x + tileWidth) % (2 * tileWidth), 
        ((int)r.y + tileHeight) % (2 * tileHeight));
}

void Rasterizer::start()
{
    // char FPS_t[16];
    // SetTargetFPS(60);
    // InitWindow(_screenWidth, _screenHeight, "Test");

    // while (!WindowShouldClose())
    // {
    //     BeginDrawing();
    //     ClearBackground(BLACK);

    //     for (int y = 0; y < _screenHeight; y++)
    //     {
    //         for (int x = 0; x < _screenWidth; x++)
    //         {
    //             DrawPixel(x, _screenHeight - y - 1, BLUE);
    //         }
    //     }

    //     snprintf(FPS_t, 16, "FPS: %d", GetFPS());
    //     DrawText(FPS_t, _screenWidth - 100, 10, 15, BLACK);

    //     EndDrawing();
    // }

    // return;

    std::ofstream logger("C:\\Users\\alexa\\Documents\\3DRenderer\\Tests\\3DRendererTemp.txt");

    float FOV = 90;
    char FPS[16];

    SetTraceLogLevel(LOG_FATAL);
    SetTargetFPS(60);
    InitWindow(_screenWidth, _screenHeight, "Based 3D Render");

    std::vector<Vec3f> verts;
    std::vector<std::array<int, 3>> polyIdxs;
    retrievePrototypeMeshData(verts, polyIdxs, _cubeSize);
    for (Vec3f& v : verts)
        v.z -= _nearZ + _cubeSize / 2; // Makes them more easily visible

    // IMPORTANT: Overflows stack if not heap allocated
    _zBuffer = new float*[_screenHeight];
    for (int i = 0; i < _screenHeight; i++)
    {
        _zBuffer[i] = new float[_screenWidth];
    }
    
    const int nVerts = verts.size();
    const int nRemVerts = nVerts % 8;
    const int nPoly = polyIdxs.size();
    const int nTiles = (_screenWidth / TILE_WIDTH) * (_screenHeight / TILE_HEIGHT);

    // Since the number of vertices may not be divisible by 8 to evenly hold all 8-float AVX vectors, we
    // allocate a small buffer space at the end.
    _proj = new Vec2f[nVerts + 8 - nRemVerts];
    _invZ = new float[nVerts + 8 - nRemVerts];
    _tileBins = new RenderTileList[nTiles];

    // Set triangles to point to projected vertices
    _triangles = new RenderTriangle[nPoly];
    for (int i = 0; i < nPoly; i++)
    {
        _triangles[i].v1 = _proj[polyIdxs[i][0]];
        _triangles[i].v2 = _proj[polyIdxs[i][1]];
        _triangles[i].v3 = _proj[polyIdxs[i][2]];
    }

    simd::Matrix44f_m256 cameraToWorld = simd::Matrix44f_m256::identity();

    while (!WindowShouldClose())
    {
        float dCameraX = 0;
        float dCameraY = 0;
        float dCameraZ = 0;
        float dCameraPitch = 0; // Around x-axis
        float dCameraYaw = 0; // Around y-axis
        float dCameraRoll = 0; // Around z-axis

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
        {
            Vec2f drag = static_cast<Vec2f>(GetGestureDragVector());

            if (IsKeyDown(KEY_LEFT_CONTROL))
            {
                dCameraRoll -= drag.y * 3;
                dCameraYaw -= drag.x * 3;
            }
            else if (IsKeyDown(KEY_LEFT_SHIFT))
            {
                dCameraZ += drag.y * 5;
            }
            else
            {
                dCameraX += drag.x * 50;
                dCameraY -= drag.y * 50;
            }
        }
        if (IsKeyDown(KEY_W)) FOV += 1;
        if (IsKeyDown(KEY_S)) FOV -= 1;

        simd::Matrix44f_m256 xAxisRotation
        {
            1, 0, 0, 0,
            0, std::cos(dCameraPitch * DEG2RAD), -std::sin(dCameraPitch * DEG2RAD), 0,
            0, std::sin(dCameraPitch * DEG2RAD), std::cos(dCameraPitch * DEG2RAD), 0,
            0, 0, 0, 1
        };

        simd::Matrix44f_m256 yAxisRotation
        {
            std::cos(dCameraYaw * DEG2RAD), 0, std::sin(dCameraYaw * DEG2RAD), 0,
            0, 1, 0, 0,
            -std::sin(dCameraYaw * DEG2RAD), 0, std::cos(dCameraYaw * DEG2RAD), 0,
            0, 0, 0, 1
        };

        simd::Matrix44f_m256 zAxisRotation
        {
            std::cos(dCameraRoll * DEG2RAD), -std::sin(dCameraRoll * DEG2RAD), 0, 0,
            std::sin(dCameraRoll * DEG2RAD), std::cos(dCameraRoll * DEG2RAD), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        };

        cameraToWorld[0][3] += dCameraX;
        cameraToWorld[1][3] += dCameraY;
        cameraToWorld[2][3] += dCameraZ;

        cameraToWorld = zAxisRotation * yAxisRotation * xAxisRotation * cameraToWorld;
        simd::Matrix44f_m256 worldToCamera = cameraToWorld.inverse();

        const float canvasSizeX = 2 * std::tan(FOV * DEG2RAD / 2) * _nearZ;
        const float canvasSizeY = canvasSizeX * (_screenHeight / _screenWidth);

        for (int i = 0; i < verts.size(); i += 8) 
        {
            perspectiveProject(i, simd::vec3fToVec3f_m256(&verts[i]), worldToCamera, canvasSizeX,
                canvasSizeY);
        }

        if (nRemVerts > 0)
        {
            perspectiveProject(nVerts - nRemVerts, simd::vec3fToVec3f_m256(&verts[nVerts - nRemVerts], 
                nRemVerts), worldToCamera, canvasSizeX, canvasSizeY);
        }

        // All z-buffer coordinates are initially set to infinity so that the first z-value encountered
        // is garuanteed to be less than the current value.
        for (int y = 0; y < _screenHeight; y++)
        {
            for (int x = 0; x < _screenWidth; x++)
            {
                _zBuffer[y][x] = INFINITY;
            }
        }

        // Sort triangles into tile-based bins
        for (int i = 0; i < nPoly; i++)
        {
            RenderTriangle& poly = _triangles[i];

            // Projected triangles may be defined clockwise, reversing sign of edge function calculations
            poly.wSign = signf(pinedaEdge(poly.v1, poly.v2, poly.v3));

            // Trivial reject and accept offsets from bttm-left corner of tile
            poly.rOff1 = getTrivialRejectOffset(poly.v1, poly.v2, TILE_WIDTH, TILE_HEIGHT);
            poly.aOff1 = getTrivialAcceptOffset(poly.rOff1, TILE_WIDTH, TILE_HEIGHT);

            poly.rOff2 = getTrivialRejectOffset(poly.v2, poly.v3, TILE_WIDTH, TILE_HEIGHT);
            poly.aOff2 = getTrivialAcceptOffset(poly.rOff2, TILE_WIDTH, TILE_HEIGHT);

            poly.rOff3 = getTrivialRejectOffset(poly.v3, poly.v1, TILE_WIDTH, TILE_HEIGHT);
            poly.aOff3 = getTrivialAcceptOffset(poly.rOff3, TILE_WIDTH, TILE_HEIGHT);



            int tileCtr = 0;
            for (int y = 0; y < _screenHeight; y += TILE_HEIGHT)
            {
                for (int x = 0; x < _screenWidth; x += TILE_WIDTH)
                {
                    Vec2f p(x, y);
                    RenderTileList& tileList = _tileBins[tileCtr++];

                    float r1 = pinedaEdge(poly.v1, poly.v2, p + poly.rOff1) * poly.wSign;
                    float r2 = pinedaEdge(poly.v2, poly.v3, p + poly.rOff2) * poly.wSign;
                    float r3 = pinedaEdge(poly.v3, poly.v1, p + poly.rOff3) * poly.wSign;

                    if (r1 < 0 || r2 < 0 || r3 < 0)
                        continue;

                    _tileBins[tileCtr].extend();
                    RenderTileNode& tile = *(_tileBins[tileCtr++].last);

                    tile.polyIdx = i;
                    tile.x = x;
                    tile.y = y;

                    tile.a1 = pinedaEdge(poly.v1, poly.v2, p + poly.aOff1) * poly.wSign;
                    tile.a2 = pinedaEdge(poly.v2, poly.v3, p + poly.aOff2) * poly.wSign;
                    tile.a3 = pinedaEdge(poly.v3, poly.v1, p + poly.aOff3) * poly.wSign;

                    tile.r1 = r1;
                    tile.r2 = r2;
                    tile.r3 = r3;
                }
            }
        }

        BeginDrawing();
        ClearBackground(BLACK);

        // Operate on partitions of the initial tiles (descent)
        for (int i = 0; i < nTiles; i++)
        {
            int binSize = _tileBins[i].size;
            RenderTileNode* tile = _tileBins[i].root;

            while (tile != nullptr)
            {
                if (tile->isTriviallyAccepted())
                {
                    // Entire tile is contained within triangle, rasterize the whole thing
                    for (int y = 0; y < TILE_HEIGHT; y++)
                    {
                        for (int x = 0; x < TILE_WIDTH; x++)
                        {
                            DrawPixel(x, _screenHeight - y - 1, BLUE);
                        }
                    }
                }
                else
                {
                    const simd::float_m256& rIncr = increments[0];
                }

                tile = tile->next;
            }
        }

        snprintf(FPS, 16, "FPS: %d", GetFPS());
        DrawText(FPS, _screenWidth - 100, 10, 15, BLUE);

        EndDrawing();

        for (int i = 0; i < nPoly; i++)
        {
            _tileBins[i].clear();
        }
    }

    for (int i = 0; i < _screenHeight; i++)
    {
        delete[] _zBuffer[i];
    }
    delete[] _zBuffer;

    delete[] _proj;
    delete[] _invZ;

    logger.close();
}