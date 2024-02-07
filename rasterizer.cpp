#include <cmath>
#include <vector>
#include <raylib.h>
#include "rasterizer.hpp"
#include "SIMD.hpp"
#include "linearMath.hpp"

#include <fstream>
#include <chrono>
#include <iostream>
#include "devUtil.hpp"

//#define LOG
#define LIGHTBLUE ((Color){ 173, 216, 230, 255 })

simd::Vec2f_m256* Rasterizer::_proj = nullptr;
float* Rasterizer::_invZ = nullptr;
float* Rasterizer::_zBuffer = nullptr;

template <typename TVec>
Vec2<TVec> tileTestCorner(const Vec2<TVec>& v1, const Vec2<TVec>& v2, const int& tSize)
{
    // Uses the direction of the surface normal of the triangle's edge to determine
    // which corner of the tile should be used in testing.
    Vec2<TVec> n = (v2 - v1).surfaceNorm();

    return Vec2<TVec>(n.x < 0 ? 0 : tSize, n.y < 0 ? 0 : tSize);
}

template <typename TVec>
void Rasterizer::perspectiveProject(const int& i, const Vec3<TVec>& p, const Matrix44<TVec>& worldToCamera, const float& canvasSize)
{
    Vec3<TVec> v = worldToCamera * p;

    TVec invZCurr = TVec(1) / -v.z;
    
    // Perform perspective divide, considering near clipping plane
    Vec2<TVec> screen;
    screen.x = v.x * _nearZ * invZCurr;
    screen.y = v.y * _nearZ * invZCurr;

    invZCurr.storeData(&_invZ[i]);

    // Normalized Device Coordinate in range [-1, 1]
    Vec2<TVec> NDC;
    NDC.x = TVec(2) * screen.x / TVec(canvasSize);
    NDC.y = TVec(2) * screen.y / TVec(canvasSize);

    _proj[i].x = (TVec(1) + NDC.x) / TVec(2) * (float)_screenSize; 
    _proj[i].y = (TVec(1) + NDC.y) / TVec(2) * (float)_screenSize;
}

template <typename TVec>
inline TVec Rasterizer::pinedaEdge(const Vec2<TVec>& v1, const Vec2<TVec>& v2, const Vec2<TVec>& p)
{
    // Using determinants, we can calculate whether a point is the the right of left of an edge
    // (positive means right and vice versa). The 2D cross product is effectively the determinant
    // of a 2x2 matrix.
    return (v2 - v1).cross(p - v1);
}

template <typename TVec>
inline TVec Rasterizer::pinedaEdgeGetInitial(const Vec2<TVec>& v1, const Vec2<TVec>& v2, const int& minx, const int& miny)
{
    return miny * (v2.x - v1.x) + v1.y * (minx - v2.x) + v2.y * (v1.x - minx);
}

template <typename TVec>
inline void Rasterizer::pinedaEdgeIncrX(TVec& a, const Vec2<TVec>& v1, const Vec2<TVec>& v2)
{
    a += v1.y - v2.y;
}

template <typename TVec>
inline void Rasterizer::pinedaEdgeSetY(TVec& a, const TVec& aInitial, const Vec2<TVec>& v1, const Vec2<TVec>& v2, const int& y, const int& miny)
{
    a = aInitial + (y - miny + 1) * (v2.x - v1.x);
}

bool Rasterizer::pointOnLine(const Vec2f& v1, const Vec2f& v2, const Vec2f& p)
{
    static const float epsilon = 1;

    Vec2f nv = (v2 - v1).normalize();
    float distSq = (p - (v1 + nv * nv.dot(p - v1))).lengthSq();
    return distSq < epsilon;
}

Vec2f Rasterizer::getSerialVec2(int idx)
{
    const simd::Vec2f_m256& vec = _proj[idx / 8];
    int vecIdx = idx % 8;

    return Vec2f(vec.x[vecIdx], vec.y[vecIdx]);
}

void Rasterizer::start()
{
    std::ofstream logger("C:\\Users\\alexa\\Documents\\3DRenderer\\Tests\\3DRendererTemp.txt");

    std::cout << "Entry" << std::endl;

    // DO NOT MODIFY IN CODE
    float cameraX = 0;
    float cameraY = 0;
    float cameraZ = 0;
    float cameraRoll = 0; // Degrees
    float FOV = 90;
    bool rasterize = false;

    char FPS[16];

    SetTraceLogLevel(LOG_FATAL);
    SetTargetFPS(60);
    InitWindow(_screenSize, _screenSize, "Based 3D Render");

    std::vector<Vec3f> verts;
    std::vector<std::array<int, 3>> polyIdxs;
    retrievePrototypeMeshData(verts, polyIdxs, _cubeSize);
    for (Vec3f& v : verts)
        v.z -= _nearZ + _cubeSize / 2; // Makes them more easily visible

    // IMPORTANT: Overflows stack if not heap allocated
    _zBuffer = new float[_screenSize * _screenSize];
    
    const int nVerts = verts.size();
    const int nVecVerts = (int)ceilf(nVerts / 8.0f);
    const int nRemVerts = nVerts % 8;

    _proj = new simd::Vec2f_m256[nVecVerts];
    _invZ = new float[nVerts + 8 - nRemVerts];

    while (!WindowShouldClose())
    {
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
        {
            Vec2f drag = static_cast<Vec2f>(GetGestureDragVector());

            if (IsKeyDown(KEY_LEFT_CONTROL))
            {
                cameraRoll -= drag.x * 3;
            }
            else if (IsKeyDown(KEY_LEFT_SHIFT))
            {
                cameraZ += drag.y * 5;
            }
            else
            {
                cameraX += drag.x * 50;
                cameraY -= drag.y * 50;
            }
        }
        if (IsKeyDown(KEY_W)) FOV += 1;
        if (IsKeyDown(KEY_S)) FOV -= 1;
        if (IsKeyPressed(KEY_R)) rasterize = !rasterize;

        simd::Matrix44f_m256 zAxisRotation
        {
            cosf(cameraRoll * DEG2RAD), -sinf(cameraRoll * DEG2RAD), 0, 0,
            sinf(cameraRoll * DEG2RAD), cosf(cameraRoll * DEG2RAD), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        };
        simd::Matrix44f_m256 cameraToWorld 
        {
            1, 0, 0, cameraX,
            0, 1, 0, cameraY,
            0, 0, 1, cameraZ,
            0, 0, 0, 1
        };
        simd::Matrix44f_m256 worldToCamera = (zAxisRotation * cameraToWorld).inverse();

        const float canvasSize = 2 * tanf(FOV * DEG2RAD / 2) * _nearZ;

        for (int i = 0; i < verts.size(); i += 8) 
        {
            perspectiveProject(i, simd::vec3fToVec3f_m256(&verts[i]), worldToCamera, canvasSize);
        }

        if (nRemVerts > 0)
        {
            perspectiveProject(nVerts - nRemVerts, simd::vec3fToVec3f_m256(&verts[nVerts - nRemVerts], 
                nRemVerts), worldToCamera, canvasSize);
        }

        // All z-buffer coordinates are initially set to infinity so that the first z-value encountered
        // is garuanteed to be less than the current value.
        for (int y = 0; y < _screenSize; y++)
        {
            for (int x = 0; x < _screenSize; x++)
            {
                _zBuffer[(y * _screenSize) + x] = INFINITY;
            }
        }
        
        BeginDrawing();
        ClearBackground(BLACK);

        for (int i = 0; i < polyIdxs.size(); i++)
        {
            if (rasterize)
            {
                const Vec2f& v1 = getSerialVec2(polyIdxs[i][0]);
                const Vec2f& v2 = getSerialVec2(polyIdxs[i][1]);
                const Vec2f& v3 = getSerialVec2(polyIdxs[i][2]);

                float invZ1 = _invZ[polyIdxs[i][0]];
                float invZ2 = _invZ[polyIdxs[i][1]];
                float invZ3 = _invZ[polyIdxs[i][2]];

                float miny = std::max(std::min(std::min(v1.y, v2.y), v3.y), 0.0f);
                float maxy = std::min(std::max(std::max(v1.y, v2.y), v3.y), _screenSize - 1.0f);

                float a1 = (v1.y - v2.y); // Since slope is -A/B in standard form, invert the y difference
                float b1 = (v2.x - v1.x);
                float c1 = b1 * v1.y + a1 * v1.x;

                float a2 = (v2.y - v3.y);
                float b2 = (v3.x - v2.x);
                float c2 = b2 * v2.y + a2 * v2.x;

                float a3 = (v3.y - v1.y);
                float b3 = (v1.x - v3.x);
                float c3 = b3 * v3.y + a3 * v3.x;

                for (int y = miny; y <= maxy; y++) 
                {
                    if (a1 == 0 && y == v1.y)
                    {
                        for (int x = std::min(v1.x, v2.x); x <= std::max(v1.x, v2.x); x++)
                            DrawPixel(x, _screenSize - y - 1, GRAY);

                        continue;
                    }
                    else if (a2 == 0 && y == v2.y)
                    {
                        for (int x = std::min(v2.x, v3.x); x <= std::max(v2.x, v3.x); x++)
                            DrawPixel(x, _screenSize - y - 1, GRAY);

                        continue;
                    }
                    else if (a3 == 0 && y == v3.y)
                    {
                        for (int x = std::min(v3.x, v1.x); x <= std::max(v3.x, v1.x); x++)
                            DrawPixel(x, _screenSize - y - 1, GRAY);

                        continue;
                    }
                
                    int x1 = std::round((c1 - b1 * y) / a1);
                    int x2 = std::round((c2 - b2 * y) / a2);
                    int x3 = std::round((c3 - b3 * y) / a3);
                    int minx = std::max(std::min(std::min(x1, x2), x3), 0);
                    int maxx = std::min(std::max(std::max(x1, x2), x3), _screenSize - 1);

                    logger << "Minx: " << minx << "\nMaxx: " << maxx << "\n";

                    for (int x = minx; x <= maxx; x++)
                    {   
                        // Barycentric coordinates allow z-values of pixels inside a triangle to be interpolated.
                        // Perspective-correct interpolation requires that we interpolate the reciprocal z-coordinates (since
                        // perspective projection preserves lines, but not distances).
                        float z = 0; //1 / (a1 * invA * invZ1 + a2 * invA * invZ2 + a3 * invA * invZ3);

                        if (z < _zBuffer[(y * _screenSize) + x])
                        {
                            _zBuffer[(y * _screenSize) + x] = z;

                            // if (pointOnLine(v1, v2, p) || pointOnLine(v2, v3, p) || pointOnLine(v3, v1, p))
                            //     DrawPixel(x, _screenSize - y - 1, LIGHTBLUE);
                            // else
                                DrawPixel(x, _screenSize - y - 1, GRAY);
                        }
                    }
                }
            }
        }

        snprintf(FPS, 16, "FPS: %d", GetFPS());
        DrawText(FPS, _screenSize - 100, 10, 15, BLUE);

        EndDrawing();
    }

    delete[] _zBuffer;
    delete[] _proj;
    delete[] _invZ;

    logger.close();
}