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
#define BARYCENTRIC_EDGE_TEST_EPSILON 5.0f

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

        if (rasterize)
        {
            for (int i = 0; i < polyIdxs.size(); i++)
            {
                const Vec2f& v1 = getSerialVec2(polyIdxs[i][0]);
                const Vec2f& v2 = getSerialVec2(polyIdxs[i][1]);
                const Vec2f& v3 = getSerialVec2(polyIdxs[i][2]);

                float invZ1 = _invZ[polyIdxs[i][0]];
                float invZ2 = _invZ[polyIdxs[i][1]];
                float invZ3 = _invZ[polyIdxs[i][2]];

                float minx = std::max(std::min(std::min(v1.x, v2.x), v3.x), 0.0f);
                float miny = std::max(std::min(std::min(v1.y, v2.y), v3.y), 0.0f);
                float maxx = std::min(std::max(std::max(v1.x, v2.x), v3.x), _screenSize - 1.0f);
                float maxy = std::min(std::max(std::max(v1.y, v2.y), v3.y), _screenSize - 1.0f);

                float a = pinedaEdge(v1, v2, v3);
                float invA = 1 / a;

                // Triangles, when projected onto the screen, or because of a 3D rotation, could be defined
                // in clockwise order relative to the camera. This means all calculations which should be
                // positive will be negative, and their signs must be reversed. We can tell by whether the
                // triangle's area is negative that this is the case ([Winding]Sign).
                float wSign = invA < 0 ? -1 : 1;
                float invSqrtA = 1 / std::sqrt(a * wSign);

                const float a1Initial = pinedaEdgeGetInitial(v2, v3, minx, miny);
                const float a2Initial = pinedaEdgeGetInitial(v3, v1, minx, miny);
                const float a3Initial = pinedaEdgeGetInitial(v1, v2, minx, miny);

                // Barycentric coordinates represent weights towards each vertex (add up to 1)
                float a1 = a1Initial; 
                float a2 = a2Initial;
                float a3 = a3Initial;

                for (int y = miny; y <= maxy; y++)
                {
                    for (int x = minx; x <= maxx; x++)
                    {   
                        // Edge function determines whether point lies inside of the triangle using
                        // the determinant (which can determine rotational relationships). Note that
                        // the barycentric coordinate for each vertex is computed using the area between its
                        // opposite edge and the point, hence the arrangment of the vertices in the below lines.
                        if (a1 * wSign < 0 || a2 * wSign < 0 || a3 * wSign < 0) 
                        {
                            pinedaEdgeIncrX(a1, v2, v3);
                            pinedaEdgeIncrX(a2, v3, v1);
                            pinedaEdgeIncrX(a3, v1, v2);
                            continue;
                        }

                        // Normalized barycentric coordinates
                        float a1n = a1 * invA;
                        float a2n = a2 * invA;
                        float a3n = a3 * invA;

                        // Barycentric coordinates allow z-values of pixels inside a triangle to be interpolated.
                        // Perspective-correct interpolation requires that we interpolate the reciprocal z-coordinates (since
                        // perspective projection preserves lines, but not distances).
                        float z = 1 / (a1n * invZ1 + a2n * invZ2 + a3n * invZ3);

                        if (z < _zBuffer[(y * _screenSize) + x])
                        {
                            _zBuffer[(y * _screenSize) + x] = z;

                            // If a barycentric coordinate is near-zero, that means it lies close to (or on) an edge.
                            // We scale the epsilon inversely proportional to the sqrt of the triangle's area because a
                            // smaller triangle's barycentric coordinates will distort absolute distance to edges (since
                            // such coordinates represent proportional relationships). Scaling by some multiple
                            // of the triangle's area is required. However, area itself scales quadratically
                            // relative to distance, which is problematic. The sqrt gives a better approximation relative
                            // to side length.
                            if (std::min(std::min(a1n, a2n), a3n) < invSqrtA * BARYCENTRIC_EDGE_TEST_EPSILON)
                                DrawPixel(x, _screenSize - y - 1, LIGHTBLUE);
                            else
                                DrawPixel(x, _screenSize - y - 1, GRAY);
                        }

                        pinedaEdgeIncrX(a1, v2, v3);
                        pinedaEdgeIncrX(a2, v3, v1);
                        pinedaEdgeIncrX(a3, v1, v2);
                    }

                    pinedaEdgeSetY(a1, a1Initial, v2, v3, y, miny);
                    pinedaEdgeSetY(a2, a2Initial, v3, v1, y, miny);
                    pinedaEdgeSetY(a3, a3Initial, v1, v2, y, miny);
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