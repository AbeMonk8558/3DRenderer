#include <cmath>
#include <vector>
#include <raylib.h>
#include "rasterizer.hpp"
#include "SIMD.hpp"
#include "linearMath.hpp"

#include <fstream>
#include <iostream>
#include "devUtil.hpp"


simd::Vec2f_m256* Rasterizer::_proj = nullptr;
float* Rasterizer::_invZ = nullptr;
float* Rasterizer::_zBuffer = nullptr;

template <typename TVec>
TVec Rasterizer::pinedaEdge(const Vec2<TVec>& v1, const Vec2<TVec>& v2, const Vec2<TVec>& p)
{
    // Using determinants, we can calculate whether a point is the the right of left of an edge
    // (positive means right and vice versa). The 2D cross product is effectively the determinant
    // of a 2x2 matrix.
    return (v2 - v1).cross(p - v1);
}

__m256 Rasterizer::pointOnLine(const simd::Vec2f_m256& v1, const simd::Vec2f_m256& v2, const simd::Vec2f_m256& p)
{
    static const simd::float_m256 epsilon(1);

    simd::Vec2f_m256 nv = (v2 - v1).normalize();
    simd::float_m256 distSq = (p - (v1 + nv * nv.dot(p - v1))).lengthSq();
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
    std::ofstream logger("C:\\Users\\alexa\\Downloads\\3DRenderer.txt");

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
    _invZ = new float[nVerts + nRemVerts];

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

        float canvasSize = 2 * tanf(FOV * DEG2RAD / 2) * _nearZ;

        for (int i = 0; i < verts.size(); i += 8) 
        {
            simd::Vec3f_m256 v = worldToCamera * simd::vec3fToVec3f_m256(&verts[i]);

            simd::float_m256 invZCurr = simd::float_m256(1) / -v.z;
            
            // Perform perspective divide, considering near clipping plane
            simd::Vec2f_m256 screen;
            screen.x = v.x * _nearZ * invZCurr;
            screen.y = v.y * _nearZ * invZCurr;

            // Normalized Device Coordinate in range [-1, 1]
            simd::Vec2f_m256 NDC;
            NDC.x = simd::float_m256(2) * screen.x / canvasSize;
            NDC.y = simd::float_m256(2) * screen.y / canvasSize;

            invZCurr.storeData(&_invZ[i]);

            _proj[i].x = (simd::float_m256(1) + NDC.x) / simd::float_m256(2) * (float)_screenSize;
            _proj[i].y = (simd::float_m256(1) + NDC.y) / simd::float_m256(2) * (float)_screenSize;
        }

        if (nRemVerts > 0)
        {
            int i = nVerts + nRemVerts - 8;

            simd::Vec3f_m256 v = worldToCamera * simd::vec3fToVec3f_m256(&verts[nVerts - nRemVerts], nRemVerts);

            simd::float_m256 invZCurr = simd::float_m256(1) / -v.z;
            
            // Perform perspective divide, considering near clipping plane
            simd::Vec2f_m256 screen;
            screen.x = v.x * _nearZ * invZCurr;
            screen.y = v.y * _nearZ * invZCurr;

            // Normalized Device Coordinate in range [-1, 1]
            simd::Vec2f_m256 NDC;
            NDC.x = simd::float_m256(2) * screen.x / canvasSize;
            NDC.y = simd::float_m256(2) * screen.y / canvasSize;

            invZCurr.storeData(&_invZ[nVerts + nRemVerts - 8]);

            _proj[i].x = (simd::float_m256(1) + NDC.x) / 2 * (float)_screenSize;
            _proj[i].y = (simd::float_m256(1) + NDC.y) / 2 * (float)_screenSize;
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
                Vec2f v1Serial = getSerialVec2(polyIdxs[i][0]);
                Vec2f v2Serial = getSerialVec2(polyIdxs[i][1]);
                Vec2f v3Serial = getSerialVec2(polyIdxs[i][2]);

                float minx = std::max(std::min(std::min(v1Serial.x, v2Serial.x), v3Serial.x), 0.0f);
                float miny = std::max(std::min(std::min(v1Serial.y, v2Serial.y), v3Serial.y), 0.0f);
                float maxx = std::min(std::max(std::max(v1Serial.x, v2Serial.x), v3Serial.x), _screenSize - 1.0f);
                float maxy = std::min(std::max(std::max(v1Serial.y, v2Serial.y), v3Serial.y), _screenSize - 1.0f);

                float aSerial = pinedaEdge(v1Serial, v2Serial, v3Serial);

                // Triangles, when projected onto the screen, or because of a 3D rotation, could be defined
                // in clockwise order relative to the camera. This means all calculations which should be
                // positive will be negative, and their signs must be reversed. We can tell by whether the
                // triangle's area is negative that this is the case.
                float windingSignSerial = aSerial < 0 ? -1 : 1;
                aSerial *= windingSignSerial;

                simd::Vec2f_m256 v1(simd::float_m256(v1Serial.x), simd::float_m256(v1Serial.y));
                simd::Vec2f_m256 v2(simd::float_m256(v2Serial.x), simd::float_m256(v2Serial.y));
                simd::Vec2f_m256 v3(simd::float_m256(v3Serial.x), simd::float_m256(v3Serial.y));

                simd::float_m256 a(aSerial);
                simd::float_m256 windingSign(windingSignSerial);

                simd::float_m256 invZ1(_invZ[polyIdxs[i][0]]);
                simd::float_m256 invZ2(_invZ[polyIdxs[i][1]]);
                simd::float_m256 invZ3(_invZ[polyIdxs[i][2]]);

                int ctr = 0;
                simd::Vec2f_m256 p;

                for (int y = miny; y <= maxy; y++)
                {
                    for (int x = minx; x <= maxx; x++)
                    {   
                        if (ctr < 8)
                        {
                            p.x[ctr] = x;
                            p.y[ctr] = y;
                            //logger << '(' << x << ',' << y << ')' << '\n';

                            ctr++;

                            if (y < maxy && x < maxx && ctr != 8)
                                continue;
                        }

                        // Edge function determines whether point lies inside of the triangle using
                        // the determinant (which can determine rotational relationships). Note that
                        // the barycentric coordinate for each vertex is computed using the area between its
                        // opposite edge and the point, hence the arrangment of the vertices in the below lines.
                        simd::float_m256 a1 = pinedaEdge(v2, v3, p) * windingSign;
                        simd::float_m256 a2 = pinedaEdge(v3, v1, p) * windingSign;
                        simd::float_m256 a3 = pinedaEdge(v1, v2, p) * windingSign;

                        // Barycentric coordinates allow z-values of pixels inside a triangle to be interpolated.
                        a1 /= a;
                        a2 /= a;
                        a3 /= a;

                        // Perspective-correct interpolation requires that we interpolate the reciprocal z-coordinates (since
                        // perspective projection preserves lines, but not distances).
                        simd::float_m256 z = simd::float_m256(1) / (a1 * invZ1 + a2 * invZ2 + a3 * invZ3);

                        __m256 inTriangle = _mm256_and_ps(_mm256_and_ps(a1 >= 0, a2 >= 0), a3 >= 0);
                        __m256 onEdge = _mm256_or_ps(_mm256_or_ps(pointOnLine(v1, v2, p), pointOnLine(v2, v3, p)), pointOnLine(v3, v1, p));

                        for (int j = 0; j < ctr; j++)
                        {
                            int px = (int)p.x[j];
                            int py = (int)p.y[j];

                            if (inTriangle[j] && z[j] < _zBuffer[(py * _screenSize) + px])
                            {
                                _zBuffer[(py * _screenSize) + px] = z[j];

                                if (onEdge[j])
                                    DrawPixel(px, _screenSize - py - 1, RED);
                                else
                                    DrawPixel(px, _screenSize - py - 1, GRAY);
                            }
                        }

                        ctr = 0;
                    }
                }
            }
        }

        snprintf(FPS, 16, "FPS: %d", GetFPS());
        DrawText(FPS, _screenSize - 100, 10, 15, BLUE);

        EndDrawing();
    }

    delete[] _zBuffer;
    //delete[] frameBuffer;
    delete[] _proj;
    delete[] _invZ;

    logger.close();
}