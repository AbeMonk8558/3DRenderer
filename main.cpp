#include <fstream>
#include <cmath>
#include <array>
#include <algorithm>
#include <raylib.h>
#include "SIMD.hpp"
#include "3DRenderer.hpp"

#include "devUtil.hpp"

/***************************** TODO LIST ******************************
1. Implement anti-aliasing (described in scratchapixel)
2. Implement SIMD optimization for matrix-vector multiplication, and for rasterizer
3. Make vertex-line stage more accurate.
4. Encase the entire main file in a class so variables can be shared amongst functions.
**********************************************************************/

constexpr static float cubeSize = 850;
constexpr static int screenSize = 900;
constexpr static float nearZ = (float)screenSize / 2;
std::ofstream logger("C:\\Users\\alexa\\Downloads\\3DRenderer.txt");

// Using determinants, we can calculate whether a point is the the right of left of an edge
// (positive means right and vice versa). The 2D cross product is effectively the determinant
// of a 2x2 matrix.
float pinedaEdge(const Vec2f& v1, const Vec2f& v2, const Vec2f& p)
{
    return (v2 - v1).cross(p - v1);
}

bool pointOnLine(const Vec2f& v1, const Vec2f& v2, const Vec2f& p)
{
    static const float epsilon = 1;

    const Vec2f nv = (v2 - v1).normalize();
    float distSq = (p - (v1 + nv * nv.dot(p - v1))).lengthSq();
    return distSq < epsilon;

    // Derived from setting the equation of a line in slope-intercept form when given a point and a slope to zero.
    //return std::abs((p.y - v1.y) * (v2.x - v1.x) - (p.x - v1.x) * (v2.y - v1.y)) < 2;
}

template <typename TSerial, typename TPacked>
TSerial getSerialValue(TPacked* arr, int idx)
{
    return arr[idx / 8][idx % 8];
}

int main(int argc, char** argv)
{
    printf("Entry\n");

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
    InitWindow(screenSize, screenSize, "Based 3D Render");

    std::vector<Vec3f> verts;
    std::vector<std::array<int, 3>> polyIdxs;
    retrievePrototypeMeshData(verts, polyIdxs, cubeSize);
    for (Vec3f& v : verts)
        v.z -= nearZ + cubeSize / 2; // Makes them more easily visible

    // IMPORTANT: Overflows stack if not heap allocated
    float* zBuffer = new float[screenSize * screenSize];
    //Color* frameBuffer = new Color[screenSize * screenSize];
    
    const int nVerts = verts.size();
    const int nVecVerts = (int)ceilf(nVerts / 8.0f);
    const int nRemVerts = nVerts % 8;

    Vec2f_m256* proj = new Vec2f_m256[nVecVerts];
    float* invZ = new float[nVerts + nRemVerts];

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

        Matrix44f_m256 cameraToWorld 
        {
            cosf(cameraRoll * DEG2RAD), -sinf(cameraRoll * DEG2RAD), 0, cameraX,
            sinf(cameraRoll * DEG2RAD), cosf(cameraRoll * DEG2RAD), 0, cameraY,
            0, 0, 1, cameraZ,
            0, 0, 0, 1
        };
        Matrix44f_m256 worldToCamera = cameraToWorld.inverse();

        float canvasSize = 2 * tanf(FOV * DEG2RAD / 2) * nearZ;

        for (int i = 0; i < verts.size(); i += 8) 
        {
            Vec3f_m256 v = worldToCamera * simd::vec3fToVec3f_m256(&verts[i]);

            simd::float_m256 invZCurr = simd::float_m256(1) / -v.z;
            
            // Perform perspective divide, considering near clipping plane
            Vec2f_m256 screen;
            screen.x = v.x * nearZ * invZCurr;
            screen.y = v.y * nearZ * invZCurr;

            // Normalized Device Coordinate in range [-1, 1]
            Vec2f_m256 NDC;
            NDC.x = simd::float_m256(2) * screen.x / canvasSize;
            NDC.y = simd::float_m256(2) * screen.y / canvasSize;

            invZCurr.storeData(&invZ[i]);

            proj[i].x = (simd::float_m256(1) + NDC.x) / simd::float_m256(2) * (float)screenSize;
            proj[i].y = (simd::float_m256(1) + NDC.y) / simd::float_m256(2) * (float)screenSize;
        }

        if (nRemVerts > 0)
        {
            int i = nVerts + nRemVerts - 8;

            Vec3f_m256 v = worldToCamera * simd::vec3fToVec3f_m256(&verts[nVerts - nRemVerts], nRemVerts);

            simd::float_m256 invZCurr = simd::float_m256(1) / -v.z;
            
            // Perform perspective divide, considering near clipping plane
            Vec2f_m256 screen;
            screen.x = v.x * nearZ * invZCurr;
            screen.y = v.y * nearZ * invZCurr;

            // Normalized Device Coordinate in range [-1, 1]
            Vec2f_m256 NDC;
            NDC.x = simd::float_m256(2) * screen.x / canvasSize;
            NDC.y = simd::float_m256(2) * screen.y / canvasSize;

            invZCurr.storeData(&invZ[nVerts + nRemVerts - 8]);

            proj[i].x = (simd::float_m256(1) + NDC.x) / simd::float_m256(2) * (float)screenSize;
            proj[i].y = (simd::float_m256(1) + NDC.y) / simd::float_m256(2) * (float)screenSize;
        }

        // All z-buffer coordinates are initially set to infinity so that the first z-value encountered
        // is garuanteed to be less than the current value.
        for (int y = 0; y < screenSize; y++)
        {
            for (int x = 0; x < screenSize; x++)
            {
                zBuffer[(y * screenSize) + x] = INFINITY;
                //frameBuffer[(y * screenSize) + x] = BLACK;
            }
        }
        
        BeginDrawing();
        ClearBackground(BLACK);

        for (int i = 0; i < polyIdxs.size(); i++)
        {
            if (rasterize)
            {
                const Vec2f& v1 = getSerialValue<Vec2f>(proj, polyIdxs[i][0]);
                const Vec2f& v2 = getSerialValue<Vec2f>(proj, polyIdxs[i][1]);
                const Vec2f& v3 = getSerialValue<Vec2f>(proj, polyIdxs[i][2]);

                float invZ1 = invZ[polyIdxs[i][0]];
                float invZ2 = invZ[polyIdxs[i][1]];
                float invZ3 = invZ[polyIdxs[i][2]];

                float minx = std::max(std::min(std::min(v1.x, v2.x), v3.x), 0.0f);
                float miny = std::max(std::min(std::min(v1.y, v2.y), v3.y), 0.0f);
                float maxx = std::min(std::max(std::max(v1.x, v2.x), v3.x), screenSize - 1.0f);
                float maxy = std::min(std::max(std::max(v1.y, v2.y), v3.y), screenSize - 1.0f);

                float a = pinedaEdge(v1, v2, v3);

                // Triangles, when projected onto the screen, or because of a 3D rotation, could be defined
                // in clockwise order relative to the camera. This means all calculations which should be
                // positive will be negative, and their signs must be reversed. We can tell by whether the
                // triangle's area is negative that this is the case.
                float windingSign = a < 0 ? -1 : 1;
                a *= windingSign;

                for (int y = miny; y <= maxy; y++)
                {
                    for (int x = minx; x <= maxx; x++)
                    {   
                        Vec2f p(x, y);

                        // Edge function determines whether point lies inside of the triangle using
                        // the determinant (which can determine rotational relationships). Note that
                        // the barycentric coordinate for each vertex is computed using the area between its
                        // opposite edge and the point, hence the arrangment of the vertices in the below lines.
                        float a1 = pinedaEdge(v2, v3, p) * windingSign;
                        float a2 = pinedaEdge(v3, v1, p) * windingSign;
                        float a3 = pinedaEdge(v1, v2, p) * windingSign;

                        if (a1 < 0 || a2 < 0 || a3 < 0) continue;

                        // Barycentric coordinates allow z-values of pixels inside a triangle to be interpolated.
                        a1 /= a;
                        a2 /= a;
                        a3 /= a;

                        // Perspective-correct interpolation requires that we interpolate the reciprocal z-coordinates (since
                        // perspective projection preserves lines, but not distances).
                        float z = 1 / (a1 * invZ1 + a2 * invZ2 + a3 * invZ3);

                        if (z < zBuffer[(y * screenSize) + x])
                        {
                            zBuffer[(y * screenSize) + x] = z;

                            if (pointOnLine(v1, v2, p) || pointOnLine(v2, v3, p) || pointOnLine(v3, v1, p))
                                DrawPixel(x, screenSize - y - 1, RED);
                                //frameBuffer[(y * screenSize) + x] = RED;
                            else
                                DrawPixel(x, screenSize - y - 1, GRAY);
                                //frameBuffer[(y * screenSize) + x] = GRAY;
                        }

                        // logger << (p.x - t.v[0].x) / (t.v[1].x - t.v[0].x) << ", " << (p.y - t.v[0].y) / (t.v[1].y - t.v[0].y) << '\n';
                    }
                }
            }
        }

        // for (int y = 0; y < screenSize; y++)
        // {
        //     for (int x = 0; x < screenSize; x++)
        //         DrawPixel(x, screenSize - y - 1, frameBuffer[(y * screenSize) + x]);
        // }

        snprintf(FPS, 16, "FPS: %d", GetFPS());
        DrawText(FPS, screenSize - 100, 10, 15, BLUE);

        EndDrawing();
    }

    delete[] zBuffer;
    //delete[] frameBuffer;
    delete[] proj;
    delete[] invZ;

    logger.close();
    return 0;
}