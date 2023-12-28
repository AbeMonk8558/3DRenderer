#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <raylib.h>
#include "3DRenderer.hpp"

#include "devUtil.hpp"

constexpr static float cubeSize = 850;
constexpr static int screenSize = 900;
constexpr static float nearZ = (float)screenSize / 2;

// Using determinants, we can calculate whether a point is the the right of left of an edge
// (positive means right and vice versa). The 2D cross product is effectively the determinant
// of a 2x2 matrix.
float pinedaEdge(const Vec2f& v1, const Vec2f& v2, const Vec2f& p)
{
    return (v2 - v1).cross(p - v1);
}

// Invert y-coordinate during conversion from screen to raster space.
Vector2 toRaylibRaster(const Vec2f& p)
{
    return (Vector2){p.x, screenSize - p.y - 1};
}

// Derived from setting the equation of a line in slope-intercept form when given a point and a slope to zero.
bool pointOnLine(const Vec2f& v1, const Vec2f& v2, const Vec2f& p)
{
    // Uses epsilon to avoid floating-point imprecision errors.
    return std::abs((p.y - v1.y) * (v2.x - v1.x) - (p.x - v1.x) * (v2.y - v1.y)) < 2;
}

int main(int argc, char** argv)
{
    printf("Entry\n");

    // Vec2f v0(1009.435345, 5674567.2142134);
    // Vec2f v1(506.12341, 4567567.2345);
    // Vec2f p = v0 + (v1 - v0) * 0.3452345;

    // std::cout << '(' << p.x << ", " << p.y << ')' << std::endl;
    // std::cout << pointOnLine(v0, v1, p) << std::endl;
    // return 0;

    // DO NOT MODIFY IN CODE
    float cameraX = 0;
    float cameraY = 0;
    float cameraZ = 0;
    float cameraRoll = 0; // Degrees
    float FOV = 90;
    bool rasterize = false;

    SetTraceLogLevel(LOG_FATAL);
    SetTargetFPS(60);
    InitWindow(screenSize, screenSize, "Based 3D Render");

    std::vector<Vec3f> verts;
    std::vector<std::vector<int>> polyIdxs;
    retrievePrototypeMeshData(verts, polyIdxs, cubeSize);
    for (Vec3f& v : verts)
        v.z -= nearZ + cubeSize / 2; // Makes them more easily visible

    std::vector<std::array<float, screenSize>> zBuffer;
    zBuffer.reserve(screenSize);
    std::vector<Vec2f> proj;
    proj.reserve(polyIdxs.size() * 3);

    //std::ofstream logger("C:\\Users\\alexa\\Downloads\\3DRenderer.txt");

    while (!WindowShouldClose())
    {
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
        {
            Vector2 drag = GetGestureDragVector();

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
                cameraX += drag.x * 10;
                cameraY -= drag.y * 10;
            }
        }
        if (IsKeyDown(KEY_W)) FOV += 1;
        if (IsKeyDown(KEY_S)) FOV -= 1;
        if (IsKeyPressed(KEY_R)) rasterize = !rasterize;

        Matrix44f zAxisRotation
        {
            cosf(cameraRoll * DEG2RAD), -sinf(cameraRoll * DEG2RAD), 0, 0,
            sinf(cameraRoll * DEG2RAD), cosf(cameraRoll * DEG2RAD), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        };
        Matrix44f cameraToWorld 
        {
            1, 0, 0, cameraX,
            0, 1, 0, cameraY,
            0, 0, 1, cameraZ,
            0, 0, 0, 1
        };
        Matrix44f worldToCamera = (zAxisRotation * cameraToWorld).inverse();

        for (int i = 0; i < verts.size(); i++) 
        {
            float canvasSize = 2 * tanf(FOV * DEG2RAD / 2) * nearZ;

            Vec3f v = worldToCamera * verts[i];

            Vec2f screen;
            screen.x = v.x * nearZ / -v.z;
            screen.y = v.y * nearZ / -v.z;

            Vec2f NDC;
            // [-1, 1]
            NDC.x = 2 * screen.x / canvasSize;
            NDC.y = 2 * screen.y / canvasSize;

            proj[i].x = (1 + NDC.x) / 2 * (float)screenSize;
            proj[i].y = (1 + NDC.y) / 2 * (float)screenSize;
        }

        BeginDrawing();

        ClearBackground(BLACK);

        char FPS[16];
        snprintf(FPS, 16, "FPS: %d", GetFPS());
        DrawText(FPS, screenSize - 100, 10, 15, BLUE);

        // All z-buffer coordinates are initially set to infinity so that the first z-value encountered
        // is garuanteed to be less than the current value.
        for (int i = 0; i < screenSize; i++)
            for (int j = 0; j < screenSize; j++)
                zBuffer[i][j] = INFINITY;

        for (int i = 0; i < polyIdxs.size(); i++)
        {
            if (rasterize)
            {
                Triangle2D t(proj[polyIdxs[i][0]], proj[polyIdxs[i][1]], proj[polyIdxs[i][2]]);
                AABB2D bbox = t.boundingBox();
                float a = t.areaDoubled();

                // Triangles, when projected onto the screen, or because of a 3D rotation, could be defined
                // in clockwise order relative to the camera. This means all calculations which should be
                // positive will be negative, and their signs must be reversed. We can tell by whether the
                // triangle's area is negative that this is the case.
                float windingSign = a < 0 ? -1 : 1;
                a *= windingSign;

                for (int y = bbox.min.y; y <= bbox.max.y; y++)
                {
                    for (int x = bbox.min.x; x <= bbox.max.x; x++)
                    {   
                        Vec2f p(x, y);

                        float a1 = pinedaEdge(t.v[0], t.v[1], p) * windingSign;
                        float a2 = pinedaEdge(t.v[1], t.v[2], p) * windingSign;
                        float a3 = pinedaEdge(t.v[2], t.v[0], p) * windingSign;

                        if (a1 < 0 || a2 < 0 || a3 < 0) continue;

                        // Barycentric coordinates allow z-values of pixels inside a triangle to be interpolated.
                        float b1 = a1 / a, b2 = a2 / a, b3 = a3 / a;
                        float z = b1 * verts[polyIdxs[i][0]].z + b2 * verts[polyIdxs[i][1]].z + b3 * verts[polyIdxs[i][2]].z;

                        if (pointOnLine(t.v[0], t.v[1], p) || pointOnLine(t.v[1], t.v[2], p) || pointOnLine(t.v[2], t.v[0], p))
                            DrawPixelV(toRaylibRaster(p), RED);
                        else
                            DrawPixelV(toRaylibRaster(p), GRAY);

                        // logger << (p.x - t.v[0].x) / (t.v[1].x - t.v[0].x) << ", " << (p.y - t.v[0].y) / (t.v[1].y - t.v[0].y) << '\n';
                    }
                }   
            }

            // DrawLineV(toRaylibRaster(proj[polyIdxs[i][0]]), toRaylibRaster(proj[polyIdxs[i][1]]), RED);
            // DrawLineV(toRaylibRaster(proj[polyIdxs[i][1]]), toRaylibRaster(proj[polyIdxs[i][2]]), RED);
            // DrawLineV(toRaylibRaster(proj[polyIdxs[i][2]]), toRaylibRaster(proj[polyIdxs[i][0]]), RED);
        }

        EndDrawing();
    }

    // logger.close();
    return 0;
}