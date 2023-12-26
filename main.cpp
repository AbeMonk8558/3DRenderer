#include <iostream>
#include <cmath>
#include <raylib.h>
#include "3DRenderer.hpp"

#include "devUtil.hpp"

static float cubeSize = 850;
static float screenSize = 900;
static float nearZ = screenSize / 2;

float pinedaEdge(const Vec2f& v1, const Vec2f& v2, const Vec2f& p)
{
    return (v2 - v1).cross(p - v1);
}

Vector2 toRaylibRaster(const Vec2f& p)
{
    return (Vector2){p.x, screenSize - p.y - 1};
}

int main(int argc, char** argv)
{
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

        std::vector<Vec2f> proj;
        proj.reserve(polyIdxs.size() * 3);

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
            v.z -= nearZ + cubeSize / 2; // Makes them more easily visible

            Vec2f screen;
            screen.x = v.x * nearZ / -v.z;
            screen.y = v.y * nearZ / -v.z;

            Vec2f NDC;
            // [-1, 1]
            NDC.x = 2 * screen.x / canvasSize;
            NDC.y = 2 * screen.y / canvasSize;

            proj[i].x = (1 + NDC.x) / 2 * screenSize;
            proj[i].y = (1 + NDC.y) / 2 * screenSize;
        }

        BeginDrawing();

        ClearBackground(BLACK);

        char FPS[16];
        sprintf(FPS, "FPS: %d", GetFPS());
        DrawText(FPS, screenSize - 100, 10, 15, BLUE);

        for (int i = 0; i < polyIdxs.size(); i++)
        {
            if (rasterize)
            {
                Triangle2D t(proj[polyIdxs[i][0]], proj[polyIdxs[i][1]], proj[polyIdxs[i][2]]);
                AABB2D bbox = t.boundingBox();

                for (int y = bbox.min.y; y <= bbox.max.y; y++)
                {
                    for (int x = bbox.min.x; x <= bbox.max.x; x++)
                    {   
                        Vec2f p = {x, y};

                        float a = t.areaDoubled();
                        float a1 = pinedaEdge(t.v[0], t.v[1], p);
                        float a2 = pinedaEdge(t.v[1], t.v[2], p);
                        float a3 = pinedaEdge(t.v[2], t.v[0], p);

                        if (a1 < 0 || a2 < 0 || a3 < 0) continue;

                        // Barycentric coordinates
                        float b1 = a1 / a, b2 = a2 / a, b3 = a3 / a;

                        DrawPixelV(toRaylibRaster(p), GRAY);
                    }
                }   
            }

            DrawLineV(toRaylibRaster(proj[polyIdxs[i][0]]), toRaylibRaster(proj[polyIdxs[i][1]]), RED);
            DrawLineV(toRaylibRaster(proj[polyIdxs[i][1]]), toRaylibRaster(proj[polyIdxs[i][2]]), RED);
            DrawLineV(toRaylibRaster(proj[polyIdxs[i][2]]), toRaylibRaster(proj[polyIdxs[i][0]]), RED);
        }

        EndDrawing();
    }

    return 0;
}