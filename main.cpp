#include <iostream>
#include <cmath>
#include <raylib.h>
#include "3DRenderer.hpp"

#include "devUtil.hpp"



int main(int argc, char** argv)
{
    float cubeSize = 850;
    float screenSize = 900;
    float nearZ = screenSize / 2;

    // DO NOT MODIFY IN CODE
    float cameraX = 0;
    float cameraY = 0;
    float cameraZ = 0;
    float cameraRoll = 0; // Degrees
    float FOV = 45;

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
            else
            {
                cameraX += drag.x * 5;
                cameraY -= drag.y * 5;
            }
        }
        if (IsKeyDown(KEY_W)) cameraZ -= 3;
        if (IsKeyDown(KEY_S)) cameraZ += 3;

        Vector2 proj[64];

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
            float halfCanvasSize = tanf(FOV * DEG2RAD / 2) * nearZ;

            Vec3f v = worldToCamera * verts[i];
            v.z -= halfCanvasSize + cubeSize / 2;

            proj[i].x = halfCanvasSize + (v.x * nearZ) / -v.z;
            proj[i].y = halfCanvasSize - (v.y * nearZ) / -v.z;
        }

        BeginDrawing();

        ClearBackground(BLACK);

        for (int i = 0; i < polyIdxs.size(); i++)
        {
            DrawLineV(proj[polyIdxs[i][0]], proj[polyIdxs[i][1]], RED);
            DrawLineV(proj[polyIdxs[i][1]], proj[polyIdxs[i][2]], RED);
            DrawLineV(proj[polyIdxs[i][2]], proj[polyIdxs[i][0]], RED);
        }

        EndDrawing();
    }

    return 0;
}