#include <iostream>
#include <raylib.h>
#include "3DRenderer.hpp"

#include "devUtil.hpp"



int main(int argc, char** argv)
{
    float xScale = 50;
    float yScale = 50;
    float zScale = 50;
    float zTranslation = 250; // Units away from apex along z-axis

    float canvasSize = 500;
    float cameraX = 0;
    float cameraY = 0;
    float cameraZ = 0;

    SetTraceLogLevel(LOG_FATAL);
    SetTargetFPS(60);
    InitWindow(canvasSize, canvasSize, "Based 3D Render");

    // ******* Auto-generates cube vertices in an autistic manner *******
    Vec3f vertices[8]; // 12 edges, 2 vertices per edge

    int vIdx = 0;
    for (int x = -xScale / 2; x <= xScale / 2; x += xScale)
    {
        for (int y = -yScale / 2; y <= yScale / 2; y += yScale)
        {
            for (int z = -zTranslation; z >= -(zScale + zTranslation); z -= zScale)
            {
                vertices[vIdx++] = Vec3f(x, y, z);
            }
        }
    }
    // ******************************************************************

    while (!WindowShouldClose())
    {
        if (IsKeyDown(KEY_UP)) cameraY++;
        else if (IsKeyDown(KEY_DOWN)) cameraY--;
        else if (IsKeyDown(KEY_LEFT)) cameraX--;
        else if (IsKeyDown(KEY_RIGHT)) cameraX++;
        else if (IsKeyDown(KEY_W)) cameraZ--;
        else if (IsKeyDown(KEY_S)) cameraZ++;

        Vector2 proj[8];
        Matrix44f cameraToWorld 
        {
            1, 0, 0, cameraX,
            0, 1, 0, cameraY,
            0, 0, 1, cameraZ,
            0, 0, 0, 1
        };
        Matrix44f worldToCamera = cameraToWorld.inverse();

        for (int i = 0; i < 8; i++) 
        {
            const Vec3f& v = worldToCamera * vertices[i];

            proj[i].x = (canvasSize / 2) + (v.x * canvasSize) / (-v.z * 2);
            proj[i].y = (canvasSize / 2) - (v.y * canvasSize) / (-v.z * 2);
        }

        BeginDrawing();

        ClearBackground(BLACK);

        // ****** Utterlly abysmal way to test for edges *****
        for (int i = 0; i < 8; i++)
        {
            const Vec3f& v1 = vertices[i];

            for (int j = 0; j < 8; j++)
            {
                const Vec3f& v2 = vertices[j];

                int sims = 0;
                if (v1.x == v2.x) sims++;
                if (v1.y == v2.y) sims++;
                if (v1.z == v2.z) sims++;

                if (sims != 2) continue;

                DrawLineV(proj[i], proj[j], RED);
            }
        }
        // **************************************************

        EndDrawing();
    }
}