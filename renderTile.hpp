#pragma once

#include "math.hpp"
#include "SIMD.hpp"

class RenderTileNode
{
public:
    RenderTileNode* next = nullptr; // Link list node

    int x = 0;
    int y = 0;
    int polyIdx = 0;
    
    float a1 = 0;
    float a2 = 0;
    float a3 = 0;
    float r1 = 0;
    float r2 = 0;
    float r3 = 0;

    RenderTileNode();

    RenderTileNode(int _x, int _y, int _polyIdx);

    ~RenderTileNode();

    bool isTriviallyAccepted() const;
};

class RenderTileList
{
public:
    RenderTileNode* root = nullptr;
    RenderTileNode* last = nullptr;
    int size = 0;

    RenderTileList();

    void clear();

    void extend();
};

class RenderTriangle
{
public:
    Vec2f& v1;
    Vec2f& v2;
    Vec2f& v3;
    float wSign = 1;

    Vec2f rOff1; // Trivial reject offset
    Vec2f aOff1; // Trivial accept offset
    Vec2f rOff2;
    Vec2f aOff2;
    Vec2f rOff3;
    Vec2f aOff3;

    simd::float_m256 e1;
    simd::float_m256 e2;
    simd::float_m256 e3;

    RenderTriangle();
    RenderTriangle(Vec2f& _v1, Vec2f& _v2, Vec2f& _v3);

private:
    static Vec2f _dummy;
};