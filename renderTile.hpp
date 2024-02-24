#pragma once

class RenderTileNode
{
public:
    RenderTileNode* next = nullptr; // Link list node

    int x = 0;
    int y = 0;
    int polyIdx = 0;
    float trivialAccept1 = 0;
    float trivialAccept2 = 0;
    float trivialAccept3 = 0;

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