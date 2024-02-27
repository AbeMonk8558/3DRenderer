#include "renderTile.hpp"


// *********************** RenderTileNode *****************************

RenderTileNode::RenderTileNode() {}

RenderTileNode::RenderTileNode(int _x, int _y, int _polyIdx) : x(_x), y(_y), polyIdx(_polyIdx) {}

RenderTileNode::~RenderTileNode()
{
    delete next;
}

bool RenderTileNode::isTriviallyAccepted() const
{
    return a1 >= 0 && a2 >= 0 && a3 >= 0;
}

// ********************************************************************

// *********************** RenderTileList *****************************

RenderTileList::RenderTileList() {}

void RenderTileList::clear()
{
    delete root;
    
    root = last = nullptr;
    size = 0;
}

void RenderTileList::extend()
{
    if (root == nullptr)
    {
        root = new RenderTileNode();
        last = root;
        return;
    }

    last->next = new RenderTileNode();
    last = last->next;

    size++;
}

// ********************************************************************

// *********************** RenderTriangle *****************************

Vec2f RenderTriangle::_dummy;

RenderTriangle::RenderTriangle() : v1(_dummy), v2(_dummy), v3(_dummy) {}

RenderTriangle::RenderTriangle(Vec2f& _v1, Vec2f& _v2, Vec2f& _v3) : v1(_v1), v2(_v2), v3(_v3) {}

// ********************************************************************



