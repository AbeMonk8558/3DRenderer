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
    return trivialAccept1 >= 0 && trivialAccept2 >= 0 && trivialAccept3 >= 0;
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

