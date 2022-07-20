#include "kdnode.h"

KdNode::KdNode()
{
    level = -1;
    direction = -1;
    left = NULL;
    right = NULL;
}

KdNode::~KdNode() {}
