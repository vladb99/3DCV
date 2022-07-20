#include "octnode.h"

OctNode::OctNode()
{
    level = -1;
    //topLeftFront = NULL;
    //bottomRightBack = NULL;
    // garantiert das 1. immer TopLeftFront und 8. BottomLeftBack
    //children.push_back(NULL); // TopLeftFront
   // children.push_back(NULL); // TopRightFront
    //children.push_back(NULL); // BottomRightFront
    //children.push_back(NULL); // BottomLeftFront
    //children.push_back(NULL); // TopLeftBack
    //children.push_back(NULL); // TopRightBack
   //children.push_back(NULL); // BottomRightBack
   // children.push_back(NULL); // BottomLeftBack


}

OctNode::~OctNode() {}
