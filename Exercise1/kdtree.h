#ifndef KDTREE_H
#define KDTREE_H

#include "kdnode.h"
#include <vector>
#include <iostream>
#include <math.h>

class KdTree
{
public:
    KdTree(std::vector<QVector3D> _X, std::vector<QVector3D> _Y, std::vector<QVector3D> _Z);
    KdTree();
    ~KdTree();

    KdNode* root;
    std::vector<QVector3D> X;
    std::vector<QVector3D> Y;
    std::vector<QVector3D> Z;

    void construct_balanced_tree();
    void construct(int l_index, int r_index, int level, KdNode* &k, int direction);
    int  increment_direction(int direction);
    void partition_field(std::vector<QVector3D> &list, int l, int r, int m, int direction);
    std::vector<KdNode> get_all_level_nodes(int level);
    void level_traversal(std::vector<KdNode> *list, KdNode *node, int level);
};

#endif // KDTREE_H
