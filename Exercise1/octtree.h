#ifndef OCTTREE_H
#define OCTTREE_H

#include "octnode.h"
#include <QVector3D>
#include <vector>

class OctTree
{
public:
    OctTree(std::vector<QVector3D> _X, std::vector<QVector3D> _Y, std::vector<QVector3D> _Z);
    OctTree();
    ~OctTree();

    OctNode* root;
    std::vector<QVector3D> X;
    std::vector<QVector3D> Y;
    std::vector<QVector3D> Z;
    std::vector<QVector3D> test;

    void generate_tree();
    void construct_tree(std::vector<QVector3D> list_of_points, QVector3D corner_1,QVector3D corner_2, OctNode* &k);
    void construct_balanced_tree();
    std::vector<OctNode> get_all_level_nodes(int level);
    void level_traversal(std::vector<OctNode> *list, OctNode *node, int level);
};

#endif // OCTTREE_H
