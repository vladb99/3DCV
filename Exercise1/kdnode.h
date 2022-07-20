#ifndef KDNODE_H
#define KDNODE_H

#include <QVector3D>

class KdNode
{
public:
    KdNode();
    ~KdNode();

    QVector3D point;
    int level;
    int direction;
    KdNode* left;
    KdNode* right;
};

#endif // KDNODE_H
