#include "octtree.h"
#include "octnode.h"
#include <iostream>

OctTree::OctTree() {}

OctTree::~OctTree() {}

OctTree::OctTree(std::vector<QVector3D> _X, std::vector<QVector3D> _Y, std::vector<QVector3D> _Z)
{
    X = _X;
    Y = _Y;
    Z = _Z;
    root = NULL;

}

void OctTree::construct_balanced_tree() {
    //wenn würfel gefragt sind => hier die breite anpassen das es würfel werden
    construct_tree(X,
                   QVector3D(X.front().x(), Y.front().y(), Z.front().z()),
                   QVector3D(X.back().x(), Y.back().y(), Z.back().z()),
                   root);
}

void OctTree::construct_tree(std::vector<QVector3D> list_of_points ,QVector3D corner_1,QVector3D corner_2, OctNode* &k){
    if (k == NULL) {
        k = new OctNode();
        k->level = 0;
    }

    if( k->level >= 10) {
        return;
    }

    if(list_of_points.size() == 0) {
        return;
    }

    k->bottomRightBack =  corner_1;
    k->topLeftFront =     corner_2;

    if(list_of_points.size() == 1) {
        k->point = list_of_points[0];
        list_of_points.pop_back();
        return;
    }

    float p1_x = k->bottomRightBack.x();
    float p1_y = k->bottomRightBack.y();
    float p1_z = k->bottomRightBack.z();
    float p2_x = k->topLeftFront.x();
    float p2_y = k->topLeftFront.y();
    float p2_z = k->topLeftFront.z();
    float x_diff =  std::abs(p1_x - p2_x)/2;
    float y_diff =  std::abs(p1_y - p2_y)/2;
    float z_diff =  std::abs(p1_z - p2_z)/2;
    int level = k->level + 1;

    //prep
    OctNode *tmp = new OctNode();
    tmp->level = level;
    //#### TopLeftFront
    tmp->bottomRightBack = QVector3D(p1_x + x_diff , p1_y + y_diff , p1_z + z_diff);
    tmp->topLeftFront =    QVector3D(p2_x          , p2_y          , p2_z);
    k->children.push_back( tmp );

    //### TopRightFront
    tmp = new OctNode();
    tmp->level = level;
    tmp->bottomRightBack = QVector3D(p1_x          , p1_y + y_diff , p1_z + z_diff);
    tmp->topLeftFront =    QVector3D(p2_x - x_diff , p2_y          , p2_z);
    k->children.push_back( tmp );

    //### BottomRightFront
    tmp = new OctNode();
    tmp->level = level;
    tmp->bottomRightBack =  QVector3D(p1_x          , p1_y          , p1_z + z_diff);
    tmp->topLeftFront =     QVector3D(p2_x - x_diff , p2_y - y_diff , p2_z);
    k->children.push_back( tmp );

    //### BottomLeftFront
    tmp = new OctNode();
    tmp->level = level;
    tmp->bottomRightBack =  QVector3D(p1_x + x_diff , p1_y          , p1_z + z_diff);
    tmp->topLeftFront =     QVector3D(p2_x          , p2_y - y_diff , p2_z);
    k->children.push_back( tmp );

    //### TopLeftBack
    tmp = new OctNode();
    tmp->level = level;
    tmp->bottomRightBack =  QVector3D(p1_x + x_diff  , p1_y + y_diff , p1_z );
    tmp->topLeftFront =     QVector3D(p2_x           , p2_y          , p2_z - z_diff);
    k->children.push_back( tmp );

    //### TopRightBack
    tmp = new OctNode();
    tmp->level = level;
    tmp->bottomRightBack =  QVector3D(p1_x          , p1_y + y_diff , p1_z );
    tmp->topLeftFront =     QVector3D(p2_x - x_diff , p2_y          , p2_z - z_diff);
    k->children.push_back( tmp );

    //### BottomRightBack
    tmp = new OctNode();
    tmp->level = level;
    tmp->bottomRightBack =  QVector3D(p1_x          , p1_y          , p1_z);
    tmp->topLeftFront =     QVector3D(p2_x - x_diff , p2_y - y_diff , p2_z - z_diff);
    k->children.push_back( tmp );

    //## BottomLeftBack
    tmp = new OctNode();
    tmp->level = level;
    tmp->bottomRightBack =  QVector3D(p1_x + x_diff , p1_y          , p1_z );
    tmp->topLeftFront =     QVector3D(p2_x          , p2_y - y_diff , p2_z - z_diff);
    k->children.push_back( tmp );

    // zuweisung punkt in octtoren
    float mx = (p1_x + p2_x) / 2;
    float my = (p1_y + p2_y) / 2;
    float mz = (p1_z + p2_z) / 2;

    std::vector<QVector3D> list_of_points_TopLeftFront;
    std::vector<QVector3D> list_of_points_TopLeftBack;
    std::vector<QVector3D> list_of_points_BottomLeftFront;
    std::vector<QVector3D> list_of_points_BottomLeftBack;
    std::vector<QVector3D> list_of_points_TopRightFront;
    std::vector<QVector3D> list_of_points_TopRightBack;
    std::vector<QVector3D> list_of_points_BottomRightFront;
    std::vector<QVector3D> list_of_points_BottomRightBack;

    for (auto point : list_of_points) {

         if (point.x() <= mx) {
               if (point.y() <= my) {
                   if (point.z() <= mz){
                       list_of_points_TopLeftFront.push_back(point);
                       //TopLeftFront
                   }
                   else{
                       list_of_points_TopLeftBack.push_back(point);
                       //pos = TopLeftBottom
                   }
               }
               else {
                   if (point.z() <= mz){
                       list_of_points_BottomLeftFront.push_back(point);
                       //BottomLeftFront
                   }
                   else{
                       list_of_points_BottomLeftBack.push_back(point);
                       //BottomLeftBack
                   }
               }
           }
           else {
               if (point.y() <= my) {
                   if (point.z() <= mz){
                       list_of_points_TopRightFront.push_back(point);
                       //TopRightFront
                   }
                   else{
                       list_of_points_TopRightBack.push_back(point);
                       // TopRightBottom
                   }
               }
               else {
                   if (point.z() <= mz){
                       list_of_points_BottomRightFront.push_back(point);
                       //BottomRightFront
                   }
                   else{
                       list_of_points_BottomRightBack.push_back(point);
                       // BottomRightBack
                   }

               }
           }
    }

    //if (k->level == 0) {
        //test = list_of_points_BottomRightBack;
    //}

    construct_tree(list_of_points_BottomLeftBack,   k->children[1]->bottomRightBack, k->children[1]->topLeftFront, k->children[1]);
    construct_tree(list_of_points_BottomRightBack,  k->children[0]->bottomRightBack, k->children[0]->topLeftFront, k->children[0]);
    construct_tree(list_of_points_TopRightBack,     k->children[3]->bottomRightBack, k->children[3]->topLeftFront, k->children[3]);
    construct_tree(list_of_points_TopLeftBack,      k->children[2]->bottomRightBack, k->children[2]->topLeftFront, k->children[2]);
    construct_tree(list_of_points_BottomLeftFront,  k->children[5]->bottomRightBack, k->children[5]->topLeftFront, k->children[5]);
    construct_tree(list_of_points_BottomRightFront, k->children[4]->bottomRightBack, k->children[4]->topLeftFront, k->children[4]);
    construct_tree(list_of_points_TopRightFront,    k->children[7]->bottomRightBack, k->children[7]->topLeftFront, k->children[7]);
    construct_tree(list_of_points_TopLeftFront,     k->children[6]->bottomRightBack, k->children[6]->topLeftFront, k->children[6]);
}

std::vector<OctNode> OctTree::get_all_level_nodes(int level) {
    std::vector<OctNode> list;
    level_traversal(&list, root, level);
    return list;
}

void OctTree::level_traversal(std::vector<OctNode> *list, OctNode *node, int level) {
    if (node == NULL) {
        return;
    }

    if (node->level <= level) {
        list->push_back(*node);

        if(! node->children.empty()){
            level_traversal(list, node->children[0], level);
            level_traversal(list, node->children[1], level);
            level_traversal(list, node->children[2], level);
            level_traversal(list, node->children[3], level);
            level_traversal(list, node->children[4], level);
            level_traversal(list, node->children[5], level);
            level_traversal(list, node->children[6], level);
            level_traversal(list, node->children[7], level);
        }
    }
}
