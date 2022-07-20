#include "kdtree.h"

KdTree::KdTree(std::vector<QVector3D> _X, std::vector<QVector3D> _Y, std::vector<QVector3D> _Z)
{
    X = _X;
    Y = _Y;
    Z = _Z;
    root = NULL;
}

KdTree::KdTree() {}

KdTree::~KdTree() {}

void KdTree::construct_balanced_tree() {
    construct(0, X.size()-1, 0, root, 0);
}

void KdTree::construct(int l_index, int r_index, int level, KdNode* &k, int direction) {
    if (l_index <= r_index) {
        int m = (l_index + r_index) / 2;

        if (k == NULL) {
            k = new KdNode();
        }
        //std::cout << "direction " << direction << "\n";
        k->level = level;
        k->direction = direction;
        if (direction == 0) { // 0 = X
            k->point = X[m];
            partition_field(Y, l_index, r_index, m, 1);
            partition_field(Z, l_index, r_index, m, 2);
        }
        if (direction == 1){ // 1 = Y
            k->point = Y[m];
            partition_field(X, l_index, r_index, m, 0);
            partition_field(Z, l_index, r_index, m, 2);
        }
        if (direction == 2){ // 2 = Z
            k->point = Z[m];
            partition_field(X, l_index, r_index, m, 0);
            partition_field(Y, l_index, r_index, m, 1);
        }
       int next_direction = increment_direction(direction);

       construct(l_index, m-1, level+1, k->left, next_direction);
       construct(m+1, r_index, level+1, k->right, next_direction);
    }
}

void KdTree::partition_field(std::vector<QVector3D> &list, int l, int r, int m, int direction){
    std::vector<QVector3D> tmp1, tmp2;
    QVector3D median_point;

    // Debug purposes
    int cnt1 = 0, cnt2 = 0;

    if (direction == 0) {
        median_point = Y[m];

        for (int i = l; i <= r; i++) {
            if (list[i].y() < median_point.y()) {
                tmp1.push_back(list[i]);
                cnt1++;
            }
            if (list[i].y() > median_point.y()) {
                tmp2.push_back(list[i]);
                cnt2++;
            }
        }
    }

    if (direction == 1) {
        median_point = Z[m];

        for (int i = l; i <= r; i++) {
            if (list[i].z() < median_point.z()) {
                tmp1.push_back(list[i]);
                cnt1++;
            }
            if (list[i].z() > median_point.z()) {
                tmp2.push_back(list[i]);
                cnt2++;
            }
        }
    }

    if (direction == 2) {
        median_point = X[m];

        for (int i = l; i <= r; i++) {
            if (list[i].x() < median_point.x()) {
                tmp1.push_back(list[i]);
                cnt1++;
            }
            if (list[i].x() > median_point.x()) {
                tmp2.push_back(list[i]);
                cnt2++;
            }
        }
    }

    list[m] = median_point;
    for (int i = 0; i<tmp1.size(); i++) list[l+i] = tmp1[i];
    // TODO still crashes when choosing other starting direction, find out why
    for (int i = 0; i<tmp2.size(); i++) {
        if (m+i+1 >= X.size()) {
            std::cout << "Something went wrong" << "\n";
        }
        list[m+i+1] = tmp2[i];
    }
}

int KdTree::increment_direction(int pp_direction){
    int p_direction = pp_direction;
    p_direction += 1;
    if (p_direction > 2)
        return 0;
    return p_direction;
}

std::vector<KdNode> KdTree::get_all_level_nodes(int level) {
    std::vector<KdNode> list;
    level_traversal(&list, root, level);
    return list;
}

void KdTree::level_traversal(std::vector<KdNode> *list, KdNode *node, int level) {
    if (node == NULL) {
        return;
    }

    if (node->level <= level) {
        list->push_back(*node);
        level_traversal(list, node->left, level);
        level_traversal(list, node->right, level);
    }
}
