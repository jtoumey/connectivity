#ifndef KDTREE_H
#define KDTREE_H

const int k = 2;

template<typename T> class KdPoint
{
public:
    T coordinates[k];
};

template<typename T, typename S> class KdNode
{
private:

    KdNode* child_left;
    KdNode* child_right;

    KdPoint<T> point;

public:

    KdNode();
    int add();

};

template<typename T> class KdTree
{
private:
    S KdNode root;

public:
    KdTree();

};


#endif // KDTREE_H