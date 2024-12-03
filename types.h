
#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
#include <iostream>
#include <stdlib.h>

typedef unsigned int ui;
typedef double decimal;

typedef uint32_t VertexID;
typedef ui LabelID;

class TreeNode {
public:
    VertexID id_;
    VertexID parent_;
    ui level_;
    ui under_level_count_;
    ui children_count_;
    ui bn_count_;
    ui fn_count_;
    VertexID* under_level_;
    VertexID* children_;
    VertexID* bn_; // back neighbor according to the bfs order
    VertexID* fn_; // front neighbor according to the bfs order
    size_t estimated_embeddings_num_;
public:
    TreeNode() {
        id_ = 0;
        under_level_ = NULL;
        bn_ = NULL;
        fn_ = NULL;
        children_ = NULL;
        parent_ = 0;
        level_ = 0;
        under_level_count_ = 0;
        children_count_ = 0;
        bn_count_ = 0;
        fn_count_ = 0;
        estimated_embeddings_num_ = 0;
    }

    ~TreeNode() {
        delete[] under_level_;
        delete[] bn_;
        delete[] fn_;
        delete[] children_;
    }

    void initialize(const ui size) {
        under_level_ = new VertexID[size];
        bn_ = new VertexID[size];
        fn_ = new VertexID[size];
        children_ = new VertexID[size];
    }

    void printTreeNode(){
        printf("TreeNode id : %d, parent : %d, level : %d\n", id_,parent_, level_);
        printf("TreeNode id : %d Children : ",id_);
        for(int i = 0; i < children_count_;i++){
            printf(" %d ",children_[i]);
        }
        printf("\n");
    }

    void printBnAndFn(){

        printf("TreeNode id : %d Back neighbors ",id_);
        for(int i = 0; i < bn_count_; i++){
            printf(" %d ",bn_[i]);
        }
        printf("\n");
        printf("TreeNode id : %d Front neighbors ",id_);
        for(int i = 0; i < fn_count_; i++){
            printf(" %d ",fn_[i]);
        }
        printf("\n");
    }
};


#endif