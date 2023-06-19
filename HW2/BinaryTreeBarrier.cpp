#include "BinaryTreeBarrier.hpp"

thread_local int thread_sense = 1;
thread_local int thread_id;

void Node::barrier() {
    if(this->parent == nullptr) return;
    int my_sense = (this->parent->sense == 0) ? 1 : 0;
    int position = ++this->parent->count;
    if(position == this->parent->size) {
        this->parent->count.store(0);
        this->parent->barrier(); // Every time just one child calls the parent's barrier
        this->parent->sense.store(my_sense);
    }
    else {
        while(this->parent->sense.load() != my_sense); // The child waits until the barrier is released from the root
    }
}

BinaryTreeBarrier::BinaryTreeBarrier(int n) {
    this->threads = 0;
    this->tree_depth = std::log2(n);
    root = new Node();
    root->parent = nullptr;
    build_tree(root, this->tree_depth);
}

BinaryTreeBarrier::~BinaryTreeBarrier() {
    freeTree(root);
}

void BinaryTreeBarrier::freeTree(Node* to_delete) {
    if(to_delete == nullptr) return;

    freeTree(to_delete->right_child);
    freeTree(to_delete->left_child);

    delete to_delete;
}

void BinaryTreeBarrier::build_tree(Node* parent, int tree_depth) {
    if(tree_depth == 0) {
        parent->right_child = nullptr;
        parent->left_child = nullptr;
        leafs.push_back(parent);
        return;
    }

    parent->right_child = new Node();
    parent->left_child = new Node();

    parent->right_child->parent = parent;
    parent->left_child->parent = parent;

    build_tree(parent->right_child, tree_depth - 1);
    build_tree(parent->left_child, tree_depth - 1);
}

void BinaryTreeBarrier::barrier() {
    thread_id = threads++;
    leafs[thread_id]->barrier();
}




