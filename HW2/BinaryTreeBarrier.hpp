#include <thread>
#include <atomic>
#include <vector>
#include <cmath>
#include <iostream>

class Node {
public:
    std::atomic<int> count;
    std::atomic<int> sense;
    int size;
    
    Node* parent;
    Node* right_child;
    Node* left_child;

    Node() : count(0), sense(0), size(2), parent(nullptr), right_child(nullptr), left_child(nullptr) {}

    void barrier();
};

class BinaryTreeBarrier {
public:
    BinaryTreeBarrier (int n); // n is total number of threads
    ~BinaryTreeBarrier();
    void barrier();
private:
    void build_tree(Node* parent, int tree_depth);
    void freeTree(Node* to_delete);

    Node* root;
    int tree_depth;
    std::vector<Node*> leafs;
    std::atomic<int> threads;
};

