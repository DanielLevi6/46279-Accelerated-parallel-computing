#include <memory>


struct Node{
	int data;
	std::unique_ptr<Node> next;
};
