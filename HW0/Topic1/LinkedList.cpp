#include "LinkedList.h"


int main()
{
	Node* root=new Node{6, nullptr};
	root->next = std::unique_ptr<Node>(new Node{5, nullptr});
	delete root;
	delete root; // Just for case 2	
	return 0;
}
