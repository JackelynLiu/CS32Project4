#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <iostream>

template<typename ValueType>
class Trie
{
public:
	Trie();
	~Trie();
	void reset();
	void insert(const std::string& key, const ValueType& value);
	std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

	// C++11 syntax for preventing copying and assignment
	Trie(const Trie&) = delete;
	Trie& operator=(const Trie&) = delete;

private:
	struct Node
	{
		Node(char key) :label(key) {}
		char label;
		std::vector<ValueType> values;
		std::vector<Node*> children;
	};
	Node* root;
};

template<typename ValueType>
Trie<ValueType>::Trie()
{
	root = new Node(' ');
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{

}

template<typename ValueType>
void Trie<ValueType>::reset()
{

}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
	Node* current = root;
	for (int pos = 0; pos < key.size(); pos++)
	{
		int curvec_size = current->children.size();
		int v = 0;
		for (; v < curvec_size; v++)
		{
			if (key[pos] == current->children[v]->label)
				break;
		}
		if (v == curvec_size)
		{
			current->children.push_back(new Node(key[pos]));
			current = current->children[curvec_size];
		}
		else current = current->children[v];
	}
	current->values.push_back(value);
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
	//for (int )
}


#endif // TRIE_INCLUDED
