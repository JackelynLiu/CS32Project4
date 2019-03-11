#include "provided.h"
#include "Trie.h"

int main()
{
	Trie<int> trie;
	//trie.insert("GATTACA", 42);		//	GATTACA	� {42}
	trie.insert("GATTACA", 17);		//	GATTACA	� {42,	17}
	//trie.insert("GATTACA", 42);		//	GATTACA	� {42,	17,	42}
	//trie.insert("GCTTACA", 30);		// GCTTACA � {30}
}

