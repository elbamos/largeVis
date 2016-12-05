// [[Rcpp::plugins(cpp11)]]
#include <vector>
#include <memory>
#include <math.h>

template<class V, class D>
class PairingHeap {
private:
  struct PairNode {
    D distance = INFINITY;
    V index;
    bool present = false;

    PairNode* leftChild = nullptr;
    PairNode* nextSibling = nullptr;
    PairNode* prev = nullptr;
  };
	V sz = 0;

	typedef PairNode* NodePointer;

	NodePointer root = NULL;
	const V MaxSize;
	std::vector< PairNode > PointerArray;

	void compareAndLink(NodePointer &first, NodePointer second) {
		if (second == NULL) return;
		if (first->distance >= second->distance) {
			second->prev = first->prev;
			first->prev = second;
			first->nextSibling = second->leftChild;
			if (first->nextSibling != NULL) first->nextSibling->prev = first;
			second->leftChild = first;
			first = second;
		} else {
			second->prev = first;
			first->nextSibling = second->nextSibling;
			if (first->nextSibling != NULL) first->nextSibling->prev = first;
			second->nextSibling = first->leftChild;
			if (second->nextSibling != NULL) second->nextSibling->prev = second;
			first->leftChild = second;
		}
	}

	NodePointer combineSiblings(NodePointer firstSibling) {
		if (firstSibling->nextSibling == NULL) {
			return firstSibling;
		}
		static std::vector< NodePointer > treeArray(5);
		unsigned int numSiblings = 0;
		for (; firstSibling != NULL; numSiblings++) {
			if (numSiblings == treeArray.size()) treeArray.resize(numSiblings * 2);
			treeArray[numSiblings] = firstSibling;
			firstSibling->prev->nextSibling = NULL;
			firstSibling = firstSibling->nextSibling;
		}
		if (numSiblings == treeArray.size()) treeArray.resize(numSiblings + 1);
		treeArray[numSiblings] = NULL;
		unsigned int i = 0;
		for (; i + 1 < numSiblings; i += 2) compareAndLink(treeArray[i], treeArray[i + 1]);
		int j = i - 2;
		if (j == numSiblings - 3) compareAndLink (treeArray[j], treeArray[j + 2]);
		for (; j >= 2; j -= 2) compareAndLink(treeArray[j - 2], treeArray[j] );
		return treeArray[0];
	}

public:
	explicit PairingHeap(const V &N) : root(NULL), MaxSize{N},
												PointerArray(std::vector< PairNode >(N)) {
	}

	const V pop() {
		NodePointer oldRoot = root;
		if (root->leftChild == NULL) root = NULL;
		else root = combineSiblings(root->leftChild);
		V ret = oldRoot -> index;
		oldRoot -> present = false;
		sz--;
		return ret;
	}

	const V size() const {
		return sz;
	}

	const bool isEmpty() const {
		return root == NULL;
	}

	const bool contains(const V& i) const {
		return PointerArray[i].present;
	}

	void insert(const V &n, const D &x) {
		PointerArray[n].present = true;
		PointerArray[n].distance = x;
		PointerArray[n].index = n;
		PointerArray[n].leftChild = PointerArray[n].nextSibling = PointerArray[n].prev = nullptr;
		if (root == NULL) root = &PointerArray[n];
		else compareAndLink(root, &PointerArray[n]);
		sz++;
	}

	void batchInsert(const V& n, const V& start) {
		for (V i = 0; i != n; i++) {
			D key = (start == i) ? -1 : INFINITY;
			insert(i, key);
		}
	};

	bool decreaseIf(const V& i, const D &newDistance) {
		const D dist = PointerArray[i].distance;
		if (dist < newDistance) return false;
		NodePointer p = & PointerArray[i];
		p->distance = newDistance;
		if (p == root) return true;
		if (p->nextSibling != NULL)  p->nextSibling->prev = p->prev;
		if (p->prev->leftChild == p) p->prev->leftChild = p->nextSibling;
		else 												 p->prev->nextSibling = p->nextSibling;
		p->nextSibling = NULL;
		compareAndLink(root, p);
		return true;
	}

  D keyOf(const V& i) const {
  	return PointerArray[i].distance;
  }

  D topKey() const {
  	if (root == NULL) return INFINITY;
  	return root -> distance;
  }
};
