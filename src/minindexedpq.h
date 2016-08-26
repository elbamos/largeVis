// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
#include "largeVis.h"
//#define DEBUG
//#define DEBUG2
//#define DEBUG3

template<class VIDX, class D>
class PQ {
public:
  virtual bool isEmpty() const = 0;
  virtual bool contains(const VIDX& i) const = 0;
  virtual void batchInsert(const VIDX& n, const VIDX& start) = 0;
  virtual bool decreaseIf(const VIDX& i, const D& newKey) = 0;
  virtual D keyOf(const VIDX& i) const = 0;
  virtual VIDX pop() = 0;
};

// MinIndexedPQ Basedon https://github.com/kartikkukreja/blog-codes/blob/master/src/Indexed%20Min%20Priority%20Queue.cpp

template<class VIDX, class D>
class MinIndexedPQ : public PQ<VIDX, D> {
private:
	VIDX N;
	unique_ptr< VIDX[] > heap, index;
	unique_ptr< D[] > keys;

	void swap(VIDX i, VIDX j) {
		VIDX t = heap[i];
		heap[i] = heap[j];
		heap[j] = t;
		index[heap[i]] = i;
		index[heap[j]] = j;
	}

	void bubbleUp(VIDX k)    {
		while(k > 1 && keys[heap[k/2]] > keys[heap[k]])   {
			swap(k, k/2);
			k = k/2;
		}
	}

	void bubbleDown(VIDX k)  {
		VIDX j;
		while(2*k <= N) {
			j = 2*k;
			if(j < N && keys[heap[j]] > keys[heap[j+1]])
				j++;
			if(keys[heap[k]] <= keys[heap[j]])
				break;
			swap(k, j);
			k = j;
		}
	}

public:
	// Create an empty MinIndexedPQ which can contain atmost NMAX elements
	MinIndexedPQ(VIDX NMAX) : N(0) {
		keys = unique_ptr< D[] >(new D[NMAX + 1]);
		heap = unique_ptr< VIDX[] >(new VIDX[NMAX + 1]);
		index = unique_ptr< VIDX[] >(new VIDX[NMAX + 1]);
		for(VIDX n = 0; n != NMAX; n++) index[n] = -1;
	}

	// check if the PQ is empty
	virtual bool isEmpty() const {
		return N == 0;
	}

	// check if i is an index on the PQ
	virtual bool contains(const VIDX& i) const {
		return index[i] != -1;
	}

	// return the number of elements in the PQ
	VIDX size() const {
		return N;
	}

	// We can pre-load the data structure quickly because we don't need to
	// bubbleUp when they're all INFINITY.
	virtual void batchInsert(const VIDX& n, const VIDX& start) {
		for (VIDX i = 0; i != n; i++) {
			D key = (i == start) ? -1 : INFINITY;
			N++;
			index[i] = N;
			heap[N] = i;
			keys[i] = key;
		}
	}

	// delete the minimal key and return its associated index
	// Warning: Don't try to read from this index after calling this function
	VIDX pop() {
		VIDX min = heap[1];
		swap(1, N--);
		bubbleDown(1);
		index[min] = -1;
		heap[N+1] = -1;
		return min;
	}

	// returns the key associated with index i
	virtual D keyOf(const VIDX& i) const {
		return keys[i];
	}

	// change the key associated with index i to the specified value
	void decreaseKey(const VIDX& i, const D &key)  {
		keys[i] = key;
		if (index[i] == -1 || keys[i] == -1) return;
		bubbleUp(index[i]);
	}

	virtual bool decreaseIf(const VIDX& i, const D& key) {
		if (key < keys[i]) {
			decreaseKey(i, key);
			return true;
		} else return false;
	}

};

template<class VIDX, class D>
class PairingHeap : public PQ<VIDX, D> {
private:
  class PairNode {
  typedef shared_ptr< PairNode > NodePointer;
  public:
    D element;
    VIDX index;
    NodePointer leftChild;
    NodePointer nextSibling;
    NodePointer prev;
    PairNode(VIDX index, D element) : element{element}, index{index} {
      leftChild = NULL;
      nextSibling = NULL;
      prev = NULL;
    }
    void reclaimMemory() {
    	leftChild = nextSibling = prev = NULL;
    }
  };

	typedef shared_ptr< PairNode > NodePointer;
	NodePointer root;

  void reclaimMemory(NodePointer t) {
  	t -> reclaimMemory();
  }

  void compareAndLink(NodePointer &first, NodePointer second) {
    if (second == NULL) return;
    if (second->element < first->element) {
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
    static vector< NodePointer > treeArray(5);
    int numSiblings = 0;
    for (; firstSibling != NULL; numSiblings++) {
      if (numSiblings == treeArray.size()) treeArray.resize(numSiblings * 2);
      treeArray[numSiblings] = firstSibling;
      firstSibling->prev->nextSibling = NULL;
      firstSibling = firstSibling->nextSibling;
    }
    if (numSiblings == treeArray.size()) treeArray.resize(numSiblings + 1);
    treeArray[numSiblings] = NULL;
    int i = 0;
    for (; i + 1 < numSiblings; i += 2) compareAndLink(treeArray[i], treeArray[i + 1]);
    int j = i - 2;
    if (j == numSiblings - 3) compareAndLink (treeArray[j], treeArray[j + 2]);
    for (; j >= 2; j -= 2) compareAndLink(treeArray[j - 2], treeArray[j] );
    return treeArray[0];
  }

  NodePointer Insert(VIDX &n, D &x) {
  	NodePointer newNode = make_shared< PairNode >(n, x);
    if (root == NULL) root = newNode;
    else compareAndLink(root, newNode);
    PointerArray.push_back(newNode);
    ContentsArray[n] = true;
    return newNode;
  }

  void makeEmpty() {
    root = NULL;
    for (VIDX n = 0; n != PointerArray.size(); n++) {
    	PointerArray[n] -> reclaimMemory();
    	PointerArray[n] = NULL;
    }
  }

   bool decreaseIf(NodePointer& p, const D &newVal) {
    if (p->element < newVal) return false;
    p->element = newVal;
    if (p != root) {
      if (p->nextSibling != NULL)
        p->nextSibling->prev = p->prev;
      if (p->prev->leftChild == p)
        p->prev->leftChild = p->nextSibling;
      else
        p->prev->nextSibling = p->nextSibling;
      p->nextSibling = NULL;
      compareAndLink(root, p);
    }
    return true;
  }

  VIDX MaxSize = 0;

  vector< NodePointer > PointerArray;
  bool* ContentsArray;

public:
  PairingHeap(VIDX N) : MaxSize{N} {
    root = NULL;
    PointerArray = vector< NodePointer >();
    PointerArray.reserve(N);
    ContentsArray = new bool[N];
  }
  /*
   * Destroy the leftist heap.
   */
  virtual ~PairingHeap() {
    makeEmpty();
//    for (auto it = PointerArray.begin();
//         it != PointerArray.end();
//         it++) delete *it;
    delete[] ContentsArray;
  }
  /*
   * Remove the smallest item from the priority queue.
   * Throws Underflow if empty.
   */
  virtual VIDX pop() {
  	NodePointer oldRoot = root;
    if (root->leftChild == NULL) root = NULL;
    else root = combineSiblings(root->leftChild);
    VIDX ret = oldRoot -> index;
    //delete oldRoot;
    ContentsArray[ret] = false;
    return ret;
  }

  virtual bool isEmpty() const {
    return root == NULL;
  }
  virtual bool contains(const VIDX& i) const {
    return ContentsArray[i];
  }

  virtual void batchInsert(const VIDX& n, const VIDX& start) {
    for (VIDX i = 0; i != n; i++) {
    	D key = (start == i) ? -1 : INFINITY;
    	Insert(i, key);
    }
  };

  virtual bool decreaseIf(const VIDX& i, const D& newKey) {
    return decreaseIf(PointerArray[i], newKey);
  };

  virtual D keyOf(const VIDX& i) const {
    return PointerArray[i] -> element;
  };
};