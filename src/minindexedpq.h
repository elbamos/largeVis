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
  virtual bool contains(VIDX i) const = 0;
  virtual void batchInsert(VIDX n) = 0;
  virtual bool decreaseIf(VIDX i, D newKey) = 0;
  virtual D keyOf(VIDX i) const = 0;
  virtual VIDX pop() = 0;
};


// MinIndexedPQ Basedon https://github.com/kartikkukreja/blog-codes/blob/master/src/Indexed%20Min%20Priority%20Queue.cpp

template<class VIDX, class D>
class MinIndexedPQ : public PQ<VIDX, D> {
private:
	VIDX N;
	std::unique_ptr< VIDX[] > heap, index;
	std::unique_ptr< D[] > keys;

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
		keys = std::unique_ptr< D[] >(new D[NMAX + 1]);
		heap = std::unique_ptr< VIDX[] >(new VIDX[NMAX + 1]);
		index = std::unique_ptr< VIDX[] >(new VIDX[NMAX + 1]);
		for(VIDX n = 0; n != NMAX; n++) index[n] = -1;
	}

	// check if the PQ is empty
	virtual bool isEmpty() const {
		return N == 0;
	}

	// check if i is an index on the PQ
	virtual bool contains(VIDX i) const {
		return index[i] != -1;
	}

	// return the number of elements in the PQ
	VIDX size() const {
		return N;
	}

	// We can pre-load the data structure quickly because we don't need to
	// bubbleUp when they're all INFINITY.
	virtual void batchInsert(VIDX n) {
		for (VIDX i = 0; i != n; i++) {
			N++;
			index[i] = N;
			heap[N] = i;
			keys[i] = INFINITY;
		}
	}
	// associate key with index i; 0 < i < NMAX
	void insert(VIDX i, D key) {
		N++;
		index[i] = N;
		heap[N] = i;
		keys[i] = key;
		bubbleUp(N);
	}

	// returns the index associated with the minimal key
	VIDX minIndex() const {
		return heap[1];
	}

	// returns the minimal key
	virtual D minKey() const {
		return keys[heap[1]];
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
	virtual D keyOf(VIDX i) const {
		return keys[i];
	}

	// change the key associated with index i to the specified value
	void decreaseKey(VIDX i, D key)  {
		keys[i] = key;
		if (index[i] == -1 || keys[i] == -1) return;
		bubbleUp(index[i]);
	}

	virtual bool decreaseIf(VIDX i, D key) {
		if (key < keys[i]) {
			decreaseKey(i, key);
			return true;
		} else return false;
	}

	// delete the key associated with index i
	void remove(VIDX i)   {
		VIDX ind = index[i];
		swap(ind, N--);
		bubbleUp(ind);
		bubbleDown(ind);
		index[i] = -1;
	}
};




/*
* Implementation based on code found at http://www.sanfoundry.com/cpp-program-implement-pairing-heap/
*/
template<class VIDX, class D>
class PairingHeap : public PQ<VIDX, D> {
private:
  class PairNode {
  typedef std::shared_ptr< PairNode > NodePointer;
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
  };
	typedef std::shared_ptr< PairNode > NodePointer;
	NodePointer root;
  void reclaimMemory(NodePointer t) {
    if (t != NULL) {
      reclaimMemory(t->leftChild);
      reclaimMemory(t->nextSibling);
      delete t;
    }
  }
  /*
   * Internal method that is the basic operation to maintain order.
   * Links first and second together to satisfy heap order.
   * first is root of tree 1, which may not be NULL.
   *    first->nextSibling MUST be NULL on entry.
   * second is root of tree 2, which may be NULL.
   * first becomes the result of the tree merge.
   */
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

  /*
  * Internal method that implements two-pass merging.
  * firstSibling the root of the conglomerate;
  *     assumed not NULL.
  */
  NodePointer combineSiblings(NodePointer firstSibling) {
  	Rcout << "c";
    if (firstSibling->nextSibling == NULL) {
    	Rcout << "c";
    	return firstSibling;
  	}
    Rcout << "a";
    static vector< NodePointer > treeArray(5);
    int numSiblings = 0;
    Rcout << "b";
    for (; firstSibling != NULL; numSiblings++) {
      if (numSiblings == treeArray.size()) treeArray.resize(numSiblings * 2);
      treeArray[numSiblings] = firstSibling;
      firstSibling->prev->nextSibling = NULL;
      firstSibling = firstSibling->nextSibling;
    }
    Rcout << "d";
    if (numSiblings == treeArray.size()) treeArray.resize(numSiblings + 1);
    treeArray[numSiblings] = NULL;
    int i = 0;
    Rcout << "e";
    for (; i + 1 < numSiblings; i += 2) compareAndLink(treeArray[i], treeArray[i + 1]);
    Rcout << "f";
    int j = i - 2;
    Rcout << "g";
    if (j == numSiblings - 3) compareAndLink (treeArray[j], treeArray[j + 2]);
    Rcout << "h";
    for (; j >= 2; j -= 2) compareAndLink(treeArray[j - 2], treeArray[j] );
    Rcout << "c";
    return treeArray[0];
  }

  /*
  * Internal method to clone subtree.
  */
  NodePointer clone(NodePointer t) {
  	stop("called clone");
  }
//  NodePointer clone(NodePointer t) {
//    if (t == NULL) return NULL;
//    else {
//    	NodePointer p = NodePointer(new PairNode(t->element));
//      if ((p->leftChild = clone( t->leftChild)) != NULL)
//        p->leftChild->prev = p;
//      if ((p->nextSibling = clone( t->nextSibling)) != NULL)
//        p->nextSibling->prev = p;
//      return p;
//   }
//  }

protected:

  PairingHeap(PairingHeap & rhs) {
  	stop("called other constructor.");
  //  root = NULL;
  //  *this = rhs;
  }
  /*
   * Find the smallest item in the priority queue.
   * Return the smallest item, or throw Underflow if empty.
   */
  D &findMin() {
    return root->element;
  }
  /*
   * Insert item x into the priority queue, maintaining heap order.
   * Return a pointer to the node containing the new item.
   */
  NodePointer Insert(VIDX &n, D &x) {
  	NodePointer newNode = NodePointer(new PairNode(n, x));
    if (root == NULL) root = newNode;
    else compareAndLink(root, newNode);
    PointerArray.push_back(newNode);
    ContentsArray[n] = true;
    return newNode;
  }
  void makeEmpty() {
    reclaimMemory(root);
    root = NULL;
  }
  /*
   * Change the value of the item stored in the pairing heap.
   * Does nothing if newVal is larger than currently stored value.
   * p points to a node returned by insert.
   * newVal is the new value, which must be smaller
   *    than the currently stored value.
   */
   bool decreaseIf(NodePointer p, D &newVal) {
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
	PairingHeap &operator=(PairingHeap & rhs) {
		stop("it got called");
	}
//  PairingHeap &operator=(PairingHeap & rhs) {
//    if (this != &rhs) {
//      makeEmpty( );
//      root = clone(rhs.root);
//    }
//    return *this;
//  }

  VIDX MaxSize = 0;

  std::vector< NodePointer > PointerArray;
  bool* ContentsArray;

public:
  PairingHeap(VIDX N) : MaxSize{N} {
    root = NULL;
    PointerArray = std::vector< NodePointer >();
    PointerArray.reserve(N);
    ContentsArray = new bool[N];
  }
  /*
   * Destroy the leftist heap.
   */
  /*virtual ~PairingHeap() {
    makeEmpty();
    for (auto it = PointerArray.begin();
         it != PointerArray.end();
         it++) delete *it;
    delete[] ContentsArray;
  }*/
  /*
   * Remove the smallest item from the priority queue.
   * Throws Underflow if empty.
   */
  virtual VIDX pop() {
  	NodePointer oldRoot = root;
  	Rcout << "l";
    if (root->leftChild == NULL) root = NULL;
    else root = combineSiblings(root->leftChild);
    Rcout << "i";
    VIDX ret = oldRoot -> index;
    //delete oldRoot;
    ContentsArray[ret] = false;
    return ret;
  }

  virtual bool isEmpty() const {
    return root == NULL;
  }
  virtual bool contains(VIDX i) const {
    return ContentsArray[i];
  }
  double inf = INFINITY;

  virtual void batchInsert(VIDX n) {
    for (VIDX i = 0; i != n; i++) Insert(i, inf);
  };

  virtual bool decreaseIf(VIDX i, D newKey) {
    return decreaseIf(PointerArray[i], newKey);
  };

  virtual D keyOf(VIDX i) const {
    return PointerArray[i] -> element;
  };
};