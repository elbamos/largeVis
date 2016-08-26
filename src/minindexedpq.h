// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
#include "largeVis.h"
//#define DEBUG
//#define DEBUG2
// MinIndexedPQ Basedon https://github.com/kartikkukreja/blog-codes/blob/master/src/Indexed%20Min%20Priority%20Queue.cpp

template<class VIDX, class D>
class MinIndexedPQ {
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
	bool isEmpty() const {
		return N == 0;
	}

	// check if i is an index on the PQ
	bool contains(VIDX i) const {
		return index[i] != -1;
	}

	// return the number of elements in the PQ
	VIDX size() const {
		return N;
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
	D minKey() const {
		return keys[heap[1]];
	}

	// delete the minimal key and return its associated index
	// Warning: Don't try to read from this index after calling this function
	VIDX deleteMin() {
		VIDX min = heap[1];
		swap(1, N--);
		bubbleDown(1);
		index[min] = -1;
		heap[N+1] = -1;
		return min;
	}

	// returns the key associated with index i
	D keyOf(VIDX i) const {
		return keys[i];
	}

	// change the key associated with index i to the specified value
	void decreaseKey(VIDX i, D key)  {
		keys[i] = key;
		if (index[i] == -1 || keys[i] == -1) return;
		bubbleUp(index[i]);
	}

	bool decreaseIf(VIDX i, D key) {
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