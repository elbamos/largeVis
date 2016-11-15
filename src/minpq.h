#include "largeVis.h"
using namespace Rcpp;

class MinIndexedPQ {
private:
	dimidxtype N, *heap, *index;
	vertexidxtype *keys;

	inline void swap(dimidxtype i, dimidxtype j) {
		std::swap(heap[i], heap[j]);
		index[heap[i]] = i;
		index[heap[j]] = j;
	}

	void bubbleUp(dimidxtype k) {
		while(k > 1 && keys[heap[k/2]] > keys[heap[k]])   {
			swap(k, k/2);
			k = k/2;
		}
	}

	void bubbleDown(dimidxtype k) {
		dimidxtype j;
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
	explicit MinIndexedPQ(dimidxtype NMAX) : N(0) {
		heap = new dimidxtype[NMAX + 1];
		index = new dimidxtype[NMAX + 1];
		keys = new vertexidxtype[NMAX + 1];
		for (dimidxtype i = 0; i != (NMAX + 1); ++i) {
			index[i] = -1;
			keys[i] = NA_INTEGER;
		}
	}

	~MinIndexedPQ() {
		delete [] keys;
		delete [] heap;
		delete [] index;
	}

	bool isEmpty() const {
		return N == 0;
	}

	void insert(dimidxtype i, vertexidxtype key) {
		++N;
		index[i] = N;
		heap[N] = i;
		keys[i] = key;
		bubbleUp(N);
	}

	dimidxtype minIndex() const {
		return heap[1];
	}

	vertexidxtype minKey() const {
		return keys[heap[1]];
	}

	dimidxtype pop() {
		const dimidxtype min = heap[1];
		swap(1, N--);
		bubbleDown(1);
		index[min] = -1;
		heap[N+1] = -1;
		return min;
	}

	void rotate(vertexidxtype key) {
		dimidxtype i = heap[1];
		keys[i] = key;
		bubbleDown(index[i]);
	}
};