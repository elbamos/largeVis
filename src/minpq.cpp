#include "minpq.h"

void MinIndexedPQ::swap(dimidxtype i, dimidxtype j) {
	std::swap(heap[i], heap[j]);
	index[heap[i]] = i;
	index[heap[j]] = j;
}

void MinIndexedPQ::bubbleUp(dimidxtype k) {
	while(k > 1 && keys[heap[k/2]] > keys[heap[k]])   {
		swap(k, k/2);
		k = k/2;
	}
}

void MinIndexedPQ::bubbleDown(dimidxtype k) {
	while(2*k <= N) {
		dimidxtype j = 2*k;
		if(j < N && keys[heap[j]] > keys[heap[j+1]]) ++j;
		if(keys[heap[k]] <= keys[heap[j]]) break;
		swap(k, j);
		k = j;
	}
}

MinIndexedPQ::MinIndexedPQ(const dimidxtype& NMAX) : N(0) {
	const dimidxtype nn = NMAX + 1;
	heap = new dimidxtype[nn];
	index = new dimidxtype[nn];
	keys = new vertexidxtype[nn];
	std::fill(index, index + nn, -1);
	std::fill(keys, keys + nn, NA_INTEGER);
}

MinIndexedPQ::~MinIndexedPQ() {
	delete [] keys;
	delete [] heap;
	delete [] index;
}

bool MinIndexedPQ::isEmpty() const {
	return N == 0;
}

void MinIndexedPQ::insert(dimidxtype i, vertexidxtype key) {
	++N;
	index[i] = N;
	heap[N] = i;
	keys[i] = key;
	bubbleUp(N);
}

dimidxtype MinIndexedPQ::minIndex() const {
	return heap[1];
}

vertexidxtype MinIndexedPQ::minKey() const {
	return keys[heap[1]];
}

dimidxtype MinIndexedPQ::pop() {
	const dimidxtype min = heap[1];
	swap(1, N--);
	bubbleDown(1);
	index[min] = -1;
	heap[N+1] = -1;
	return min;
}

void MinIndexedPQ::rotate(const vertexidxtype& key) {
	dimidxtype i = heap[1];
	keys[i] = key;
	bubbleDown(index[i]);
}
