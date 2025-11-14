#include "largeVis.h"
using namespace Rcpp;

class MinIndexedPQ {
private:
	dimidxtype N, *heap, *index;
	vertexidxtype *keys;

	void swap(dimidxtype i, dimidxtype j);
	void bubbleUp(dimidxtype k);
	void bubbleDown(dimidxtype k);

public:
	explicit MinIndexedPQ(const dimidxtype& NMAX);
	~MinIndexedPQ();
	bool isEmpty() const;
	void insert(dimidxtype i, vertexidxtype key);
	dimidxtype minIndex() const;
	vertexidxtype minKey() const;
	dimidxtype pop();
	void rotate(const vertexidxtype& key);
};
