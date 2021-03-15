#ifndef _LARGEVISNEIGHBORS
#define _LARGEVISNEIGHBORS
#include "largeVis.h"
#include <vector>
#include <memory>
#include <progress.hpp>
#include "minpq.h"
#include <RcppParallel.h>
#define ANNOYLIB_MULTITHREADED_BUILD 1
#include <RcppAnnoy.h>

using namespace Rcpp;
using namespace std;

typedef vector< vertexidxtype > Neighborhood;
typedef shared_ptr<arma::ivec> Neighborholder;
/*
 * Helper class for n-way merge sort
 */
class Position : public std::pair<arma::imat::const_col_iterator, arma::imat::const_col_iterator > {
public:
	Position(const arma::imat& matrix, const vertexidxtype& column) :
	pair<arma::imat::const_col_iterator, arma::imat::const_col_iterator>(matrix.begin_col(column), matrix.end_col(column)) {}

	vertexidxtype advance() {
		++first;
		return (first >= second) ? -1 : *first;
	}

	vertexidxtype get() const {
		return *first;
	}
};

// V is the type of arma vector e.g., vec
// M is the type of arma matrix e.g., mat, sp_mat
template<class M, class V, typename Distance>
class AnnoySearch {
private:
	Neighborhood* treeNeighborhoods;
	arma::imat knns;
	int storedThreads = 0;
	uniform_real_distribution<double> rnd;
	mt19937_64 mt;

protected:
	void advanceHeap(MinIndexedPQ& positionHeap, vector< Position>& positionVector) const;

	inline void addHeap(vector< std::pair<distancetype, vertexidxtype> >& heap, const vertexidxtype& i, const vertexidxtype& j) const;
	inline void addToNeighborhood(const vertexidxtype& i, const vertexidxtype& j,
                         vector< std::pair<distancetype, vertexidxtype> >& neighborhood) const;



public:
	const M& data;
	AnnoyIndex<vertexidxtype, distancetype, Distance, Kiss64Random, RcppAnnoyIndexThreadPolicy> annoy_index;
	const kidxtype K;
	const vertexidxtype N;
	Progress p;

public:

	AnnoySearch(const M& data, const kidxtype& K, const bool &verbose, const int &maxIter, const int&n_trees) :
			data{data},
			annoy_index(data.n_rows),
			K{K},
			N(data.n_cols),
			p((2 * N * n_trees) + (3 * N) + (N * maxIter), verbose) {	}

	void setSeed(Rcpp::Nullable< NumericVector >& seed);

	void trees(const unsigned int& n_trees);
	void reduce();

	void exploreOne(const vertexidxtype& i, const arma::imat& old_knns,
                 vector< std::pair<distancetype, vertexidxtype> >& nodeHeap,
                 MinIndexedPQ& positionHeap, vector< Position >& positionVector);

	void reduceOne(const vertexidxtype& i);
	void sortCopyOne(vector< std::pair<distancetype, vertexidxtype>>& holder, const vertexidxtype& i);

	void exploreNeighborhood(const unsigned int& maxIter);

	arma::imat sortAndReturn();
};

template<class M, class V, typename Distance>
class SortCopyWorker : public RcppParallel::Worker  {
public:
	AnnoySearch<M, V,Distance> *searcher;

	SortCopyWorker<M,V,Distance>(AnnoySearch<M,V,Distance> *searcher) : searcher {searcher} {}

	void operator()(std::size_t begin, std::size_t end) {
		vector< std::pair<distancetype, vertexidxtype>> holder;
		holder.reserve(searcher->K);
		for (vertexidxtype i = begin; i != end; ++i) if (searcher->p.increment()) {
			searcher->sortCopyOne(holder, i);
		}
	}
};

template<class M, class V, typename Distance>
class ReduceWorker : public RcppParallel::Worker  {
public:
	AnnoySearch<M,V,Distance> *searcher;

	ReduceWorker<M,V,Distance>(AnnoySearch<M,V,Distance> *searcher) : searcher {searcher} {}

	void operator()(std::size_t begin, std::size_t end) {
		for (vertexidxtype i = begin; i < end; ++i) if (searcher->p.increment()) {
			searcher->reduceOne(i);
		}
	}
};

template<class M, class V, typename Distance>
class ExploreWorker : public RcppParallel::Worker {
public:
	AnnoySearch<M, V,Distance> *searcher;
	arma::imat *old_knns;

	ExploreWorker<M,V,Distance>(AnnoySearch<M,V,Distance> *searcher, arma::imat *old_knns) : searcher {searcher}, old_knns {old_knns} {}

	void operator()(std::size_t begin, std::size_t end) {
		/*
		 * The goal here is to maintain a size-K minHeap of the points with the shortest distances
		 * to the target point. This is a merge sort with more than two sorted arrays being merged.
		 * We can use a simple priority queue because the number of entries in the queue, which equals
		 * K + 1, is small and well-controlled.
		 */
		vector< std::pair<distancetype, vertexidxtype> > nodeHeap;
		nodeHeap.reserve(searcher->K);
		MinIndexedPQ positionHeap(searcher->K + 1);
		vector< Position > positionVector;
		positionVector.reserve(searcher->K + 1);

		for (vertexidxtype i = begin; i < end; ++i) if (searcher->p.increment()) {
			searcher->exploreOne(i, *old_knns, nodeHeap, positionHeap, positionVector);
		}
	}
};

#endif
