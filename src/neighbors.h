#ifndef _LARGEVISNEIGHBORS
#define _LARGEVISNEIGHBORS
#include "largeVis.h"
#include <vector>
#include <memory>
#include <progress.hpp>
#include "minpq.h"
#include <RcppParallel.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

typedef vector< vertexidxtype > Neighborhood;
typedef shared_ptr<ivec> Neighborholder;
/*
 * Helper class for n-way merge sort
 */
class Position : public std::pair<imat::const_col_iterator, imat::const_col_iterator > {
public:
	Position(const imat& matrix, const vertexidxtype& column) :
	pair<imat::const_col_iterator, imat::const_col_iterator>(matrix.begin_col(column), matrix.end_col(column)) {}

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
template<class M, class V>
class AnnoySearch {
private:
	Neighborhood* treeNeighborhoods;
	imat knns;
	int storedThreads = 0;
	uniform_real_distribution<double> rnd;
	mt19937_64 mt;
public:
	mutex trees_mutex;

protected:
	void advanceHeap(MinIndexedPQ& positionHeap, vector< Position>& positionVector) const;

	inline void addHeap(vector< std::pair<distancetype, vertexidxtype> >& heap, const V& x_i, const vertexidxtype& j) const;
	inline void addToNeighborhood(const V& x_i, const vertexidxtype& j,
                         vector< std::pair<distancetype, vertexidxtype> >& neighborhood) const;

public:
	const M& data;
	const kidxtype K;
	const vertexidxtype N;
	Progress p;

protected:
	int threshold2 = 0;

	virtual double distanceFunction(const V& x_i, const V& x_j) const = 0;
	virtual vec hyperplane(const ivec& indices) = 0;

	inline long sample(const long& i) {
		return (long) (rnd(mt) * (i - 1));
	}

public:
	unsigned int threshold = 0;

	AnnoySearch(const M& data, const kidxtype& K, const bool &verbose, const int &maxIter, const int&n_trees) :
			trees_mutex(), data{data}, K{K},
			N(data.n_cols),
			p((N * n_trees) + (3 * N) + (N * maxIter), verbose) {
		treeNeighborhoods = new Neighborhood[N];
		for (vertexidxtype i = 0; i != N; ++i) treeNeighborhoods[i] = Neighborhood();
	}

	virtual ~AnnoySearch() {
		delete[] treeNeighborhoods;
	}

	void setSeed(Rcpp::Nullable< NumericVector >& seed);

	void trees(const unsigned int& n_trees, const unsigned int& newThreshold);
	void reduce();

	void exploreOne(const vertexidxtype& i, const imat& old_knns,
                 vector< std::pair<distancetype, vertexidxtype> >& nodeHeap,
                 MinIndexedPQ& positionHeap, vector< Position >& positionVector);

	void reduceOne(const vertexidxtype& i, vector< std::pair<distancetype, vertexidxtype> >& newNeighborhood);
	void sortCopyOne(vector< std::pair<distancetype, vertexidxtype>>& holder, const vertexidxtype& i);

	void exploreNeighborhood(const unsigned int& maxIter);

	void recurse(const Neighborholder& indices, list< Neighborholder >& localNeighborhood);
	void mergeNeighbors(const list< Neighborholder >& neighbors);

	imat sortAndReturn();
};

template<class M, class V>
class TreesWorker : public RcppParallel::Worker  {
public:
	AnnoySearch<M, V> *searcher;
	Neighborholder *indices;

	TreesWorker<M,V>(AnnoySearch<M,V> *searcher, Neighborholder *indices) : searcher {searcher}, indices{indices} {}

	void operator()(std::size_t begin, std::size_t end) {
		for (vertexidxtype i = begin; i != end; ++i) if (! searcher->p.check_abort()) {
			list< Neighborholder > local;
			lock_guard<mutex> local_mutex(searcher->trees_mutex);
			searcher->recurse(*indices, local);
			searcher->mergeNeighbors(local);
		}
	}
};

template<class M, class V>
class SortCopyWorker : public RcppParallel::Worker  {
public:
	AnnoySearch<M, V> *searcher;

	SortCopyWorker<M,V>(AnnoySearch<M,V> *searcher) : searcher {searcher} {}

	void operator()(std::size_t begin, std::size_t end) {
		vector< std::pair<distancetype, vertexidxtype>> holder;
		holder.reserve(searcher->K);
		for (vertexidxtype i = begin; i != end; ++i) if (searcher->p.increment()) {
			searcher->sortCopyOne(holder, i);
		}
	}
};

template<class M, class V>
class ReduceWorker : public RcppParallel::Worker  {
public:
	AnnoySearch<M, V> *searcher;

	ReduceWorker<M,V>(AnnoySearch<M,V> *searcher) : searcher {searcher} {}

	void operator()(std::size_t begin, std::size_t end) {
		vector< std::pair<distancetype, vertexidxtype> > newNeighborhood;
		newNeighborhood.reserve(searcher->K * searcher->threshold);
		for (vertexidxtype i = begin; i < end; ++i) if (searcher->p.increment()) {
			searcher->reduceOne(i, newNeighborhood);
		}
	}
};

template<class M, class V>
class ExploreWorker : public RcppParallel::Worker {
public:
	AnnoySearch<M, V> *searcher;
	imat *old_knns;

	ExploreWorker<M,V>(AnnoySearch<M,V> *searcher, imat *old_knns) : searcher {searcher}, old_knns {old_knns} {}

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
