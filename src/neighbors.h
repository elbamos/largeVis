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
using namespace Annoy;

typedef float annoy_distance;

typedef vector< vertexidxtype > Neighborhood;

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

template <typename distancemetric, typename Distance>
using LVAnnoyIndex = AnnoyIndex<vertexidxtype, distancemetric, Distance, Kiss64Random, RcppAnnoyIndexThreadPolicy>;


template<typename distancemetric, typename Distance>
class AnnoySearch {
protected:
	Neighborhood* treeNeighborhoods;
	arma::imat knns;
	void advanceHeap(MinIndexedPQ& positionHeap, vector< Position>& positionVector) const;

	inline void addHeap(vector< std::pair<distancemetric, vertexidxtype> >& heap, const vertexidxtype& i, const vertexidxtype& j) const;
	inline void addToNeighborhood(const vertexidxtype& i, const vertexidxtype& j,
                               vector< std::pair<distancemetric, vertexidxtype> >& neighborhood) const;

public:
	LVAnnoyIndex<distancemetric, Distance>* annoy_index;
	const kidxtype K;
	const vertexidxtype N;
	Progress p;
	void reduce();

	void exploreOne(const vertexidxtype& i, const arma::imat& old_knns,
                 vector< std::pair<distancemetric, vertexidxtype> >& nodeHeap,
                 MinIndexedPQ& positionHeap, vector< Position >& positionVector);
	void reduceOne(const vertexidxtype& i);
	void sortCopyOne(vector< std::pair<distancemetric, vertexidxtype>>& holder, const vertexidxtype& i);
	void exploreNeighborhood(const unsigned int& maxIter);
	arma::imat sortAndReturn();

	AnnoySearch(LVAnnoyIndex<distancemetric, Distance>* annoy_index, const vertexidxtype& N, const kidxtype& K, const bool &verbose, const int &maxIter, const long &pCount) :
		annoy_index{annoy_index},
		K{K},
		N{N},
		p(pCount, verbose) {}
};


template<class M, typename distancemetric, typename Distance>
void trees(
		LVAnnoyIndex<distancemetric, Distance>& annoy_index,
		const M& data,
		const unsigned int& n_trees,
		const Rcpp::Nullable< Rcpp::String > &savefile,
		Progress& p
	);

template<typename distancemetric, typename Distance>
class SortCopyWorker : public RcppParallel::Worker  {
public:
	AnnoySearch<distancemetric, Distance> *searcher;

	SortCopyWorker(AnnoySearch<distancemetric, Distance> *searcher) : searcher {searcher} {}

	void operator()(std::size_t begin, std::size_t end) {
		vector< std::pair<distancemetric, vertexidxtype>> holder;
		holder.reserve(searcher->K);
		for (vertexidxtype i = begin; i != end; ++i) if (searcher->p.increment()) {
			searcher->sortCopyOne(holder, i);
		}
	}
};

template<typename distancemetric, typename Distance>
class ReduceWorker : public RcppParallel::Worker  {
public:
	AnnoySearch<distancemetric, Distance> *searcher;

	ReduceWorker(AnnoySearch<distancemetric, Distance> *searcher) : searcher {searcher} {}

	void operator()(std::size_t begin, std::size_t end) {
		for (vertexidxtype i = begin; i < end; ++i) if (searcher->p.increment()) {
			searcher->reduceOne(i);
		}
	}
};

template<typename distancemetric, typename Distance>
class ExploreWorker : public RcppParallel::Worker {
public:
	AnnoySearch<distancemetric, Distance> *searcher;
	arma::imat *old_knns;

	ExploreWorker(AnnoySearch<distancemetric, Distance> *searcher, arma::imat *old_knns) : searcher {searcher}, old_knns {old_knns} {}

	void operator()(std::size_t begin, std::size_t end) {
		/*
		 * The goal here is to maintain a size-K minHeap of the points with the shortest distances
		 * to the target point. This is a merge sort with more than two sorted arrays being merged.
		 * We can use a simple priority queue because the number of entries in the queue, which equals
		 * K + 1, is small and well-controlled.
		 */
		vector< std::pair<distancemetric, vertexidxtype> > nodeHeap;
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
