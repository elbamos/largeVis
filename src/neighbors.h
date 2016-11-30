#ifndef _LARGEVISNEIGHBORS
#define _LARGEVISNEIGHBORS
#include "largeVis.h"
#include <vector>
#include <memory>
#include "progress.hpp"
#include "minpq.h"

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

	void recurse(const Neighborholder& indices, list< Neighborholder >& localNeighborhood);
	void mergeNeighbors(const list< Neighborholder >& neighbors);

	void reduceOne(const vertexidxtype& i, vector< std::pair<distancetype, vertexidxtype> >& newNeighborhood);
	void reduceThread(const vertexidxtype& loopstart, const vertexidxtype& end);

	void exploreThread(const imat& old_knns, const vertexidxtype& loopstart, const vertexidxtype& end);
	void exploreOne(const vertexidxtype& i, const imat& old_knns,
                  vector< std::pair<distancetype, vertexidxtype> >& nodeHeap,
                  MinIndexedPQ& positionHeap, vector< Position >& positionVector);
	void advanceHeap(MinIndexedPQ& positionHeap, vector< Position>& positionVector) const;

	void sortCopyOne(vector< std::pair<distancetype, vertexidxtype>>& holder, const vertexidxtype& i);
	void sortCopyThread(const vertexidxtype& start, const vertexidxtype& end);

	inline void addHeap(vector< std::pair<distancetype, vertexidxtype> >& heap, const V& x_i, const vertexidxtype& j) const;
	inline void addToNeighborhood(const V& x_i, const vertexidxtype& j,
                         vector< std::pair<distancetype, vertexidxtype> >& neighborhood) const;

protected:
	const M& data;
	const kidxtype K;
	const vertexidxtype N;
	Progress& p;
	int threshold = 0;

	virtual double distanceFunction(const V& x_i, const V& x_j) const = 0;
	virtual vec hyperplane(const ivec& indices) = 0;

	inline long sample(const long& i) {
		return (long) (rnd(mt) * (i - 1));
	}

public:
	AnnoySearch(const M& data, const kidxtype& K, Progress& p) : data{data}, K{K}, N(data.n_cols), p{p} {
		treeNeighborhoods = new Neighborhood[N];
		for (vertexidxtype i = 0; i != N; ++i) treeNeighborhoods[i] = Neighborhood();
	}

	AnnoySearch(const AnnoySearch& other) : AnnoySearch(other.data, other.K, other.p) {}

	virtual ~AnnoySearch() {
		delete[] treeNeighborhoods;
	}

	void setSeed(Rcpp::Nullable< NumericVector >& seed);

	void trees(const unsigned int& n_trees, const unsigned int& newThreshold);
	void reduce();
	void exploreNeighborhood(const unsigned int& maxIter);
	imat sortAndReturn();
};
#endif
