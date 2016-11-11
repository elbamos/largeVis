#ifndef _LARGEVISNEIGHBORS
#define _LARGEVISNEIGHBORS
#include "largeVis.h"
#include <vector>
#include "progress.hpp"
#include "minpq.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

typedef vector< vertexidxtype > Neighborhood;

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

	inline void reduceOne(const vertexidxtype& i,
                 vector< std::pair<distancetype, vertexidxtype> >& newNeighborhood);

	inline void reduceThread(const vertexidxtype& loopstart, const vertexidxtype& end);

	inline void exploreThread(const imat& old_knns, const vertexidxtype& loopstart, const vertexidxtype& end);

	inline void exploreOne(const vertexidxtype& i, const imat& old_knns,
                  vector< std::pair<distancetype, vertexidxtype> >& nodeHeap,
                  MinIndexedPQ& positionHeap,
                  vector< Position >& positionVector);

	inline void sortCopyOne(vector< std::pair<distancetype, vertexidxtype>>& holder, const vertexidxtype& i);
	inline void sortCopyThread(const vertexidxtype& start, const vertexidxtype& end);

	inline void add(vector< std::pair<distancetype, vertexidxtype> >& heap,
          const V& x_i, const vertexidxtype& j) const;

	inline void addHeap(vector< std::pair<distancetype, vertexidxtype> >& heap,
              const V& x_i, const vertexidxtype& j) const;

	inline void addToNeighborhood(const V& x_i, const vertexidxtype& j,
                         vector< std::pair<distancetype, vertexidxtype> >& neighborhood) const;

	inline void advanceHeap(MinIndexedPQ& positionHeap, vector< Position>& positionVector) const;
	void recurse(const ivec& indices);
	inline void addNeighbors(const ivec& indices);

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
	}

	AnnoySearch(const AnnoySearch& other) : AnnoySearch(other.data, other.K, other.p) {}

	virtual ~AnnoySearch() {
		delete[] treeNeighborhoods;
	}

	void setSeed(Rcpp::Nullable< NumericVector >& seed);

	void trees(const int& n_trees, const int& newThreshold);
	void reduce();
	void exploreNeighborhood(const unsigned int& maxIter);
	imat sortAndReturn();
};
#endif
