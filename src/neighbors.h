#ifndef _LARGEVISNEIGHBORS
#define _LARGEVISNEIGHBORS
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

struct HeapObject {
	distancetype d;
	vertexidxtype n;
	HeapObject(distancetype d, vertexidxtype n) : d(d), n(n) {}
	bool operator<(const struct HeapObject& other) const {
		return d < other.d;
	}
};
typedef priority_queue<HeapObject> MaxHeap;
typedef vector< imat::col_iterator > PositionVector;
typedef vector< vertexidxtype > Neighborhood;

template<class M, class V>
class DistanceAdder {
protected:
	const M& data;
	const kidxtype K;
	virtual double distanceFunction(const V& x_i, const V& x_j) = 0;
	DistanceAdder(const M& data,
               const kidxtype K) :
		data{data}, K{K} {}
public:
	void add(MaxHeap& thisHeap,
          const V& x_i,
          const vertexidxtype j) {
		const distancetype d = distanceFunction(x_i, data.col(j));
	//	if (d != 0) {
			thisHeap.emplace(d, j);
			if (thisHeap.size() > K) thisHeap.pop();
	//	}
	}
};

template<class M, class V>
class AnnoySearch {
protected:
	const M& data;
	const vertexidxtype N;
	Progress& p;
	Neighborhood* treeNeighborhoods;
	set< vertexidxtype >** treeHolder;
	imat knns;

	int threshold = 0;

	std::uniform_real_distribution<double> rnd;
	std::mt19937_64 mt;

	virtual arma::vec hyperplane(const arma::ivec& indices) = 0;

	void addNeighbors(const arma::ivec& indices) {
		ivec neighbors = ivec(indices);
		Neighborhood tmpStorage = Neighborhood();
		ivec::iterator newEnd = neighbors.end();
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			for (ivec::iterator it = neighbors.begin();
		      it != newEnd;
		      it++) {
				tmpStorage.clear();
				tmpStorage.swap(treeNeighborhoods[*it]);
				treeNeighborhoods[*it].reserve(tmpStorage.size() + indices.n_elem);
				ivec::iterator newIt = neighbors.begin();
				vector< vertexidxtype >::iterator oldIt = tmpStorage.begin();
				vector< vertexidxtype >::iterator oldEnd = tmpStorage.end();
				vertexidxtype last;
				vertexidxtype best = -1;
				while (oldIt != oldEnd || newIt != newEnd) {
					if (oldIt == oldEnd) best = *newIt++;
					else if (newIt == newEnd) best = *oldIt++;
					else best = (*newIt < *oldIt) ? *newIt++ : *oldIt++;
					if (best == last || best == *it) continue;
					treeNeighborhoods[*it].push_back(best);
					last = best;
				}
			}
		}
	}

	void recurse(const arma::ivec& indices) {
		const vertexidxtype I = indices.n_elem;
		// const int D = data.n_rows;
		if (p.check_abort()) return;
		if (I < 2) stop("Tree split failure.");
		if (I <= threshold) {
			addNeighbors(indices);
			p.increment(I);
			return;
		}

		vec direction = hyperplane(indices);
		const distancetype middle = median(direction);
		const uvec left = find(direction > middle);
		const uvec right = find(direction <= middle);

		if (left.size() >= 2 && right.size() >= 2) {
			recurse(indices(left));
			recurse(indices(right));
		} else { // Handles the rare case where the split fails because of equidistant points
			recurse(indices.subvec(0, indices.size() / 2));
			recurse(indices.subvec(indices.size() / 2, indices.size() - 1));
		}
	};

	inline long sample(long i) {
		return (long) (rnd(mt) * (i - 1));
	}

	inline void copyHeapToMatrix(set< vertexidxtype >* tree,
                              const kidxtype K,
                              const vertexidxtype i) {
		set< vertexidxtype >::iterator sortIterator = tree -> begin();
		set< vertexidxtype >::iterator end = tree -> end();
		vertexidxtype j = 0;
		while (sortIterator != end) knns(j++, i) = *sortIterator++;
		if (j == 0) stop("Tree failure.");
		while (j < K) knns(j++, i) = -1;
	}

	inline void heapToSet(MaxHeap& thisHeap,
                       set< vertexidxtype >* theSet) const {
		while (! thisHeap.empty()) {
			theSet -> emplace(thisHeap.top().n);
			thisHeap.pop();
		}
	}

public:
	AnnoySearch(const M& data, Progress& p) : data{data}, N(data.n_cols), p{p} { }

	void setSeed(Rcpp::Nullable< Rcpp::NumericVector > seed) {
		long innerSeed;
		if (seed.isNotNull()) {
#ifdef _OPENMP
			omp_set_num_threads(1);
#endif
			innerSeed = NumericVector(seed)[0];
		} else {
			std::random_device hardseed;
			innerSeed = hardseed();
		}
		mt = mt19937_64(innerSeed);
	}

	void trees(const int n_trees, const int newThreshold) {
		threshold = newThreshold;
		treeNeighborhoods = new Neighborhood[N];
		for (vertexidxtype i = 0; i < N; i++) {
			treeNeighborhoods[i].push_back(i);
		}
		const ivec indices = regspace<ivec>(0, data.n_cols - 1);
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int t = 0; t < n_trees; t++) if (! p.check_abort()) {
			recurse(indices);
		}
	}

	void reduce(const kidxtype K,
             DistanceAdder<M, V>* adder) {
		treeHolder = new set< vertexidxtype >*[N];
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (vertexidxtype i = 0; i < N; i++) if (p.increment()) {
			const V x_i = data.col(i);
			MaxHeap thisHeap = MaxHeap();
			vector< vertexidxtype > neighborhood = treeNeighborhoods[i];
			for (vector< vertexidxtype >::iterator j = neighborhood.begin();
        j != neighborhood.end();
        j++)
				adder -> add(thisHeap, x_i, *j);
			treeNeighborhoods[i].clear();
			treeHolder[i] = new set< vertexidxtype >();
			heapToSet(thisHeap, treeHolder[i]);
		}
		delete[] treeNeighborhoods;
	}

	void convertToMatrix(const kidxtype K) {
		knns = imat(K,N);
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (vertexidxtype i = 0; i < N; i++) if (p.increment()) {
			copyHeapToMatrix(treeHolder[i], K, i);
			delete treeHolder[i];
		}
		delete[] treeHolder;
	}

	/*
	 * Re-sort by distance.
	 */
	arma::imat getMatrix(DistanceAdder<M, V>* adder) {
	const kidxtype K = knns.n_rows;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (vertexidxtype i = 0; i < N; i++) if (p.increment()) {
			const V x_i = data.col(i);
			MaxHeap thisHeap = MaxHeap();
			for (imat::col_iterator it = knns.begin_col(i);
        it != knns.end_col(i);
        it++) {
				adder -> add( thisHeap, x_i, *it);
			}
			vertexidxtype j = K - 1;
			while (j >= thisHeap.size()) knns(j--, i) = -1;
			while (! thisHeap.empty()) {
				knns(j--, i) = thisHeap.top().n;
				thisHeap.pop();
			}
		}
		return knns;
	}

	arma::imat exploreNeighborhood(const int maxIter,
                                  DistanceAdder<M, V>* adder) {
		const kidxtype K = knns.n_rows;
		imat old_knns  = imat(K,N);
		for (int T = 0; T < maxIter; T++) if (! p.check_abort()) {
			imat tmp = old_knns;
			old_knns = knns;
			knns = tmp;
			MaxHeap thisHeap = MaxHeap();
			set< vertexidxtype > sorter = set< vertexidxtype >();
#ifdef _OPENMP
#pragma omp parallel for shared(old_knns) private(thisHeap, sorter)
#endif
			for (vertexidxtype i = 0; i < N; i++) if (p.increment()) {
				const V x_i = data.col(i);

				PositionVector positions = PositionVector(), ends = PositionVector();
				positions.reserve(K + 1); ends.reserve(K + 1);

				positions.push_back(old_knns.begin_col(i));
				ends.push_back(old_knns.end_col(i));

				for (imat::col_iterator it = old_knns.begin_col(i);
         it != ends[0] && *it != -1;
         it++) {
					positions.push_back(old_knns.begin_col(*it));
					ends.push_back(old_knns.end_col(*it));
				}

				vertexidxtype lastOne = N + 1;
				// This is a K + 1 vector merge sort running in O(K * N)
				PositionVector::iterator theEnd = positions.end();
				while (true) {
					imat::col_iterator whch = 0;

					for (pair< PositionVector::iterator,
          PositionVector::iterator > it(positions.begin(),
                                        ends.begin());
          it.first != theEnd;
          it.first++, it.second++) while (*it.first != *it.second) { // For each neighborhood, keep going until
          	// we find a non-dupe or get to the end

          	if (**it.first == -1) advance(*it.first, distance(*it.first, *it.second));
          	else if (**it.first == i || **it.first == lastOne) advance(*it.first, 1);
          	else if (whch == 0 || **it.first < *whch) {
          		whch = *it.first;
          		break;
          	} else break;
          }
          if (whch == 0) break;
          lastOne = *whch;
          advance(whch, 1);

          adder ->add( thisHeap, x_i, lastOne);
				}
				/*
				 * Before the last iteration, we keep the matrix sorted by vertexid, which makes the merge above
				 * more efficient.  In the last iteration, sort by distance.
				 */
				if (T != maxIter - 1) {
					sorter.clear();
					heapToSet(thisHeap, &sorter);

					set< vertexidxtype >::iterator sortIterator = sorter.begin();
					vertexidxtype j = 0;
					while (sortIterator != sorter.end()) knns(j++, i) = *sortIterator++;
					if (j == 0) stop("Neighbor exploration failure.");
					while (j < K) knns(j++,i) = -1;
				} else {
					vertexidxtype j = K - 1;
					while (j >= thisHeap.size()) knns(j--, i) = -1;
					while (! thisHeap.empty()) {
						knns(j--, i) = thisHeap.top().n;
						thisHeap.pop();
					}
				}
			}
		}
		return knns;
	}
};
#endif