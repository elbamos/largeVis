#ifndef _LARGEVISNEIGHBORS
#define _LARGEVISNEIGHBORS
#include "largeVis.h"
#include <vector>
#include "progress.hpp"

using namespace Rcpp;
using namespace std;
using namespace arma;

typedef vector< vertexidxtype > Neighborhood;

// V is the type of arma vector e.g., vec
// M is the type of arma matrix e.g., mat, sp_mat
template<class M, class V>
class DistanceAdder {
protected:
	const M& data;
	const kidxtype K;
	virtual double distanceFunction(const V& x_i, const V& x_j) const = 0;
	DistanceAdder(const M& data,
                const kidxtype& K) : data{data}, K{K} {}
public:
	virtual ~DistanceAdder() {}
	void add(const V& x_i,
           std::pair<distancetype, vertexidxtype>& j) const {
		j.first = distanceFunction(x_i, data.col(j.second));
	}
	void add(vector< std::pair<distancetype, vertexidxtype> >& heap,
           const V& x_i,
           const vertexidxtype& j) const {
		const distancetype d = distanceFunction(x_i, data.col(j));
		heap.emplace_back(d, j);
	}
	void addHeap(vector< std::pair<distancetype, vertexidxtype> >& heap,
          const V& x_i,
          const vertexidxtype& j) const {
		const distancetype d = distanceFunction(x_i, data.col(j));
		heap.emplace_back(d, j);
		push_heap(heap.begin(), heap.end());
		if (heap.size() > K) {
			pop_heap(heap.begin(), heap.end());
			heap.pop_back();
		}
	}
};

template<class M, class V>
class AnnoySearch {
private:
	Neighborhood* treeNeighborhoods;
	imat knns;
	int storedThreads = 0;
	uniform_real_distribution<double> rnd;
	mt19937_64 mt;
	inline static bool compareIndices(const std::pair<distancetype, vertexidxtype>& lhs,
                          	 const std::pair<distancetype, vertexidxtype>& rhs) {
		return lhs.second < rhs.second;
	}
	inline vertexidxtype getPosition(const std::pair<imat::const_col_iterator, imat::const_col_iterator>& position) const {
		return *position.first;
	}
	inline bool advancePosition(std::pair<imat::const_col_iterator, imat::const_col_iterator>& position) const {
		++position.first;
		return (position.first == position.second) || (*position.first == -1);
	}
	/*
	 * Helper class for n-way merge sort
	 */
	class Position : public std::pair<imat::const_col_iterator, imat::const_col_iterator > {
	public:
		Position(const imat& matrix, const vertexidxtype& column) :
			pair<imat::const_col_iterator, imat::const_col_iterator>(matrix.begin_col(column), matrix.end_col(column)) {}
	};

	/*
	 * Adapter for copying a vector of pairs into a vector of the
	 * second element of each pair.
	 */
	class CopyToMatrixIterator : public std::iterator<std::input_iterator_tag,
                                                        vertexidxtype,
                                                        long,
                                                        const vertexidxtype*,
                                                        vertexidxtype> {
	private:
		const vector< std::pair<distancetype, vertexidxtype> > theVector;
		int position;
	public:
		CopyToMatrixIterator(const vector< std::pair<distancetype, vertexidxtype> >& theVector, int  position) : theVector{theVector}, position{position} {}
		explicit CopyToMatrixIterator(const vector< std::pair<distancetype, vertexidxtype> >& theVector) : CopyToMatrixIterator(theVector, 0) {}
		CopyToMatrixIterator(const CopyToMatrixIterator& other) : CopyToMatrixIterator(other.theVector, other.position) {}
		CopyToMatrixIterator& operator++() {
			++position;
			return *this;
		}
		CopyToMatrixIterator operator++(int) {
			CopyToMatrixIterator tmp = CopyToMatrixIterator(*this);
			return ++(*this);
		}
		bool operator==(CopyToMatrixIterator& other) const {return other.position == position;}
		bool operator!=(CopyToMatrixIterator& other) const {return other.position != position;}
		reference operator*() const {return theVector[position].second;}
	};

protected:
	const M& data;
	const vertexidxtype N;
	Progress& p;
	int threshold = 0;

	virtual vec hyperplane(const ivec& indices) = 0;
/*
 * During the annoy-tree phase, used to copy the elements of a leaf
 * into the neighborhood for each point in the leaf.
 * The neighborhood is maintained in vertex-index order.
 */
	inline void addNeighbors(const ivec& indices) {
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			static vector< vertexidxtype > tmpStorage;
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				tmpStorage.clear();
				Neighborhood neighborhood = treeNeighborhoods[*it];
				tmpStorage.swap(neighborhood);
				neighborhood.reserve(tmpStorage.size() + indices.size());
				set_union(indices.begin(), indices.end(),
              	  tmpStorage.begin(), tmpStorage.end(),
              	  back_inserter(neighborhood));
				treeNeighborhoods[*it] = neighborhood;
			}
		}
	}

	/*
	 * The key function of the annoy-trees phase.
	 */
	void recurse(const ivec& indices) {
		const vertexidxtype I = indices.n_elem;
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

public:
	AnnoySearch(const M& data, Progress& p) : data{data}, N(data.n_cols), p{p} {
		treeNeighborhoods = new Neighborhood[N];
	}
	AnnoySearch(const AnnoySearch& other) : AnnoySearch(other.data, other.p) {}

	~AnnoySearch() {
		delete[] treeNeighborhoods;
	}

	void setSeed(Rcpp::Nullable< NumericVector > seed) {
		long innerSeed;
		if (seed.isNotNull()) {
#ifdef _OPENMP
			storedThreads = omp_get_max_threads();
			omp_set_num_threads(1);
#endif
			innerSeed = NumericVector(seed)[0];
		} else {
			random_device hardseed;
			innerSeed = hardseed();
		}
		mt = mt19937_64(innerSeed);
	}

	void trees(const int& n_trees, const int& newThreshold) {
		threshold = newThreshold;
		const ivec indices = regspace<ivec>(0, data.n_cols - 1);
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int t = 0; t < n_trees; t++) if (! p.check_abort()) {
			recurse(indices);
		}
#ifdef _OPENMP
		if (storedThreads > 0) omp_set_num_threads(storedThreads);
#endif
	}

	/*
	 * After the annoy-tree phase, each point has a neighborhood of candidate points
	 * generated by each tree.  This function finds the K-shortest-distance points
	 * for each point, and copies them into the knns matrix, sorted by index.
	 */
	void reduce(const kidxtype& K,
              DistanceAdder<M, V>  *adder) {
		knns = imat(K,N);
		vector< std::pair<distancetype, vertexidxtype> > newNeighborhood;

#ifdef _OPENMP
#pragma omp parallel for private(newNeighborhood)
#endif
		for (vertexidxtype i = 0; i < N; i++) if (p.increment()) {
			const V x_i = data.col(i);
			Neighborhood neighborhood = treeNeighborhoods[i];
			newNeighborhood.clear();
			newNeighborhood.reserve(K * threshold);
			for (auto j = neighborhood.begin(); j != neighborhood.end(); ++j) {
				newNeighborhood.emplace_back(-1, *j);
				adder -> add(x_i, *--(newNeighborhood.end()));
			}
			if (newNeighborhood.size() <= (K + 1)) sort(newNeighborhood.begin(), newNeighborhood.end());
			else {
				partial_sort(newNeighborhood.begin(), newNeighborhood.begin() + K + 1, newNeighborhood.end());
				newNeighborhood.resize(K + 1);
			}
			/*
			 * Get rid of identity element
			 */
			if (newNeighborhood[0].second != i) {
				auto pred = [&i](const std::pair<distancetype, vertexidxtype>& testMe) {
					return i == testMe.second;
				};
				auto identityElement = find_if(newNeighborhood.begin(), newNeighborhood.end(), pred);
				swap(*identityElement, newNeighborhood.front());
			}
			sort(++(newNeighborhood.begin()), newNeighborhood.end(), compareIndices);
			auto continueWriting = copy(CopyToMatrixIterator(newNeighborhood, 1),
                               		CopyToMatrixIterator(newNeighborhood, newNeighborhood.size()),
                              		knns.begin_col(i));
			while (continueWriting < knns.end_col(i)) {
				*continueWriting = -1;
				continueWriting++;
			}
			treeNeighborhoods[i].clear();
			treeNeighborhoods[i].shrink_to_fit();
		}
	}

public:
	/*
	 * Neighborhood exploration is a key innovation of the LargeVis algorithm.
	 */
	imat exploreNeighborhood(const unsigned int& maxIter,
                            DistanceAdder<M, V>* adder) {
		const kidxtype K = knns.n_rows;
		imat old_knns  = imat(K,N);

		auto cmp = [](const Position& a, const Position& b) {
			return *a.first > *b.first;
		};

		for (unsigned int T = 0; T != maxIter; T++) if (! p.check_abort()) {
			swap(knns, old_knns);

			vector< std::pair<distancetype, vertexidxtype> > thisHeap;
			thisHeap.reserve(K * K + K);
			vector< Position > positionHeap;
			positionHeap.reserve(K + 1);

#ifdef _OPENMP
#pragma omp parallel for shared(old_knns) private(thisHeap, positionHeap)
#endif
			for (vertexidxtype i = 0; i < N; i++) if (p.increment()) {
				const V x_i = data.col(i);

				positionHeap.emplace_back(old_knns, i);
				for (auto it = old_knns.begin_col(i); it != old_knns.end_col(i); ++it) {
					if (*it == -1) break;
					positionHeap.emplace_back(old_knns, *it);
				}
				make_heap(positionHeap.begin(), positionHeap.end(), cmp);

				vertexidxtype lastOne = -1;
				thisHeap.clear();

				while (! positionHeap.empty()) {
					const vertexidxtype nextOne = getPosition(positionHeap.front());

					if (nextOne > lastOne && nextOne != i) {
						adder -> addHeap(thisHeap, x_i, nextOne);
						lastOne = nextOne;
					}

					pop_heap(positionHeap.begin(), positionHeap.end(), cmp);

					if (advancePosition(positionHeap.back())) positionHeap.pop_back();
					else push_heap(positionHeap.begin(), positionHeap.end(), cmp);
				}

				/*
				 * Before the last iteration, we keep the matrix sorted by vertexid, which makes the merge above
				 * more efficient.  In the last iteration, sort by distance.
				 */
				sort(thisHeap.begin(), thisHeap.end(), compareIndices);
				auto copyContinuation = copy(CopyToMatrixIterator(thisHeap, 0), CopyToMatrixIterator(thisHeap, thisHeap.size()), knns.begin_col(i));
				if (copyContinuation == knns.begin_col(i)) stop("Neighbor exploration failure.");
				while (copyContinuation < knns.end_col(i)) {
					*copyContinuation = -1;
					++copyContinuation;
				}
			}

		}
		// Return the points sorted by distance
		vector< std::pair<distancetype, vertexidxtype>> holder;
		holder.reserve(K);
#ifdef _OPENMP
#pragma omp parallel for private(holder)
#endif
		for (vertexidxtype i = 0; i < N; i++) if (p.increment()) {
			const V x_i = data.col(i);
			holder.clear();
			for (auto it = knns.begin_col(i); it != knns.end_col(i) && *it != -1; ++it) {
				adder -> add(holder, x_i, *it);
			}
			sort(holder.begin(), holder.end());
			auto continueWriting = copy(CopyToMatrixIterator(holder, 0), CopyToMatrixIterator(holder, holder.size()), knns.begin_col(i));
			while (continueWriting < knns.end_col(i)) {
				*continueWriting = -1;
				++continueWriting;
			}
		}
		return knns;
	}
};
#endif
