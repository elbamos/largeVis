#include "neighbors.h"

#ifdef _OPENMP
#define LOOPSTART omp_get_thread_num()
#define CHUNK omp_get_num_threads()
#else
#define LOOPSTART 0
#define CHUNK 1
#endif

template<class M, class V>
void AnnoySearch<M, V>::advanceHeap(MinIndexedPQ& positionHeap, vector< Position>& positionVector) const {
	dimidxtype whichColumn = positionHeap.minIndex();
	Position& iterators = positionVector[whichColumn];
	vertexidxtype adv = iterators.advance();
	if (adv == -1) positionHeap.pop();
	else positionHeap.rotate(adv);
}


template<class M, class V>
void AnnoySearch<M, V>::add(vector< std::pair<distancetype, vertexidxtype> >& heap,
          const V& x_i,
          const vertexidxtype& j) const {
		const distancetype d = distanceFunction(x_i, data.col(j));
		heap.emplace_back(d, j);
	}

template<class M, class V>
void AnnoySearch<M, V>::addHeap(vector< std::pair<distancetype, vertexidxtype> >& heap,
              const V& x_i,
              const vertexidxtype& j) const {
		const distancetype d = distanceFunction(x_i, data.col(j));
		heap.emplace_back(d, j);
		push_heap(heap.begin(), heap.end(), std::less<std::pair<distancetype, vertexidxtype>>());
		if (heap.size() > K) {
			pop_heap(heap.begin(), heap.end(), std::less<std::pair<distancetype, vertexidxtype>>());
			heap.pop_back();
		}
	}

template<class M, class V>
void AnnoySearch<M, V>::addToNeighborhood(const V& x_i,
                        const vertexidxtype& j,
                        vector< std::pair<distancetype, vertexidxtype> >& neighborhood) const {
		const distancetype d = distanceFunction(x_i, data.col(j));
		neighborhood.emplace_back(d, j);
		push_heap(neighborhood.begin(), neighborhood.end(), std::less<std::pair<distancetype, vertexidxtype>>());
		if (neighborhood.size() > K) {
			pop_heap(neighborhood.begin(), neighborhood.end(), std::less<std::pair<distancetype, vertexidxtype>>());
			neighborhood.pop_back();
		}
	}

/*
 * During the annoy-tree phase, used to copy the elements of a leaf
 * into the neighborhood for each point in the leaf.
 * The neighborhood is maintained in vertex-index order.
 */
template<class M, class V>
	inline void AnnoySearch<M, V>::addNeighbors(const ivec& indices) {
#ifdef _OPENMP
#pragma omp critical
#endif
{
	static vector< vertexidxtype > tmpStorage;
	const auto indPtr = indices.memptr();
	auto end = indPtr + indices.n_elem;
	for (auto it = indPtr; it != end; ++it) {
		vertexidxtype cur = *it;
		tmpStorage.clear();
		Neighborhood& neighborhood = treeNeighborhoods[cur];
		tmpStorage.swap(neighborhood);
		neighborhood.reserve(neighborhood.size() + indices.n_elem);

		auto it2 = indPtr;
		auto tmp = tmpStorage.begin();
		auto tmpe = tmpStorage.end();
		while (it2 != end && tmp != tmpe) {
			vertexidxtype newone = *it2;
			if (newone < *tmp) ++it2;
			else if (*tmp < newone) {
				newone = *tmp;
				++tmp;
			} else {
				++tmp; ++it2;
			}
			if (newone != cur) neighborhood.emplace_back(newone);
		}
		for ( ; it2 != end; ++it2) if (*it2 != cur) neighborhood.emplace_back(*it2);
		while (tmp < tmpe) {
			neighborhood.emplace_back(*tmp);
			++tmp;
		}
		//        treeNeighborhoods[indPtr[it]] = neighborhood;
	}
}
	}

	/*
	* The key function of the annoy-trees phase.
	*/
template<class M, class V>
void AnnoySearch<M, V>::recurse(const ivec& indices) {
	const vertexidxtype I = indices.n_elem;
	if (p.check_abort()) return;
	if (I < 2) stop("Tree split failure.");
	if (I <= threshold) {
		addNeighbors(indices);
		p.increment(I);
		return;
	}

	const vec direction = hyperplane(indices);
	const distancetype middle = median(direction);
	const uvec left = find(direction > middle);
	const uvec right = find(direction <= middle);

	if (left.n_elem >= 2 && right.n_elem >= 2) {
		recurse(indices(left));
		recurse(indices(right));
	} else { // Handles the rare case where the split fails because of equidistant points
		recurse(indices.subvec(0, I / 2));
		recurse(indices.subvec(I / 2, I - 1));
	}
};

template<class M, class V>
void AnnoySearch<M, V>::setSeed(Rcpp::Nullable< NumericVector >& seed) {
	long innerSeed;
	if (seed.isNotNull()) {
#ifdef _OPENMP
		storedThreads = omp_get_max_threads();
		omp_set_num_threads(1);
		omp_set_dynamic(0);
#endif
		innerSeed = NumericVector(seed)[0];
	} else {
		random_device hardseed;
		innerSeed = hardseed();
	}
	mt = mt19937_64(innerSeed);
}

template<class M, class V>
void AnnoySearch<M, V>::trees(const int& n_trees, const int& newThreshold) {
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

template<class M, class V>
inline void AnnoySearch<M, V>::reduceOne(const vertexidxtype& i,
                                  vector< std::pair<distancetype, vertexidxtype> >& newNeighborhood) {
	const V x_i = data.col(i);
	newNeighborhood.clear();
	Neighborhood& neighborhood = treeNeighborhoods[i];

	/*
	* Sort by distance the first K items, by first assembling into a heap. using compGreater
	* transforms it from a max heap to a min heap.
	*/
	for (auto j = neighborhood.begin(); j != neighborhood.end(); ++j) {
		addToNeighborhood(x_i, *j, newNeighborhood);
	}
	sort_heap(newNeighborhood.begin(), newNeighborhood.end(), std::less<std::pair<distancetype, vertexidxtype>>());

	/*
	 * Copy the remainder (max K elements) into a column
	 * of the matrix. Sort those K by vertex index. Pad the column with -1's, if necessary.
	 */
	auto continueWriting = knns.begin_col(i);
	auto otherEnd = knns.end_col(i);

	for (auto it = newNeighborhood.begin(); it != newNeighborhood.end(); ++it, ++continueWriting) *continueWriting = it->second;
	sort(knns.begin_col(i), continueWriting);
	for (; continueWriting != otherEnd; ++continueWriting) *continueWriting = -1;

	treeNeighborhoods[i].resize(0);
}


template<class M, class V>
inline void AnnoySearch<M, V>::reduceThread(const vertexidxtype& loopstart,
                                     const vertexidxtype& end) {
	vector< std::pair<distancetype, vertexidxtype> > newNeighborhood;
	newNeighborhood.reserve(K * threshold);
	for (vertexidxtype i = loopstart; i != end; ++i) if (p.increment()) {
		reduceOne(i, newNeighborhood);
	}
}

	/*
	* After the annoy-tree phase, each point has a neighborhood of candidate points
	* generated by each tree.  This function finds the K-shortest-distance points
	* for each point, and copies them into the knns matrix, sorted by index.
	*/
template<class M, class V>
void AnnoySearch<M, V>::reduce() {
	knns = imat(K,N);
	for (auto it = knns.begin(); it != knns.end(); ++it) *it = 1e7;
#ifdef _OPENMP
	int dynamo = omp_get_dynamic();
	omp_set_dynamic(0);
	vertexidxtype chunk = (N / omp_get_max_threads()) + 1;
#pragma omp parallel for
#else
	vertexidxtype chunk = N;
#endif
	for (vertexidxtype i = 0; i <= N; i += chunk) {
		reduceThread(i, min(i + chunk, N));
	}
#ifdef _OPENMP
	omp_set_dynamic(dynamo);
#endif
}

template<class M, class V>
inline void AnnoySearch<M, V>::exploreThread(const imat& old_knns,
                                      const vertexidxtype& loopstart,
                                      const vertexidxtype& end) {
	/*
	 * The goal here is to maintain a size-K minHeap of the points with the shortest distances
	 * to the target point.
	 */
	vector< std::pair<distancetype, vertexidxtype> > nodeHeap;
	nodeHeap.reserve(K);
	MinIndexedPQ positionHeap(K + 1);
	vector< Position > positionVector;
	positionVector.reserve(K + 1);

	for (vertexidxtype i = loopstart; i != end; ++i) if (p.increment()) {
		exploreOne(i, old_knns, nodeHeap, positionHeap, positionVector);
	}
}

template<class M, class V>
inline void AnnoySearch<M,V>::exploreOne(const vertexidxtype& i,
									                const imat& old_knns,
									                vector< std::pair<distancetype, vertexidxtype> >& nodeHeap,
									                MinIndexedPQ& positionHeap,
									                vector< Position >& positionVector) {
	const V& x_i = data.col(i);

	positionVector.clear();
	nodeHeap.clear();

	positionVector.emplace_back(old_knns, i);

	positionHeap.insert(0, old_knns(0, i));
	int posVecCnt = 1;
	auto oldEnd = old_knns.end_col(i);
	for (auto it = old_knns.begin_col(i); it != oldEnd; ++it) {
		if (*it == -1) break;
		positionVector.emplace_back(old_knns, *it);
		vertexidxtype id = * (positionVector.back().first);
		positionHeap.insert(posVecCnt++, id);
	}

	vertexidxtype lastOne = -1;
	while (! positionHeap.isEmpty()) {
		const vertexidxtype nextOne = positionHeap.minKey();

		if (nextOne != lastOne && nextOne != i) {
			addHeap(nodeHeap, x_i, nextOne);
			lastOne = nextOne;
		}
		advanceHeap(positionHeap, positionVector);
	}

	/*
	* Before the last iteration, we keep the matrix sorted by vertexid, which makes the merge above
	* more efficient.  In the last iteration, sort by distance.
	*/
	sort_heap(nodeHeap.begin(), nodeHeap.end(), std::less<std::pair<distancetype, vertexidxtype>>());
	//if (! is_sorted(nodeHeap->begin(), nodeHeap->end())) stop ("Node heap not reverse sorted");
	auto copyContinuation = knns.begin_col(i);
	auto nend = nodeHeap.end();
	for (auto it = nodeHeap.begin(); it != nend; ++it, ++copyContinuation) *copyContinuation = it->second;
	sort(knns.begin_col(i), copyContinuation);
	if (copyContinuation == knns.begin_col(i)) stop("Neighbor exploration failure.");
	auto kend = knns.end_col(i);
	for ( ; copyContinuation < kend; ++copyContinuation) *copyContinuation = -1;
}

template<class M, class V>
void AnnoySearch<M,V>::exploreNeighborhood(const unsigned int& maxIter) {
	const kidxtype K = knns.n_rows;
	imat old_knns = imat(K,N);

	for (unsigned int T = 0; T != maxIter; ++T) if (! p.check_abort()) {
		swap(knns, old_knns);

#ifdef _OPENMP
		int dynamo = omp_get_dynamic();
		omp_set_dynamic(0);
		vertexidxtype chunk = (N / omp_get_max_threads()) + 1;
#pragma omp parallel for shared(old_knns)
#else
		vertexidxtype chunk = N;
#endif
		for (vertexidxtype i = 0; i <= N; i += chunk) {
			exploreThread(old_knns, i, min(i + chunk, N));
		}
#ifdef _OPENMP
		omp_set_dynamic(dynamo);
#endif
	}
}

/*
 * Resort the matrix so in each column the neighbors are sorted by distance
 */
template<class M, class V>
imat AnnoySearch<M, V>::sortAndReturn() {
#ifdef _OPENMP
	int dynamo = omp_get_dynamic();
	if (omp_get_num_threads() > 1) omp_set_dynamic(0);
	vertexidxtype chunk = (N / omp_get_max_threads()) + 1;
#pragma omp parallel for
#else
	vertexidxtype chunk = N;
#endif
	for (vertexidxtype i = 0; i <= N; i += chunk) {
		sortCopyThread(i, min(i + chunk, N));
	}
#ifdef _OPENMP
	if (omp_get_num_threads() > 1) omp_set_dynamic(dynamo);
#endif
	return knns;
}

template<class M, class V>
void AnnoySearch<M,V>::sortCopyThread(const vertexidxtype& start,
                                      const vertexidxtype& end) {
	vector< std::pair<distancetype, vertexidxtype>> holder;
	holder.reserve(K);
	for (vertexidxtype i = start; i != end; ++i) if (p.increment()) {
		sortCopyOne(holder, i);
	}
}

template<class M, class V>
void AnnoySearch<M,V>::sortCopyOne(vector< std::pair<distancetype, vertexidxtype>>& holder,
                                   const vertexidxtype& i) {
	holder.clear();
	const V& x_i = data.col(i);
	/*
	* Its cheaper to not maintain a heap and instead just sort because we'll never have more entries than we need.
	*/
	for (auto it = knns.begin_col(i); it != knns.end_col(i) && *it != -1; ++it) {
		add(holder, x_i, *it);
	}
	sort(holder.begin(), holder.end());
	auto copyContinuation = knns.begin_col(i);
	auto hend = holder.end();
	for (auto it = holder.begin(); it != hend; ++it, ++copyContinuation) *copyContinuation = it->second;
	auto kend = knns.end_col(i);
	for ( ; copyContinuation < kend; ++copyContinuation) *copyContinuation = -1;
}

template class AnnoySearch<Mat<double>, Col<double>>;
template class AnnoySearch<SpMat<double>, SpMat<double>>;
