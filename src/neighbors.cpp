#include "neighbors.h"

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
									          const V& x_i, const vertexidxtype& j) const {
		const distancetype d = distanceFunction(x_i, data.col(j));
		heap.emplace_back(d, j);
	}

template<class M, class V>
void AnnoySearch<M, V>::addHeap(vector< std::pair<distancetype, vertexidxtype> >& heap,
              const V& x_i, const vertexidxtype& j) const {
		const distancetype d = distanceFunction(x_i, data.col(j));
		heap.emplace_back(d, j);
		push_heap(heap.begin(), heap.end(), std::less<std::pair<distancetype, vertexidxtype>>());
		if (heap.size() > K) {
			pop_heap(heap.begin(), heap.end(), std::less<std::pair<distancetype, vertexidxtype>>());
			heap.pop_back();
		}
	}

template<class M, class V>
void AnnoySearch<M, V>::addToNeighborhood(const V& x_i, const vertexidxtype& j,
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
inline void AnnoySearch<M, V>::mergeNeighbors(const list< ivec >& localNeighborhoods) {
	Neighborhood tmp;
#ifdef _OPENMP
#pragma omp critical
#endif
{
	for (auto it = localNeighborhoods.begin(); it != localNeighborhoods.end(); ++it) {
		const ivec& indices = *it;
		const auto indicesEnd = indices.end();
		for (auto it2 = indices.begin(); it2 != indicesEnd; ++it2) {
			const vertexidxtype cur = *it2;
		  Neighborhood& neighborhood = treeNeighborhoods[cur];
		  tmp.clear();
		  tmp.swap(neighborhood);
		  neighborhood.reserve(tmp.size() + indices.size() - 1);

		  auto it3 = indices.begin();
		  auto neighboriterator = tmp.begin();
		  auto tmpe = tmp.end();

		  while (neighboriterator != tmpe && it3 != indicesEnd) {
		  	vertexidxtype newone = *neighboriterator;
		  	if (newone < *it3) ++neighboriterator;
		  	else if (*it3 < newone) {
		  		if (*it3 == cur) {
		  			++it3;
		  			continue;
		  		}
		  		newone = *it3;
		  		++it3;
		  	} else {
		  		++neighboriterator; ++it3;
		  	}
		  	neighborhood.emplace_back(newone);
		  }
		  for ( ; neighboriterator != tmpe; ++neighboriterator) neighborhood.emplace_back(*neighboriterator);
		  for ( ; it3 != indicesEnd; ++it3) if (*it3 != cur) neighborhood.emplace_back(*it3);
	  }
	}
}
}

	/*
	* The key function of the annoy-trees phase.
	*/
template<class M, class V>
void AnnoySearch<M, V>::recurse(const ivec& indices, list< ivec >& localNeighborhood) {
	const vertexidxtype I = indices.n_elem;
	if (p.check_abort()) return;
	if (I < 2) stop("Tree split failure.");
	if (I <= threshold) {
		localNeighborhood.push_back(indices);
		p.increment(I);
		return;
	}

	const vec direction = hyperplane(indices);
	const distancetype middle = median(direction);
	const uvec left = find(direction > middle);
	const uvec right = find(direction <= middle);

	if (left.n_elem >= 2 && right.n_elem >= 2) {
		recurse(indices(left), localNeighborhood);
		recurse(indices(right), localNeighborhood);
	} else { // Handles the rare case where the split fails because of equidistant points
		recurse(indices.subvec(0, I / 2), localNeighborhood);
		recurse(indices.subvec(I / 2, I - 1), localNeighborhood);
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
void AnnoySearch<M, V>::trees(const unsigned int& n_trees, const unsigned int& newThreshold) {
	threshold = newThreshold;
	const ivec indices = regspace<ivec>(0, data.n_cols - 1);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int t = 0; t < n_trees; t++) if (! p.check_abort()) {
		list< ivec > local;
		recurse(indices, local);
		mergeNeighbors(local);
	}
#ifdef _OPENMP
	if (storedThreads > 0) omp_set_num_threads(storedThreads);
#endif
}

template<class M, class V>
inline void AnnoySearch<M, V>::reduceOne(const vertexidxtype& i,
                                  vector< std::pair<distancetype, vertexidxtype> >& newNeighborhood) {
	const V& x_i = data.col(i);
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

	for (auto it = newNeighborhood.begin(); it != newNeighborhood.end(); ++it, ++continueWriting) *continueWriting = it->second;
	sort(knns.begin_col(i), continueWriting);
	std::fill(continueWriting, knns.end_col(i), -1);

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
#ifdef _OPENMP
	const unsigned int dynamo = omp_get_dynamic();
	omp_set_dynamic(0);
	const vertexidxtype chunk = (N / omp_get_max_threads()) + 1;
#pragma omp parallel for
#else
	const vertexidxtype chunk = N;
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

	auto copyContinuation = knns.begin_col(i);
	auto nend = nodeHeap.end();
	for (auto it = nodeHeap.begin(); it != nend; ++it, ++copyContinuation) *copyContinuation = it->second;
	sort(knns.begin_col(i), copyContinuation);
	if (copyContinuation == knns.begin_col(i)) stop("Neighbor exploration failure.");
	std::fill(copyContinuation, knns.end_col(i), -1);
}

template<class M, class V>
void AnnoySearch<M,V>::exploreNeighborhood(const unsigned int& maxIter) {
	const kidxtype K = knns.n_rows;
	imat old_knns = imat(K,N);

	for (unsigned int T = 0; T != maxIter; ++T) if (! p.check_abort()) {
		swap(knns, old_knns);

#ifdef _OPENMP
		const unsigned int dynamo = omp_get_dynamic();
		omp_set_dynamic(0);
		const vertexidxtype chunk = (N / omp_get_max_threads()) + 1;
#pragma omp parallel for shared(old_knns)
#else
		const vertexidxtype chunk = N;
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
	const unsigned int dynamo = omp_get_dynamic();
	if (omp_get_num_threads() > 1) omp_set_dynamic(0);
	const vertexidxtype chunk = (N / omp_get_max_threads()) + 1;
#pragma omp parallel for
#else
	const vertexidxtype chunk = N;
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
	std::fill(copyContinuation, knns.end_col(i), -1);
}

template class AnnoySearch<Mat<double>, Col<double>>;
template class AnnoySearch<SpMat<double>, SpMat<double>>;
