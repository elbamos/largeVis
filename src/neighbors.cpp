#include "neighbors.h"

template<class M, class V>
void AnnoySearch<M, V>::advanceHeap(MinIndexedPQ& positionHeap,
                                    vector< Position>& positionVector) const {
	dimidxtype whichColumn = positionHeap.minIndex();
	Position& iterators = positionVector[whichColumn];
	vertexidxtype adv = iterators.advance();
	if (adv == -1) positionHeap.pop();
	else positionHeap.rotate(adv);
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
void AnnoySearch<M, V>::mergeNeighbors(const list< Neighborholder >& localNeighborhoods) {
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	Neighborhood tmp;
	for (auto it = localNeighborhoods.begin(); it != localNeighborhoods.end(); ++it) {
		const ivec& indices = **it;
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
		  auto back = std::back_inserter(neighborhood);
		  copy(neighboriterator, tmpe, back);
		  copy_if(it3, indicesEnd, back, [&cur](const vertexidxtype& tst) {return tst != cur;});
	  }
	}
}
}

Neighborholder copyTo(const Neighborholder& indices, const uvec& selections) {
	Neighborholder out = make_shared<ivec>(selections.n_elem);
	std::transform(selections.begin(), selections.end(), out->begin(),
                [&indices](const uword& it) {return (*indices)[it];});
	return out;
}

	/*
	* The key function of the annoy-trees phase.
	*/
template<class M, class V>
void AnnoySearch<M, V>::recurse(const Neighborholder& indices, list< Neighborholder >& localNeighborhood) {
	const arma::uword I = indices->n_elem;
	if (I <= threshold) {
		localNeighborhood.emplace_back(indices);
		p.increment(I);
	} else {
		vec direction = hyperplane(*indices);
		distancetype middle = median(direction);
		uvec left = find(direction > middle);
		if (left.n_elem > (I - 2) || left.n_elem < 2) {
			direction.randu();
			middle = 0.5;
			left = find(direction > middle);
		}
		const uvec right = find(direction <= middle);
		recurse(copyTo(indices, left), localNeighborhood);
		recurse(copyTo(indices, right), localNeighborhood);
	}
}

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
	threshold2 = threshold * 4;
	Neighborholder indices = make_shared<ivec>(regspace<ivec>(0, data.n_cols - 1));
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (unsigned int t = 0; t < n_trees; ++t) if (! p.check_abort()) {
		list< Neighborholder > local;
		recurse(indices, local);
		mergeNeighbors(local);
	}
#ifdef _OPENMP
	if (storedThreads > 0) omp_set_num_threads(storedThreads);
#endif
}

template<class M, class V>
void AnnoySearch<M, V>::reduceOne(const vertexidxtype& i,
                                  vector< std::pair<distancetype, vertexidxtype> >& newNeighborhood) {
	const V& x_i = data.col(i);
	newNeighborhood.clear();
	Neighborhood& neighborhood = treeNeighborhoods[i];

	/*
	* Sort by distance the first K items, by assembling into a heap.
	*/
	for (auto j = neighborhood.begin(); j != neighborhood.end(); ++j) {
		addToNeighborhood(x_i, *j, newNeighborhood);
	}

	/*
	 * Copy the remainder (max K elements) into a column
	 * of the matrix. Sort those K by vertex index. Pad the column with -1's, if necessary.
	 */
	auto continueWriting = std::transform(newNeighborhood.begin(), newNeighborhood.end(), knns.begin_col(i),
                                        [](const std::pair<distancetype, vertexidxtype>& input) {return input.second;});
	if (continueWriting == knns.begin_col(i)) throw Rcpp::exception("At reduction, no neighbors for vertex.");
	sort(knns.begin_col(i), continueWriting);
	std::fill(continueWriting, knns.end_col(i), -1);

	treeNeighborhoods[i].resize(0);
}

	/*
	* After the annoy-tree phase, each point has a neighborhood of candidate points
	* generated by each tree.  This function finds the K-shortest-distance points
	* for each point, and copies them into the knns matrix, sorted by index.
	*/
template<class M, class V>
void AnnoySearch<M, V>::reduce() {
	knns = imat(K,N);

	ReduceWorker<M,V> worker(this);
	parallelFor(0, N, worker);
}


template<class M, class V>
void AnnoySearch<M,V>::exploreOne(const vertexidxtype& i,
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
	* more efficient.
	*
	* We can't use std:copy because we're copying from a vector of pairs
	*/
	auto copyContinuation = std::transform(nodeHeap.begin(), nodeHeap.end(), knns.begin_col(i),
                                        [](const std::pair<distancetype, vertexidxtype>& input) {return input.second;});
	if (copyContinuation == knns.begin_col(i)) throw Rcpp::exception("No neighbors after exploration - this is a bug.");
	sort(knns.begin_col(i), copyContinuation);
	std::fill(copyContinuation, knns.end_col(i), -1);
}

template<class M, class V>
void AnnoySearch<M,V>::exploreNeighborhood(const unsigned int& maxIter) {
	const kidxtype K = knns.n_rows;
	imat old_knns = imat(K,N);

	for (unsigned int T = 0; T != maxIter; ++T) if (! p.check_abort()) {
		swap(knns, old_knns);
		ExploreWorker<M, V> worker(this, &old_knns);

		parallelFor(0, N, worker);
	}
}

/*
 * Resort the matrix so in each column the neighbors are sorted by distance
 */
template<class M, class V>
imat AnnoySearch<M, V>::sortAndReturn() {
	SortCopyWorker<M, V> worker(this);
	parallelFor(0, N, worker);
	return knns;
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
		const distancetype d = distanceFunction(x_i, data.col(*it));
		holder.emplace_back(d, *it);
	}
	sort(holder.begin(), holder.end());
	auto copyContinuation = std::transform(holder.begin(), holder.end(), knns.begin_col(i),
                                        [](const std::pair<distancetype, vertexidxtype>& input) {return input.second;});
	std::fill(copyContinuation, knns.end_col(i), -1);
}

template class AnnoySearch<Mat<double>, Col<double>>;
template class AnnoySearch<SpMat<double>, SpMat<double>>;
