#include "neighbors.h"

int getNumThreads() {
	char* num_threads = std::getenv("RCPP_PARALLEL_NUM_THREADS");
	if (strcmp(num_threads, "") == 0) return -1;
	else return atoi(num_threads);
}

using namespace arma;

template<class M, class V, typename Distance>
void AnnoySearch<M, V, Distance>::advanceHeap(MinIndexedPQ& positionHeap,
                                    vector< Position>& positionVector) const {
	dimidxtype whichColumn = positionHeap.minIndex();
	Position& iterators = positionVector[whichColumn];
	vertexidxtype adv = iterators.advance();
	if (adv == -1) positionHeap.pop();
	else positionHeap.rotate(adv);
}

template<class M, class V, typename Distance>
void AnnoySearch<M, V, Distance>::addHeap(vector< std::pair<annoy_distance, vertexidxtype> >& heap,
              									const vertexidxtype& i, const vertexidxtype& j) const {
		const annoy_distance d = annoy_index.get_distance(i, j);
		heap.emplace_back(d, j);
		push_heap(heap.begin(), heap.end(), std::less<std::pair<annoy_distance, vertexidxtype>>());
		if (heap.size() > K) {
			pop_heap(heap.begin(), heap.end(), std::less<std::pair<annoy_distance, vertexidxtype>>());
			heap.pop_back();
		}
	}

template<class M, class V, typename Distance>
void AnnoySearch<M, V, Distance>::addToNeighborhood(const vertexidxtype& i, const vertexidxtype& j,
									                        vector< std::pair<annoy_distance, vertexidxtype> >& neighborhood) const {
		const annoy_distance d = annoy_index.get_distance(i, j);
		neighborhood.emplace_back(d, j);
		push_heap(neighborhood.begin(), neighborhood.end(), std::less<std::pair<annoy_distance, vertexidxtype>>());
		if (neighborhood.size() > K) {
			pop_heap(neighborhood.begin(), neighborhood.end(), std::less<std::pair<annoy_distance, vertexidxtype>>());
			neighborhood.pop_back();
		}
	}

template<class M, class V, typename Distance>
void AnnoySearch<M, V, Distance>::trees(const unsigned int& n_trees, const Rcpp::Nullable< Rcpp::String > &savefile) {
	if (savefile.isNotNull()) {
		String save_file_path = String(savefile);
		std::string file_path = save_file_path.get_cstring();
		annoy_index.on_disk_build(file_path.c_str(), NULL);
	}
	vector<annoy_distance> tmp(data.n_rows);
	for (size_t i = 0; i < N; ++i) {
		copy(data.col(i).begin(), data.col(i).end(), tmp.begin());
		annoy_index.add_item(i, tmp.data());
	}
	p.increment(N * n_trees);
	int threads = getNumThreads();
	annoy_index.build(n_trees, threads);
	p.increment(N * n_trees);
}

template<class M, class V, typename Distance>
void AnnoySearch<M, V, Distance>::reduceOne(const vertexidxtype& i) {
	vector<vertexidxtype> neighbor_index;
	annoy_index.get_nns_by_item(i, K + 1, -1, &neighbor_index, NULL);

	sort(neighbor_index.begin(), neighbor_index.end());

	copy_if(neighbor_index.begin(), neighbor_index.end(),
         knns.begin_col(i),
         [&i](int x) { return x != i;});
}

	/*
	* After the annoy-tree phase, each point has a neighborhood of candidate points
	* generated by each tree.  This function finds the K-shortest-distance points
	* for each point, and copies them into the knns matrix, sorted by index.
	*/
	template<class M, class V, typename Distance>
	void AnnoySearch<M, V, Distance>::reduce() {
	knns = imat(K,N);

	ReduceWorker<M,V,Distance> worker(this);
	parallelFor(0, N, worker);
}


/*
 * Given a neighborhood for a point (which may include an index for the point as well), use a
 * constant size heap to find the K nearest neighbors of the point.
 */
template<class M, class V, typename Distance>
void AnnoySearch<M, V, Distance>::exploreOne( const vertexidxtype& i,
									                 const imat& old_knns,
									                 vector< std::pair<annoy_distance, vertexidxtype> >& nodeHeap,
									                 MinIndexedPQ& positionHeap,
									                 vector< Position >& positionVector) {
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
			addHeap(nodeHeap, i, nextOne);
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
                                        [](const std::pair<annoy_distance, vertexidxtype>& input) {return input.second;});
	if (copyContinuation == knns.begin_col(i)) throw Rcpp::exception("No neighbors after exploration - this is a bug.");
	sort(knns.begin_col(i), copyContinuation);
	std::fill(copyContinuation, knns.end_col(i), -1);
}

template<class M, class V, typename Distance>
void AnnoySearch<M, V, Distance>::exploreNeighborhood(const unsigned int& maxIter) {
	const kidxtype K = knns.n_rows;
	imat old_knns = imat(K,N);

	for (unsigned int T = 0; T != maxIter; ++T) if (! p.check_abort()) {
		swap(knns, old_knns);
		ExploreWorker<M, V, Distance> worker(this, &old_knns);

		parallelFor(0, N, worker);
	}
}

/*
 * Resort the matrix so in each column the neighbors are sorted by distance
 */
template<class M, class V, typename Distance>
imat AnnoySearch<M, V, Distance>::sortAndReturn() {
	SortCopyWorker<M, V, Distance> worker(this);
	parallelFor(0, N, worker);
	return knns;
}

template<class M, class V, typename Distance>
void AnnoySearch<M, V, Distance>::sortCopyOne(vector< std::pair<annoy_distance, vertexidxtype>>& holder,
                                   const vertexidxtype& i) {
	holder.clear();
	/*
	* Its cheaper to not maintain a heap and instead just sort because we'll never have more entries than we need.
	*/
	for (auto it = knns.begin_col(i); it != knns.end_col(i) && *it != -1; ++it) {
		const annoy_distance d = annoy_index.get_distance(i, *it);
		holder.emplace_back(d, *it);
	}
	sort(holder.begin(), holder.end());
	auto copyContinuation = std::transform(holder.begin(), holder.end(), knns.begin_col(i),
                                        [](const std::pair<annoy_distance, vertexidxtype>& input) {return input.second;});
	std::fill(copyContinuation, knns.end_col(i), -1);
}

template class AnnoySearch<Mat<double>, Col<double>, Euclidean>;
template class AnnoySearch<SpMat<double>, SpMat<double>, Euclidean>;
template class AnnoySearch<Mat<double>, Col<double>, Angular>;
template class AnnoySearch<SpMat<double>, SpMat<double>, Angular>;
