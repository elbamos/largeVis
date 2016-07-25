// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

void searchTree(const int& threshold,
                const arma::ivec& indices,
                const sp_mat& data,
                Neighborhood* heap[],
                Progress& progress) {
  const int I = indices.size();
  // const int D = data.n_rows;
  if (progress.check_abort()) return;
  if (I < 2) stop("Tree split failure.");
  if (I <= threshold) {
    addNeighbors(indices, heap, I);
    progress.increment(I);
    return;
  }
  vec direction = vec(indices.size());
  {
    int x1idx, x2idx;
    sp_mat v;
    sp_mat m;
    do {
      const vec selections = randu(2) * (I - 1);
      x1idx = indices[selections[0]];
      x2idx = indices[selections[1]];
      if (x1idx == x2idx) x2idx = indices[((int)selections[1] + 1) % indices.size()];
      const SpSubview<double> x2 = data.col(x2idx);
      const SpSubview<double> x1 = data.col(x1idx);
      // Get hyperplane
      m =  (x1 + x2) / 2; // Base point of hyperplane
      const sp_mat d = x1 - x2;
      const double dn = as_scalar(norm(d, 2));
      v =  d / dn; // unit vector
    } while (x1idx == x2idx);

    for (int i = 0; i < indices.size(); i++) {
      const int I = indices[i];
      const SpSubview<double> X = data.col(I);
      direction[i] = dot((X - m), v);
    }
  }
  // Normalize direction
  const double middle = median(direction);

  const uvec left = find(direction > middle);
  const uvec right = find(direction <= middle);
  if (left.size() >= 2 && right.size() >= 2) {
    searchTree(threshold, indices(left), data, heap, progress);
    searchTree(threshold, indices(right), data, heap, progress);
  } else {
    searchTree(threshold, indices.subvec(0, indices.size() / 2), data, heap, progress);
    searchTree(threshold, indices.subvec(indices.size() / 2, indices.size() - 1), data, heap, progress);
  }
};


arma::mat searchTreesSparse(const int& threshold,
                            const int& n_trees,
                            const int& K,
                            const int& maxIter,
                            const sp_mat& data,
                            const std::string& distMethod,
                            bool verbose) {

  const int N = data.n_cols;

  double (*distanceFunction)(const sp_mat& x_i, const sp_mat& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = sparseRelDist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = sparseCosDist;
  else distanceFunction = sparseRelDist;

  Progress p((N * n_trees) + (N) + (N * maxIter), verbose);

  Neighborhood** treeNeighborhoods = createNeighborhood(N);

  { // Artificial scope to destroy indices
  	sp_mat dataMat;
  	if (distMethod.compare(std::string("Cosine")) == 0) {
  		dataMat = sp_mat(data);
  		for (int d = 0; d < dataMat.n_cols; d++) dataMat.col(d) /= norm(dataMat.col(d));
  	} else {
  		dataMat = data;
  	}
    ivec indices = regspace<ivec>(0, N - 1);
#ifdef _OPENMP
#pragma omp parallel for shared(indices,treeNeighborhoods)
#endif
    for (int t = 0; t < n_trees; t++) if (! p.check_abort()) {
      searchTree(threshold,
                 indices,
                 dataMat,
                 treeNeighborhoods,
                 p
      );

      if (t > 0 && ! p.check_abort())
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        for (int i = 0; i < N; i++) {
          vector<int>* neighbors = treeNeighborhoods[i];
          sort(neighbors -> begin(), neighbors -> end());
          vector<int>::iterator theEnd = unique(neighbors -> begin(), neighbors -> end());
          neighbors -> erase(theEnd, neighbors -> end());
          if (neighbors -> size() < 3) stop("Tree failure.");
        }
      }
    }
  }

  if (p.check_abort()) return mat(0);

  // Initialize the knn matrix, and reduce the number of candidate neighbors per node
  // to K.  Otherwise the first neighborhood exploration pass takes N * trees * (threshold + 1),
  // instead of (N * K), which is prohibitive of large thresholds.
  mat knns = mat(threshold,N);
  knns.fill(-1);
#ifdef _OPENMP
#pragma omp parallel for shared(knns)
#endif
  for (int i = 0; i < N; i++) if (p.increment()){
    const SpSubview<double> x_i = data.col(i);
    priority_queue<HeapObject> MaxHeap = priority_queue<HeapObject>();
    vector<int>* stack = treeNeighborhoods[i];
    for (vector<int>::iterator it = stack -> begin(); it != stack -> end(); it++) {
      const double d = distanceFunction(x_i, data.col(*it));
      MaxHeap.push(HeapObject(d, *it));
      if (MaxHeap.size() > threshold) MaxHeap.pop();
    }
    int j = 0;
    do {
      knns(j,i) = MaxHeap.top().n;
      MaxHeap.pop();
      j++;
    } while (j < threshold && ! MaxHeap.empty());
    if (j == 1 && knns(0,i) == -1) stop("Bad neighbor matrix.");
  }
  if (p.check_abort()) return mat(0);

  for (int T = 0; T < maxIter; T++) {
    mat old_knns = knns;
    knns = mat(K,N);
    knns.fill(-1);
#ifdef _OPENMP
#pragma omp parallel for shared(knns, treeNeighborhoods)
#endif
    for (int i = 0; i < N; i++) if (p.increment()) {
      double d;

      const vec neighborhood = old_knns.col(i);
      const SpSubview<double> x_i = data.col(i);

      priority_queue<HeapObject> heap;
      vector<int> pastVisitors = *(treeNeighborhoods[i]);
      pastVisitors.reserve((K + 1) * K);
      // Loop through immediate neighbors of i
      for (int jidx = 0; jidx < old_knns.n_rows; jidx++) {
        const int j = neighborhood[jidx];
        if (j == -1) break;
        if (j == i) continue; // This should never happen
        d = distanceFunction(x_i, data.col(j));
        if (d == 0) continue; // duplicate
        heap.push(HeapObject(d, j));
        if (heap.size() > K) heap.pop();

        // For each immediate neighbor j, loop through its neighbors
        const vec locality = old_knns.col(j);
        for (int kidx = 0; kidx < old_knns.n_rows; kidx++) {
          const int k = locality[kidx];
          if (k == -1) break;
          if (k == i) continue;
          // Check if this is a neighbor we've already seen.  O(log k)
          pair<vector<int>::iterator,
               vector<int>::iterator > firstlast = equal_range(pastVisitors.begin(),
                                                               pastVisitors.end(),
                                                               k);
          if (*(firstlast.first) == k) continue; // Found

          if (firstlast.second == pastVisitors.end()) pastVisitors.push_back(k);
          else pastVisitors.insert(firstlast.second, k);

          d = distanceFunction(x_i, data.col(k));
          if (d == 0) continue;
          if (heap.size() < K) heap.push(HeapObject(d,k));
          else if (d < heap.top().d) {
            heap.push(HeapObject(d, k));
            if (heap.size() > K) heap.pop();
          }
        }
      }
      int j = 0;
      while (j < K && ! heap.empty()) {
        knns(j, i) = heap.top().n;
        heap.pop();
        j++;
      }
      if (j == 0) stop("Failure in neighborhood exploration - this should never happen.");
      vector<int>(pastVisitors).swap(pastVisitors); // pre-C++11 shrink
    }
  }
  return knns;
};


// [[Rcpp::export]]
arma::mat searchTreesCSparse(const int& threshold,
                             const int& n_trees,
                             const int& K,
                             const int& maxIter,
                             const arma::uvec& i,
                             const arma::uvec& p,
                             const arma::vec& x,
                             const std::string& distMethod,
                             bool verbose) {
  const int N = p.size() -1;
  const sp_mat data = sp_mat(i,p,x,N,N);
  return searchTreesSparse(threshold,n_trees,K,maxIter,data,distMethod,verbose);
}

// [[Rcpp::export]]
arma::mat searchTreesTSparse(const int& threshold,
                             const int& n_trees,
                             const int& K,
                             const int& maxIter,
                             const arma::uvec& i,
                             const arma::uvec& j,
                             const arma::vec& x,
                             const std::string& distMethod,
                             bool verbose) {
  const umat locations = join_cols(i,j);
  const sp_mat data = sp_mat(locations,x);
  return searchTreesSparse(threshold,n_trees,K,maxIter,data,distMethod,verbose);
}
