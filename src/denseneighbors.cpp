// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

/*
* When a leaf node is found, store the identity of the neighbors
* along with each vertex.
*/
void addNeighbors(const arma::ivec& indices,
                  Neighborhood* heap[],
                                    const int I) {
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
    tmpStorage.swap(*heap[*it]);
    heap[*it] -> reserve(tmpStorage.size() + I);
    ivec::iterator newIt = neighbors.begin();
    vector<int>::iterator oldIt = tmpStorage.begin();
    vector<int>::iterator oldEnd = tmpStorage.end();
    int last;
    int best = -1;
    while (oldIt != oldEnd || newIt != newEnd) {
      if (oldIt == oldEnd) best = *newIt++;
      else if (newIt == newEnd) best = *oldIt++;
      else best = (*newIt < *oldIt) ? *newIt++ : *oldIt++;
      if (best == last || best == *it) continue;
      heap[*it] -> push_back(best);
      last = best;
    }
  }
}
}

arma::vec hyperplane(const arma::ivec& indices,
                     const arma::mat& data,
                     const int I) {
  vec direction = vec(indices.size());
  int x1idx, x2idx;
  vec v;
  vec m;
  do {
    const vec selections = randu(2) * (I - 1);
    x1idx = indices[selections[0]];
    x2idx = indices[selections[1]];
    if (x1idx == x2idx) x2idx = indices[((int)selections[1] + 1) % indices.size()];
    const vec x2 = data.col(x2idx);
    const vec x1 = data.col(x1idx);
    // Get hyperplane
    m =  (x1 + x2) / 2; // Base point of hyperplane
    const vec d = x1 - x2;
    v =  d / as_scalar(norm(d, 2)); // unit vector
  } while (x1idx == x2idx);

  for (int i = 0; i < indices.size(); i++) {
    const vec X = data.col(indices[i]);
    direction[i] = dot((X - m), v);
  }
  return direction;
}

/*
* The recursive function for the annoy neighbor search
* algorithm. Partitions space by a random hyperplane,
* and calls itself recursively (twice) on each side.
*
* If called with fewer nodes than the threshold,
*
*/
void searchTree(const int& threshold,
                const arma::ivec& indices,
                const arma::mat& data,
                Neighborhood* heap[],
                                  const int& iterations,
                                  Progress& progress) {
  const int I = indices.size();
  // const int D = data.n_rows;
  if (progress.check_abort()) return;
  if (I < 2) stop("Tree split failure.");
  if (I <= threshold || iterations == 0) {
    addNeighbors(indices, heap, I);
    progress.increment(I);
    return;
  }
  vec direction = hyperplane(indices, data, I);
  const double middle = median(direction);
  const uvec left = find(direction > middle);
  const uvec right = find(direction <= middle);

  if (left.size() >= 2 && right.size() >= 2) {
    searchTree(threshold, indices(left), data, heap, iterations - 1, progress);
    searchTree(threshold, indices(right), data, heap, iterations - 1, progress);
  } else { // Handles the rare case where the split fails because of equidistant points
    searchTree(threshold, indices.subvec(0, indices.size() / 2), data, heap, iterations - 1, progress);
    searchTree(threshold, indices.subvec(indices.size() / 2, indices.size() - 1), data, heap, iterations - 1, progress);
  }
};

Neighborhood** createNeighborhood(int N) {
  Neighborhood** treeNeighborhoods = new Neighborhood*[N];
  for (int i = 0; i < N; i++) {
    int seed[] = {i};
    treeNeighborhoods[i] = new vector<int>(seed, seed + sizeof(seed) / sizeof(int));
  }
  return treeNeighborhoods;
}

void copyHeapToMatrix(set<int>* tree,
                      const int K,
                      const int i,
                      arma::imat& knns) {
  set<int>::iterator sortIterator = tree -> begin();
  set<int>::iterator end = tree -> end();
  int j = 0;
  while (sortIterator != end) knns(j++, i) = *sortIterator++;
  if (j == 0) stop("Tree failure.");
  while (j < K) knns(j++, i) = -1;
}

void addDistance(const arma::vec& x_i,
                 const arma::mat& data,
                 const int j,
                 MaxHeap& heap,
                 const int K,
                 double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j)) {
  const double d = distanceFunction(x_i, data.col(j));
  if (d != 0) {
    heap.emplace(d, j);
    if (heap.size() > K) heap.pop();
  }
}

void heapToSet(MaxHeap& heap, set<int>* set) {
  while (! heap.empty()) {
    set -> emplace(heap.top().n);
    heap.pop();
  }
}

arma::imat annoy(const int n_trees,
                 const int threshold,
                 const arma::mat& data,
                 const int max_recursion_degree,
                 const int N,
                 const int K,
                 double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j),
                 Progress& p) {
  set<int>** treeHolder = new set<int>*[N];
  Neighborhood** treeNeighborhoods = createNeighborhood(N);
  const ivec indices = regspace<ivec>(0, N - 1);

#ifdef _OPENMP
#pragma omp parallel for shared(treeNeighborhoods)
#endif
  for (int t = 0; t < n_trees; t++) if (! p.check_abort())
    searchTree(threshold,
               indices,
               data,
               treeNeighborhoods,
               max_recursion_degree, // maximum permitted level of recursion
               p
    );
  if (p.check_abort()) return imat(K, N);
  // Reduce size from threshold * n_trees to top K, and sort
  MaxHeap thisHeap = MaxHeap();
#ifdef _OPENMP
#pragma omp parallel for shared(treeHolder, treeNeighborhoods) private(thisHeap)
#endif
  for (int i = 0; i < N; i++) if (p.increment()) {
    const vec x_i = data.col(i);
    vector<int> *neighborhood = treeNeighborhoods[i];
    for (vector<int>::iterator j = neighborhood -> begin();
         j != neighborhood -> end();
         j++)
      addDistance(x_i, data, *j, thisHeap, K, distanceFunction);
    delete treeNeighborhoods[i];
    treeHolder[i] = new set<int>();
    heapToSet(thisHeap, treeHolder[i]);
  }
  // Copy sorted neighborhoods into matrix. This is faster than
  // sorting in-place.
  imat knns = imat(K,N);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < N; i++) if (p.increment()) {
    copyHeapToMatrix(treeHolder[i], K, i, knns);
    delete treeHolder[i];
  }
  return knns;
}

// [[Rcpp::export]]
arma::imat searchTrees(const int& threshold,
                       const int& n_trees,
                       const int& K,
                       const int& max_recursion_degree,
                       const int& maxIter,
                       const arma::mat& data,
                       const std::string& distMethod,
                       bool verbose) {

  const int N = data.n_cols;

  double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = relDist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = cosDist;
  else distanceFunction = relDist;

  Progress p((N * n_trees) + (2 * N) + (N * maxIter), verbose);

  imat knns;
  {
  	mat dataMat;
  	if (distMethod.compare(std::string("Cosine")) == 0) dataMat = normalise(data);
  	else dataMat = data;
  	knns = annoy(n_trees,
                    threshold,
                    dataMat,
                    max_recursion_degree,
                    N,
                    K,
                    distanceFunction,
                    p);
  }

  if (p.check_abort()) return imat(0);
  imat old_knns  = imat(K,N);
  for (int T = 0; T < maxIter; T++) if (! p.check_abort()) {
    imat tmp = old_knns;
    old_knns = knns;
    knns = tmp;
    MaxHeap thisHeap = MaxHeap();
    set<int> sorter = set<int>();
#ifdef _OPENMP
#pragma omp parallel for shared(old_knns, knns) private(thisHeap, sorter)
#endif
    for (int i = 0; i < N; i++) if (p.increment()) {
      const vec x_i = data.col(i);

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

      int lastOne = N + 1;
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

             addDistance(x_i, data, lastOne, thisHeap, K, distanceFunction);
      }

      sorter.clear();
      heapToSet(thisHeap, &sorter);

      set<int>::iterator sortIterator = sorter.begin();
      int j = 0;
      while (sortIterator != sorter.end()) knns(j++, i) = *sortIterator++;
      if (j == 0) stop("Neighbor exploration failure.");
      while (j < K) knns(j++,i) = -1;
    }
  }
  return knns;
};

