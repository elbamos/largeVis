#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]
#include "progress.hpp"
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <set>
#include "helpers.h"
using namespace Rcpp;
using namespace std;

/*
* Functions for identifying candidate nearest neighbors using random projection trees and neighborhood exploration.
*/

struct heapObject {
  double d;
  int n;

  heapObject(double d, int n) : d(d), n(n) {}

  bool operator<(const struct heapObject& other) const {
    return d < other.d;
  }
};

double relDist(const arma::vec& i, const arma::vec& j) {
  return sum(square(i - j));
}

void searchTree(const int& threshold,
                const arma::vec& indices,
                const arma::mat& data,
                std::vector<std::vector<int>* >& heap,
                const int& iterations,
                Progress& progress) {
  const int I = indices.size();
  const int D = data.n_rows;
  if (progress.check_abort()) return;
  if (I < 2) return;
  if (I == 2) {
    #pragma omp critical
    {
      heap[indices[0]] -> push_back(indices[1]);
      heap[indices[1]] -> push_back(indices[0]);
    }
    return;
  }
  if (I < threshold || iterations == 0) {
    #pragma omp critical
    {
      for (int i = 0; i < I; i++) {
        heap[indices[i]] -> reserve(I);
        for (int j = 0; j < I; j++)
          if (i != j) heap[indices[i]] -> push_back(indices[j]);
      }
    }
    progress.increment(I);
    return;
  }
  arma::vec direction = arma::vec(indices.size());
  int x1idx, x2idx;
  do {
    const arma::vec selections = arma::randu(2) * (I - 1);
    x1idx = indices[selections[0]];
    x2idx = indices[selections[1]];
  } while (x1idx == x2idx);

  {
    const arma::vec x2 = data.col(x2idx);
    const arma::vec x1 = data.col(x1idx);
    // Get hyperplane
    const arma::vec m =  (x1 + x2) / 2; // Base point of hyperplane
    const arma::vec v =  ((x1 - x2) / sqrt(sum(pow(x1 - x2,2)))); // Unit vector

    for (int i = 0; i < indices.size(); i++) {
      const int I = indices[i];
      const arma::vec X = data.col(I);
      direction[i] = sum((X - m).t() * v);
    }
    // Normalize direction
    const double middle = arma::median(direction);
    direction = direction - middle;
  }

  searchTree(threshold, indices(arma::find(direction > 0)), data, heap, iterations - 1, progress);
  searchTree(threshold, indices(arma::find(direction <= 0)), data, heap, iterations - 1, progress);
};



// [[Rcpp::export]]
arma::mat searchTrees(const int& threshold,
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

  Progress p((N * n_trees) + (N) + (N * maxIter), verbose);

  std::vector<std::vector<int>* > treeNeighborhoods = std::vector<std::vector<int>* >(N);
  for (int i = 0; i < N; i++) {
    int seed[] = {i};
    treeNeighborhoods[i] = new std::vector<int>(seed, seed + sizeof(seed) / sizeof(int));
  }

  { // Artificial scope to destroy indices
    arma::vec indices = arma::regspace<arma::vec>(0, N - 1);

    #pragma omp parallel for shared(indices,treeNeighborhoods)
    for (int t = 0; t < n_trees; t++) if (! p.check_abort()) {
      searchTree(threshold,
                 indices,
                 data,
                 treeNeighborhoods,
                 max_recursion_degree, // maximum permitted level of recursion
                 p
      );

      if (t > 0 && ! p.check_abort())
      #pragma omp critical
      {
        for (int i = 0; i < N; i++) {
          std::vector<int>* neighbors = treeNeighborhoods[i];
          std::sort(neighbors -> begin(), neighbors -> end());
          std::vector<int>::iterator theEnd = std::unique(neighbors -> begin(), neighbors -> end());
          neighbors -> erase(theEnd, neighbors -> end());
        }
      }
    }
  }

  if (p.check_abort()) return arma::mat(0);

  // Initialize the knn matrix, and reduce the number of candidate neighbors per node
  // to K.  Otherwise the first neighborhood exploration pass takes N * trees * (threshold + 1),
  // instead of (N * K), which is prohibitive of large thresholds.
  arma::mat knns = arma::mat(threshold,N);
  knns.fill(-1);
  #pragma omp parallel for shared(knns)
  for (int i = 0; i < N; i++) if (p.increment()){
    const arma::vec x_i = data.col(i);
    std::priority_queue<heapObject> maxHeap = std::priority_queue<heapObject>();
    std::vector<int>* stack = treeNeighborhoods[i];
    for (std::vector<int>::iterator it = stack -> begin(); it != stack -> end(); it++) {
      const double d = distanceFunction(x_i, data.col(*it));
      maxHeap.push(heapObject(d, *it));
      if (maxHeap.size() > threshold) maxHeap.pop();
    }
    int j = 0;
    do {
      knns(j,i) = maxHeap.top().n;
      maxHeap.pop();
      j++;
    } while (j < threshold && ! maxHeap.empty());
  }
  if (p.check_abort()) return arma::mat(0);

  for (int T = 0; T < maxIter; T++) {
    arma::mat old_knns = knns;
    knns = arma::mat(K,N);
    knns.fill(-1);
    #pragma omp parallel for shared(knns, treeNeighborhoods)
    for (int i = 0; i < N; i++) if (p.increment()) {
      double d;

      const arma::vec neighborhood = old_knns.col(i);
      const arma::vec x_i = data.col(i);

      std::priority_queue<heapObject> heap;
      std::vector<int> pastVisitors = *(treeNeighborhoods[i]);
      // Loop through immediate neighbors of i
      for (int jidx = 0; jidx < old_knns.n_rows; jidx++) {
        const int j = neighborhood[jidx];
        if (j == -1) break;
        if (j == i) continue; // This should never happen
        d = distanceFunction(x_i, data.col(j));
        if (d != 0) {
          heap.push(heapObject(d, j));
          if (heap.size() > K) heap.pop();
        }

        pastVisitors.reserve(K);
        // For each immediate neighbor j, loop through its neighbors
        const arma::vec locality = old_knns.col(j);
        for (int kidx = 0; kidx < old_knns.n_rows; kidx++) {
          const int k = locality[kidx];
          if (k == -1) break;
          if (k == i) continue;
          // Check if this is a neighbor we've already seen.  O(log k)
          std::pair<std::vector<int>::iterator,
                    std::vector<int>::iterator > firstlast = std::equal_range(pastVisitors.begin(),
                                                                              pastVisitors.end(),
                                                                              k);
          if (*(firstlast.first) == k) continue; // Found

          if (firstlast.second == pastVisitors.end()) pastVisitors.push_back(k);
          else pastVisitors.insert(firstlast.second, k);

          d = distanceFunction(x_i, data.col(k));
          if (d != 0 && d < heap.top().d) {
            heap.push(heapObject(d, k));
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
      std::vector<int>(pastVisitors).swap(pastVisitors); // pre-C++11 shrink
    }
  }
  return knns;
};
