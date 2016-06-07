#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp,cpp11)]]
#include "progress.hpp"
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <set>
#include <memory>
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

typedef std::priority_queue<heapObject> maxHeap;
typedef std::vector< arma::imat::col_iterator > positionVector;
typedef std::vector<int> neighborhood;

void searchTree(const int& threshold,
                const arma::ivec& indices,
                const arma::mat& data,
                neighborhood* heap[],
                const int& iterations,
                Progress& progress) {
  const int I = indices.size();
  const int D = data.n_rows;
  if (progress.check_abort()) return;
  if (I < 2) stop("Tree split failure.");
  if (I <= threshold || iterations == 0) {
    arma::ivec neighbors = arma::ivec(indices);
    neighborhood tmpStorage = neighborhood();
    arma::ivec::iterator newEnd = neighbors.end();
    #pragma omp critical
    {
        for (arma::ivec::iterator it = neighbors.begin();
             it != newEnd;
             it++) {
          tmpStorage.clear();
          tmpStorage.swap(*heap[*it]);
          heap[*it] -> reserve(tmpStorage.size() + I);
          arma::ivec::iterator newIt = neighbors.begin();
          std::vector<int>::iterator oldIt = tmpStorage.begin();
          std::vector<int>::iterator oldEnd = tmpStorage.end();
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
    progress.increment(I);
    return;
  }
  arma::vec direction = arma::vec(indices.size());
  {
    int x1idx, x2idx;
    arma::vec v;
    arma::vec m;
    do {
      const arma::vec selections = arma::randu(2) * (I - 1);
      x1idx = indices[selections[0]];
      x2idx = indices[selections[1]];
      if (x1idx == x2idx) x2idx = indices[((int)selections[1] + 1) % indices.size()];
      const arma::vec x2 = data.col(x2idx);
      const arma::vec x1 = data.col(x1idx);
      // Get hyperplane
      m =  (x1 + x2) / 2; // Base point of hyperplane
      const arma::vec d = x1 - x2;
      v =  d / arma::as_scalar(arma::norm(d, 2)); // unit vector
    } while (x1idx == x2idx);

    for (int i = 0; i < indices.size(); i++) {
      const arma::vec X = data.col(indices[i]);
      direction[i] = dot((X - m), v);
    }
  }
  // Normalize direction
  const double middle = arma::median(direction);

  const arma::uvec left = arma::find(direction > middle);
  const arma::uvec right = arma::find(direction <= middle);
  if (left.size() >= 2 && right.size() >= 2) {
    searchTree(threshold, indices(left), data, heap, iterations - 1, progress);
    searchTree(threshold, indices(right), data, heap, iterations - 1, progress);
  } else { // Handles the rare case where the split fails because of equidistant points
    searchTree(threshold, indices.subvec(0, indices.size() / 2), data, heap, iterations - 1, progress);
    searchTree(threshold, indices.subvec(indices.size() / 2, indices.size() - 1), data, heap, iterations - 1, progress);
  }
};

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

  std::set<int>** treeHolder = new std::set<int>*[N];
  { // Artificial scope to destroy indices and treeneighborhoods
    // array of vectors to hold lists of neighbors for each n from all trees
    neighborhood** treeNeighborhoods = new neighborhood*[N];
    for (int i = 0; i < N; i++) {
      int seed[] = {i};
      treeNeighborhoods[i] = new std::vector<int>(seed, seed + sizeof(seed) / sizeof(int));
    }
    const arma::ivec indices = arma::regspace<arma::ivec>(0, N - 1);

    #pragma omp parallel for shared(treeNeighborhoods)
    for (int t = 0; t < n_trees; t++) if (! p.check_abort()) {
      searchTree(threshold,
                 indices,
                 data,
                 treeNeighborhoods,
                 max_recursion_degree, // maximum permitted level of recursion
                 p
      );
    }
    if (p.check_abort()) return arma::imat(0);
    // Reduce size from threshold * n_trees to top K, and sort
    treeHolder = new std::set<int>*[N];
    maxHeap thisHeap = maxHeap();
    #pragma omp parallel for shared(treeHolder, treeNeighborhoods) private(thisHeap)
    for (int i = 0; i < N; i++) if (p.increment()) {
      const arma::vec x_i = data.col(i);
      std::vector<int> *neighborhood = treeNeighborhoods[i];
      for (std::vector<int>::iterator j = neighborhood -> begin();
           j != neighborhood -> end();
           j++) {

        const double d = distanceFunction(x_i, data.col(*j));
        if (d != 0) {
          thisHeap.emplace(d, *j);
          if (thisHeap.size() > K) thisHeap.pop();
        }
      }
      delete treeNeighborhoods[i];
      treeHolder[i] = new std::set<int>();
      while (! thisHeap.empty()) {
        treeHolder[i] -> emplace(thisHeap.top().n);
        thisHeap.pop();
      }
    }
    delete treeNeighborhoods;
  }
  // Copy sorted neighborhoods into matrix. This is faster than
  // sorting in-place.
  arma::imat knns = arma::imat(K,N);
  #pragma omp parallel for
  for (int i = 0; i < N; i++) if (p.increment()) {
    std::set<int>::iterator sortIterator = treeHolder[i] -> begin();
    std::set<int>::iterator end = treeHolder[i] -> end();
    int j = 0;
    while (sortIterator != end) knns(j++, i) = *sortIterator++;
    if (j == 0) stop("Tree failure.");
    while (j < K) knns(j++,i) = -1;
    delete treeHolder[i];
  }

  if (p.check_abort()) return arma::imat(0);
  arma::imat old_knns  = arma::imat(K,N);
  for (int T = 0; T < maxIter; T++) if (! p.check_abort()) {
    arma::imat tmp = old_knns;
    old_knns = knns;
    knns = tmp;
    maxHeap thisHeap = maxHeap();
    std::set<int> sorter = std::set<int>();
    #pragma omp parallel for shared(old_knns, knns) private(thisHeap, sorter)
    for (int i = 0; i < N; i++) if (p.increment()) {
      double d;

      const arma::vec x_i = data.col(i);

      positionVector positions = positionVector();
      positionVector ends = positionVector();

      positions.reserve(K + 1);
      ends.reserve(K + 1);

      positions.push_back(old_knns.begin_col(i));
      ends.push_back(old_knns.end_col(i));

      for (arma::imat::col_iterator it = old_knns.begin_col(i);
           it != ends[0] && *it != -1;
           it++) {
        positions.push_back(old_knns.begin_col(*it));
        ends.push_back(old_knns.end_col(*it));
      }

      int lastOne = N + 1;
      // This is a K + 1 vector merge sort running in O(K * N)
      positionVector::iterator theEnd = positions.end();
      while (true) {
        arma::imat::col_iterator whch = 0;

        for (std::pair< positionVector::iterator,
                        positionVector::iterator > it(positions.begin(),
                                                      ends.begin());
             it.first != theEnd;
             it.first++, it.second++) while (*it.first != *it.second) { // For each neighborhood, keep going until
                                                                        // we find a non-dupe or get to the end

          if (**it.first == -1) std::advance(*it.first, std::distance(*it.first, *it.second));
          else if (**it.first == i || **it.first == lastOne) std::advance(*it.first, 1);
          else if (whch == 0 || **it.first < *whch) {
            whch = *it.first;
            break;
          } else break;
        }

        if (whch == 0) break;
        lastOne = *whch;
        std::advance(whch, 1);

        d = distanceFunction(x_i, data.col(lastOne));
        if (d > 0) {
          thisHeap.emplace(d, lastOne);
          if (thisHeap.size() > K) thisHeap.pop();
        }
      }

      sorter.clear();
      while (!thisHeap.empty()) {
        sorter.emplace(thisHeap.top().n);
        thisHeap.pop();
      }
      std::set<int>::iterator sortIterator = sorter.begin();
      int j = 0;
      while (sortIterator != sorter.end()) knns(j++, i) = *sortIterator++;
      if (j == 0) stop("Neighbor exploration failure.");
      while (j < K) knns(j++,i) = -1;
    }
  }
  return knns;
};
