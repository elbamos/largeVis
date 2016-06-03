#include <RcppArmadillo.h>
#include "progress.hpp"
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends("RcppArmadillo")]]
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <string>
#include "helpers.h"

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
  if (I < 2) stop("Tree split failure.");
  if (I == 2) {
  //  #pragma omp critical
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
        std::vector<int>* thisHeap = heap[indices[i]];
        thisHeap -> reserve(I);
        thisHeap -> insert(thisHeap -> end(), indices.begin(), indices.end());
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
      const int I = indices[i];
      const arma::vec X = data.col(I);
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
  } else {
    searchTree(threshold, indices.subvec(0, indices.size() / 2), data, heap, iterations - 1, progress);
    searchTree(threshold, indices.subvec(indices.size() / 2, indices.size() - 1), data, heap, iterations - 1, progress);
  }
};

template <class T, class S, class C>
S& Container(std::priority_queue<T, S, C>& q) {
  struct HackedQueue : private std::priority_queue<T, S, C> {
    static S& Container(std::priority_queue<T, S, C>& q) {
      return q.*&HackedQueue::c;
    }
  };
  return HackedQueue::Container(q);
}

template <class T>
T front(std::set<T>& q) {
  return * q.begin();
}

typedef std::priority_queue<heapObject, std::set<heapObject> > maxHeap;
typedef std::vector< maxHeap* > heapVector;

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
  heapVector setOfHeaps = heapVector(N);

  { // Artificial scope to destroy indices
    std::vector<std::vector<int>* > treeNeighborhoods = std::vector<std::vector<int>* >(N);
    for (int i = 0; i < N; i++) {
      int seed[] = {i};
      treeNeighborhoods[i] = new std::vector<int>(seed, seed + sizeof(seed) / sizeof(int));
    }
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
    }
    if (p.check_abort()) return arma::mat(0);

    #pragma omp parallel for shared(treeNeighborhoods)
    for (int i = 0; i < N; i++) if (p.increment()){
      std::vector<int>* neighbors = treeNeighborhoods[i];
      const arma::vec x_i = data.col(i);
      maxHeap *thisHeap = new maxHeap();

      for (std::vector<int>::iterator it = neighbors -> begin(); it != neighbors -> end(); it++) {
        if (*it == i) continue;
        const double d = distanceFunction(x_i, data.col(*it));
        if (d == 0) continue;
        if (d < thisHeap -> top().d) {
          thisHeap -> push(heapObject(d, *it));
          if (thisHeap-> size() > K) thisHeap -> pop();
        }
      }
      setOfHeaps[i] = thisHeap;
    }
  }

  // Initialize the knn matrix, and reduce the number of candidate neighbors per node
  // to K.  Otherwise the first neighborhood exploration pass takes N * trees * (threshold + 1),
  // instead of (N * K), which is prohibitive of large thresholds.

  if (p.check_abort()) return arma::mat(0);

  for (int T = 0; T < maxIter; T++) {
#pragma omp parallel for shared(setOfHeaps)
    for (int i = 0; i < N; i++) if (p.increment()) {
      double d;

      std::set<heapObject> &oldNeighborhood = Container(*setOfHeaps[i]);
      std::set<heapObject>::iterator oldIt = oldNeighborhood.end();
      maxHeap *newiHeap = new maxHeap(oldNeighborhood.begin(),
                                      oldIt);
      const arma::vec x_i = data.col(i);

      for (std::set<heapObject>::iterator j = oldNeighborhood.begin();
           j != oldIt;
           j++) {
        int jidx = j -> n;

        std::set<heapObject> &nextNeighbors = Container(*setOfHeaps[jidx]);
        for (std::set<heapObject>::iterator k = nextNeighbors.begin();
             k != nextNeighbors.end();
             k++) {
          int kidx = k -> n;
          if (kidx == i) continue;
          d = distanceFunction(x_i, data.col(jidx));
          if (d < newiHeap -> top().d) {
            newiHeap -> push(heapObject(d, kidx));
            if (newiHeap -> size() > K) newiHeap -> pop();
          }
        }
      }
      setOfHeaps[i] = newiHeap;
    }
  }
  arma::mat knns = arma::mat(K, N);
#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    int j = 0;
    maxHeap *thisHeap = setOfHeaps[i];
    do {
      knns(j, i) = thisHeap -> top().n;
      j++;
    } while (! thisHeap -> empty());
  }

  return knns;
};
