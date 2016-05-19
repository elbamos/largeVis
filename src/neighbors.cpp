#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp,cpp11)]]
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <set>
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

// The Euclidean distance between two vectors
inline double dist(NumericVector i, NumericVector j) {
  return sum(pow(i - j, 2));
}

inline double dist(arma::vec i, arma::vec j) {
  return sum(pow(i - j, 2));
}

void searchTree(int threshold,
                const std::vector<int>& indices,
                const NumericMatrix& data,
                std::vector<std::set<int> >& heap,
                const int& iterations,
                Function& callback,
                List& progressParams) {
  const int I = indices.size();
  if (I < 2) return;
  if (I == 2) {
      heap[indices[0]].insert(indices[1]);
      heap[indices[1]].insert(indices[0]);
      return;
  }
  if (I < threshold || iterations == 0) {
    int i = 0;
    const int mx = (threshold < I) ? threshold : I;
    do {
      int j = i + 1;
      do {
        heap[indices[i]].insert(indices[j]);
        heap[indices[j]].insert(indices[i]);
        j++;
      } while (j < I && (j % threshold) < mx);
      i++;
    } while(i < I - 1);
    callback(I, _["tokens"] = progressParams);
    return;
  }
  std::vector<int> left = std::vector<int>();
  std::vector<int> right = std::vector<int>();
  int x1idx, x2idx;
  do {
    const arma::vec selections = arma::randu(2) * (I - 1);
    x1idx = indices[selections[0]];
    x2idx = indices[selections[1]];
  } while (x1idx == x2idx);

  {
    const NumericVector x2 = data.row(x2idx);
    const NumericVector x1 = data.row(x1idx);
    // Get hyperplane
    const NumericVector m =  (x1 + x2) / 2; // Base point of hyperplane
    const NumericVector v =  (x1 - x2) / sqrt(sum(pow(x1 - x2,2))); // Unit vector

    for (int i = 0; i < indices.size(); i++) {
      const int I = indices[i];
      if (sum((data(I,_) - m) * v) > 0) {
        left.push_back(I);
      } else {
        right.push_back(I);
      };
    }
  }
  left.shrink_to_fit();
  right.shrink_to_fit();

  #pragma omp parallel sections
  {
    #pragma omp section
    {
      searchTree(threshold, left, data, heap, iterations - 1, callback, progressParams);
    }
    #pragma omp section
    {
      searchTree(threshold, right, data, heap, iterations - 1, callback, progressParams);
    }
  }
};

// [[Rcpp::export]]
arma::mat searchTrees(const int& threshold,
                      const int& n_trees,
                      const int& K,
                      const int& max_recursion_degree,
                      const int& maxIter,
                      NumericMatrix data,
                      Function callback) {
  const int N = data.nrow();
  std::vector<std::set<int> > treeNeighborhoods = std::vector<std::set<int> >(N);

  std::vector<int> indices = std::vector<int>(N);
  std::iota(indices.begin(), indices.end(), 0);

  Rcpp::List progressParams = Rcpp::List::create();

  for (int t = 0; t < n_trees; t++) {
    progressParams["phase"] = "Random Projection Tree " +  std::to_string(t + 1) + "/" + std::to_string(n_trees);
    searchTree(threshold,
               indices,
               data,
               treeNeighborhoods,
               max_recursion_degree, // maximum permitted level of recursion
               callback,
               progressParams
               );
  }
  arma::mat knns = arma::mat(threshold,N);
  knns.fill(-1);
  progressParams["phase"] = "Reducing Tree Leaves";
  #pragma omp parallel for shared(knns)
  for (int i = 0; i < N; i++) {
    const NumericVector x_i = data.row(i);
    std::priority_queue<heapObject> maxHeap = std::priority_queue<heapObject>();
    std::set<int> stack = treeNeighborhoods[i];
    for (std::set<int>::iterator it = stack.begin(); it != stack.end(); it++) {
      double d = dist(x_i, data.row(*it));
      maxHeap.push(heapObject(d, *it));
      if (maxHeap.size() > threshold) maxHeap.pop();
    }
    int j = 0;
    do {
      knns(j,i) = maxHeap.top().n;
      maxHeap.pop();
      j++;
    } while (j < threshold && ! maxHeap.empty());
    stack.insert(i);
    if (i > 0 && i % 1000 == 0) callback(1000, _["tokens"] = progressParams);
  }

  for (int T = 0; T < maxIter; T++) {
    progressParams["phase"] = "Exploring Neighborhood " + std::to_string(T + 1) + "/" +  std::to_string(maxIter);
    arma::mat old_knns = knns;
    knns = arma::mat(K,N);
    knns.fill(-1);
    #pragma omp parallel for shared(knns, treeNeighborhoods)
    for (int i = 0; i < N; i++) {
      double d;
      if (i > 0 && i % 1000 == 0) callback(1000, _["tokens"] = progressParams);

      const NumericVector x_i = data.row(i);
      std::priority_queue<heapObject> heap;
      std::set<int> pastVisitors = treeNeighborhoods[i];

      const arma::vec neighborhood = old_knns.col(i);
      int cnt = sum(neighborhood == -1);
      arma::vec locality;
      for (int jidx = 0; jidx < old_knns.n_rows; jidx++) {
        const int j = neighborhood[jidx];
        if (j == -1) break;
        if (j != i) {
          d = dist(x_i, data.row(j));
          if (d != 0) {
            heap.push(heapObject(d, j));
            if (heap.size() > K) heap.pop();
          }
        }
        locality = old_knns.col(j);
        for (int kidx = 0; kidx < old_knns.n_rows; kidx++) {
          const int k = locality[kidx];
          if (k == -1) break;
          if (k != i && pastVisitors.insert(k).second) {
            d = dist(x_i, data.row(k));
            if (d != 0 && d < heap.top().d) {
              heap.push(heapObject(d, k));
              if (heap.size() > K) heap.pop();
            }
          }
        }
      }
      int j = 0;
      while (j < K && ! heap.empty()) {
        knns(j, i) = heap.top().n;
        heap.pop();
        j++;
      }
    }
  }
  return knns;
};
