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

/*
 * Neighborhood exploration.
 */
// [[Rcpp::export]]
void neighbors_inner( int maxIter,
                      NumericMatrix old_knns,
                      NumericMatrix data,
                      NumericMatrix outputKnns,
                      bool prefilter,
                      Function callback) {
  int N = old_knns.ncol();
  int K = outputKnns.nrow();

  arma::mat knns = as<arma::mat>(old_knns);

  // Pre-filter to reduce the number of comparisons
  if (prefilter) {
    arma::mat old_knns = knns;
    knns = arma::mat(K,N);
    knns.fill(-1);

    #pragma omp parallel for shared(knns)
    for (int i = 0; i < N; i++) {
      if (i > 0 && i % 1000 == 0) callback(10);
      NumericVector x_i = data.row(i);

      std::set<int> seen;
      std::priority_queue<heapObject> heap;

      seen.insert(i);
      arma::vec neighbors = old_knns.col(i);
      for (int jidx = 0; jidx < old_knns.n_rows; jidx++) {
        const int j = neighbors[jidx];
        if (j == -1) break;
        if (j != i && seen.insert(j).second) {
          const double d = dist(x_i, data.row(j));
          if (d != 0) {
            heap.push(heapObject(d, j));
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
    }
  }

  for (int T = 0; T < maxIter; T++) {
    arma::mat old_knns = knns;
    knns = arma::mat(K,N);
    knns.fill(-1);
    #pragma omp parallel for shared(knns)
    for (int i = 0; i < N; i++) {
      double d;
      if (i > 0 && i % 1000 == 0) callback(1000);
      NumericVector x_i = data.row(i);

      std::set<int> seen;
      std::priority_queue<heapObject> heap;

      seen.insert(i);
      const arma::vec neighborhood = old_knns.col(i);
      for (int jidx = 0; jidx < old_knns.n_rows; jidx++) {
        const int j = neighborhood[jidx];
        if (j == -1) break;
        if (j != i && seen.insert(j).second) {
          d = dist(x_i, data.row(j));
          if (d != 0) {
            heap.push(heapObject(d, j));
            if (heap.size() > K) heap.pop();
          }
        }
        const arma::vec locality = old_knns.col(j);
        for (int kidx = 0; kidx < old_knns.n_rows; kidx++) {
          const int k = locality[kidx];
          if (k == -1) break;
          if (k != i && seen.insert(k).second) {
            d = dist(x_i, data.row(k));
            if (d != 0) {
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
  for (int i = 0; i < N; i++) for (int j = 0; j < K; j++) outputKnns(j,i) = knns(j,i);
};

void searchTree(int threshold,
                const arma::uvec& indices,
                const arma::mat& data,
                arma::imat& output,
                Function callback) {
  if (indices.size() < 2) return;
#pragma omp critical
{
  if (indices.size() == 2) {
    output(0, indices[0]) = indices[1];
    output(0, indices[1]) = indices[0];
    return;
  }
  if (indices.size() < threshold) {
    int i = 0;
    do {
      int j = i + 1;
      do {
        output(j, indices[i]) = indices[j];
        output(i, indices[j]) = indices[i];
        j++;
      } while (j < indices.size());
      i++;
    } while(i < indices.size() - 1);
    callback(indices.size());
    return;
  }
}
  // Get hyperplane
  arma::uvec selections = indices.elem(arma::randi<arma::uvec>(2, arma::distr_param(0, indices.size() - 1)));
  arma::vec v =  data.col(selections[1]) - data.col(selections[0]);
  arma::vec m = sum(data.cols(selections), 1) / 2;
  double mv = dot(m,v); // This is the hyperplane
  arma::vec direction = arma::vec(indices.size());
  for (int i = 0; i < indices.size(); i++) direction[i] = sum(data.col(indices[i]) % v) - mv;

  #pragma omp parallel sections
  {
        #pragma omp section
    {
      searchTree(threshold, indices.elem(arma::find(direction > 0)), data, output, callback);
    }
        #pragma omp section
    {
      searchTree(threshold, indices.elem(arma::find(direction <= 0)), data, output, callback);
    }
  }
};

// [[Rcpp::export]]
arma::mat searchTrees(int threshold,
                      int n_trees,
                      NumericMatrix data, Function callback) {
  arma::mat inputData = as<arma::mat>(data).t();
  std::vector<std::set<int> > heap = std::vector<std::set<int> >(inputData.n_cols);
#pragma omp parallel for shared(heap)
  for (int t = 0; t < n_trees; t++) {
    arma::imat output = arma::imat(threshold + 1, inputData.n_cols);
    output.fill(-1);
    searchTree(threshold,
               arma::regspace<arma::uvec>(0, inputData.n_cols - 1),
               inputData,
               output,
               callback);
    for (int i = 0; i < inputData.n_cols; i++)
      for (int j = 0; j <= threshold; j++)
        if (output(j,i) != -1) heap[i].insert(output(j,i));
  }

  int sz = 0;
  for (int i = 0; i < heap.size(); i++) if (heap[i].size() > sz) sz = heap[i].size();
  arma::mat output = arma::mat(sz, inputData.n_cols).fill(-1);

  for (int i = 0; i < heap.size(); i++) {
    std::set<int> neighbors = heap[i];
    int j = 0;
    for (std::set<int>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
      output(j, i) = *it;
      j++;
    }
  }
  return output;
}
