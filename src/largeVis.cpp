#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <set>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(openmp)]]

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

// [[Rcpp::export]]
void neighbors_inner( int maxIter,
                      NumericMatrix old_knns,
                      NumericMatrix data,
                      NumericMatrix outputKnns,
                      Function callback) {
  int N = old_knns.ncol();
  int K = outputKnns.nrow();
  int oldK;

  NumericMatrix nextKnns;
  for (int T = 0; T < maxIter; T++) {
    if (T > 0) old_knns = nextKnns;
    oldK = old_knns.nrow();

    nextKnns = NumericMatrix(K, N);
    for (int i = 0; i < K; i++) for (int j = 0; j < N; j++) nextKnns(i,j) = 0; // Initialize matrix

    for (int i = 0; i < N; i++) {
      int j, k;
      double d;
      if (i > 0 && i % 1000 == 0) callback(1000);
      NumericVector x_i = data.row(i);

      std::set<int> seen;
      std::priority_queue<heapObject> heap;

      seen.insert(i);

      for (int jidx = 0; jidx < oldK; jidx++) {
        j = old_knns.column(i)[jidx];
        if (j == 0) break;
        if (j - 1 != i && seen.insert(j - 1).second) {
          d = dist(x_i, data.row(j-1));
          heap.push(heapObject(d, j));
          if (heap.size() > K) heap.pop();
        }

        for (int kidx = 0; kidx < oldK; kidx++) {
          k = old_knns.column(j - 1)[kidx];
          if (k == 0) break;
          if (k - 1 != i && seen.insert(k - 1).second) {
            d = dist(x_i, data.row(k-1));
            heap.push(heapObject(d, k));
            if (heap.size() > K) heap.pop();
          }
        }
      }
      j = 0;
      while (j < K && ! heap.empty()) {
        nextKnns(j, i) = heap.top().n;
        heap.pop();
        j++;
      }
    }
  }
  for (int i = 0; i < N; i++) for (int j = 0; j < K; j++) outputKnns(j,i) = nextKnns(j,i);
};

/*
 * Note changes from the paper:
 * 1. \rho is set to decline non-linearly
 * 2. The gradient applied to update y_i for the positive sample is doubled (2 * \rho)
 *
 *
 */

void checkVector(NumericVector x, std::string label) {
  if (sum(is_na(x)) + sum(is_nan(x)) + sum(is_infinite(x)) > 0)
    Rcout << "\n Failure at " << label;
}

// [[Rcpp::export]]
void sgd(NumericMatrix coords,
         NumericVector positiveEdges,
         NumericVector is, // vary randomly
         NumericVector js, // ordered
         NumericVector ps, // N+1 length vector of indices to start of each row j in vector is
         NumericVector ws, // w{ij}
         int gamma,
         int rho,
         int minRho,
         bool useWeights,
         int M,
         int alpha,
         Function callback) {

  int N = ps.length() - 1;
  // Calculate negative sample weights, d_{i}^0.75.
  // Stored as a vector of cumulative sums, normalized, so it can
  // be readily searched using binary searches.
  // TODO:  CREATE A FORM FOR WHEN USEWEIGHTS = TRUE
  arma::vec negativeSampleWeights = pow(diff(ps), 0.75);
  double scale = sum(negativeSampleWeights);
  negativeSampleWeights = negativeSampleWeights / scale;
  double acc = 0;
  for (int idx = 0; idx < N; idx++) {
    acc+= negativeSampleWeights[idx];
    negativeSampleWeights[idx] = acc;
  }

  // Iterate through the edges in the positiveEdges vector
  #pragma omp parallel for shared(coords, rho)
  for (int eIdx=0; eIdx < positiveEdges.length(); eIdx++) {
    const int e_ij = positiveEdges[eIdx];
    const int i = is[e_ij];
    const int j = js[e_ij];

    NumericVector y_i = coords.row(i);
    NumericVector y_j = coords.row(j);

    // wij
    const double w = 1;// (useWeights) ? ws[e_ij] : 1;

    NumericVector yigrad = - alpha * (y_i-y_j)/(1 + (alpha * dist(y_i, y_j)));
    checkVector(yigrad, "yigrad");

    arma::vec samples = arma::randu<arma::vec>(M * 2);
    arma::vec::iterator targetIt = samples.begin();
    int sampleIdx = 1;
    int m = 0;
   // Rcout << "\n" << i << " " << j << "\t";
    while (m < M) {
      if (sampleIdx % (M * 2) == 0) samples.randu();
      // binary search for lowest number greater than the sampled number
      const double target = targetIt[sampleIdx++ % (M * 2)];
      arma::vec::iterator loc = std::upper_bound(negativeSampleWeights.begin(),
                                                 negativeSampleWeights.end(), target);
      const int k = std::distance(negativeSampleWeights.begin(), loc);
      if (k == i || k == j) continue;
      bool tst = false;
      for (int tstidx = ps(i); tstidx < ps(i + 1); tstidx++) if (is[tstidx] == k) {
        tst = true;
        break;
      }
      if (tst) continue;

      // Calculate gradient for a single negative sample
      NumericVector y_k = coords.row(k);
      const double dk = dist(y_i, y_k);
      const double alphadk1 = (alpha * dk) + 1;
      yigrad = yigrad + (
        2 * gamma * alpha * (y_i - y_k) / (
            pow((alpha * dk) + 1, 2) *
               (1 - (1 / ((alpha * dk) + 1)))
        ));
      checkVector(yigrad, "yigradneg");
      NumericVector kgrad = -2 * gamma * alpha * (y_i - y_k) /
        ((1 - (1/alphadk1))*pow(alphadk1,2));
      checkVector(kgrad, "kgrad");
      // coords(k,_) = y_k + (kgrad * rho * w);
      m++;
    }
    NumericVector yjgrad = 2 * alpha * (y_i - y_j) / (1 + (alpha * dist(y_i,y_j)));
    checkVector(yjgrad, "jgrad");

    coords(i,_) = y_i + (yigrad * w * rho);
    coords(j,_) = y_j + (yjgrad * w * rho);

    #pragma omp atomic
    rho = rho - ((rho - minRho) / (positiveEdges.length() + 1));
    if (eIdx > 0 && eIdx % 10000 == 0) callback(10000);
  }
};

// Take a matrix of data and two vectors of row indices, compute the pairwise Euclidean distance,
// and store the results in a third vector.
// [[Rcpp::export]]
void distance(NumericVector is, NumericVector js, NumericVector xs, NumericMatrix data) {
  for (int i=0; i < is.length(); i++) {
    xs[i] = sqrt(sum(pow(data.row(is[i] - 1) - data.row(js[i] - 1), 2)));
  }
};

// Take four vectors (i indices, j indices, edge distances, and sigmas), and calculate
// p(j|i) and then w_{ij}.
// [[Rcpp::export]]
arma::sp_mat distMatrixTowij(
  NumericVector is,
  NumericVector js,
  NumericVector xs,
  NumericVector sigmas,
  int N,
  Function callback
) {
  NumericVector rowSums = NumericVector(N);
  NumericVector pjis = NumericVector(is.length());
  for (int idx=0; idx < N; idx++) rowSums[idx] = 0;
  int i, j;
  double pji;
  // Compute pji, accumulate rowSums at the same time
#pragma omp parallel for shared(pjis, rowSums) private(pji, i)
  for (int e=0; e < pjis.length(); e++) {
    i = is[e];
    pji = exp(- pow(xs[e], 2)) / sigmas[i - 1];
    pjis[e] = pji;
#pragma omp atomic
    rowSums[i - 1] = rowSums[i - 1] + pji;
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
  // Now convert p(j|i) to w_{ij} by symmetrizing.
  arma::sp_mat wij = arma::sp_mat(N, N);
  for (int e=0; e < pjis.length(); e++) {
    int newi = is[e] - 1, newj = js[e] - 1;
    if (newi < newj) {
      int t = newi;
      newi = newj;
      newj = t;
    }
    double oldw = wij(newi, newj);
    wij(newi, newj) = (oldw > 0) ?
                                  (oldw / 2) + ((pjis[e] / rowSums[is[e] - 1]) / (2 * N )) :
                                  (pjis[e] / rowSums[is[e] - 1]) / N;
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
  return wij;
};

// [[Rcpp::export]]
void searchTree(int threshold, NumericVector indices,
                NumericMatrix data, NumericMatrix output,
                Function callback) {
  if (indices.length() <= threshold) {
    for (int i = 0; i < indices.length() - 1; i++) {
      for (int j = i + 1; j < indices.length(); j++) {
          output(i + j- 1, indices[i] - 1) = indices[j];
          output(i, indices[j] - 1) = indices[i];
      }
    }
    callback(indices.length());
    return;
  }
  // Get hyperplane
  NumericVector selections = RcppArmadillo::sample(indices, 2, false);
  NumericVector x1, x2, v, m;
  x1 = data.row(selections[0] - 1);
  x2 = data.row(selections[1] - 1);
  v = x2 - x1;
  m = (x1 + x2) / 2;
  double mv = sum(m * v); // This is the hyperplane
  NumericVector direction = NumericVector(indices.length());
  // Calculate distances from hyperplane - partition based on which distances are > 0
  for (int idx = 0; idx < indices.length();idx++) direction[idx] = sum(v * data.row(indices[idx] - 1)) - mv;
  int branch = 0;
  for (int i = 0; i < indices.length(); i++) if (direction[i] > 0) branch++;
  // Don't create branches that have only 2 nodes; if the split is that lopsided, recurse and try again
  if (branch < threshold / 3 || branch > indices.length() - (threshold / 3)) {
    searchTree(threshold, indices, data, output, callback);
    return;
  }
  // Recurse, attempt to work around segfault
  NumericVector left = NumericVector(branch);
  NumericVector right = NumericVector(indices.length() - branch);
  for (int leftidx = 0, rightidx = 0;
       leftidx + rightidx < indices.length(); ) {
    if (direction[leftidx + rightidx] > 0) left[leftidx++] = indices[leftidx + rightidx - 1];
    else right[rightidx++] = indices[leftidx + rightidx - 1];
  }
#pragma omp parallel sections
{
#pragma omp section
  searchTree(threshold, left, data, output, callback);
#pragma omp section
  searchTree(threshold, right, data, output, callback);
}
};

// [[Rcpp::export]]
double sigFunc(double sigma, NumericVector x_i, double perplexity) {
  NumericVector xs = exp(- pow(x_i, 2) / sigma);
  NumericVector softxs = xs / sum(xs);
  double p2 = - sum(log(softxs) / log(2)) / xs.length();
  return perplexity - p2;
};


