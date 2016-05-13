#include <omp.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <algorithm>
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

    #pragma omp parallel for shared(nextKnns)
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
 * Note three changes from the paper:
 * 1. \rho is set to decline non-linearly
 * 2. The gradient applied to update y_i for the positive sample is doubled (2 * \rho)
 *
 *
 */

// [[Rcpp::export]]
void sgd(NumericMatrix coords,
         NumericVector positiveEdges,
         NumericVector is,
         NumericVector js,
         NumericVector ws, // w{ij}
         int gamma,
         int rho,
         int minRho,
         bool useWeights,
         int M,
         int alpha,
         Function callback) {

  NumericVector negativeSampleWeights = NumericVector(is.length());
  NumericVector ps = NumericVector(is.length() + 1);
  int jidx = -1;
  // Cache a vector of negative weights for sampling
  // At the same time, build a p-vector of indices to the start of each row (column)
  for (int eidx = 0; eidx < negativeSampleWeights.length(); eidx++) {
    negativeSampleWeights[eidx] = pow(ps[eidx + 1] - ps[eidx], 0.75);
    if (is[eidx] != jidx) {
      ps[++jidx] = eidx + 1;
    }
  }
  // Normalize the negative sample weights.  This will allow us to sample without looping
  // through the whole vector.
  negativeSampleWeights = negativeSampleWeights / sum(negativeSampleWeights);

  int i,j,e_ij,k;
  double w = 1;
  IntegerVector availjs, negjs;
  NumericVector y_i, y_j, grads;

  #pragma omp parallel for shared(coords, rho) private(i,j,e_ij,y_i,y_j,grads,availjs,negjs,k)
  for (int eIdx=0; eIdx < positiveEdges.length(); eIdx++) {
    e_ij = positiveEdges[eIdx] - 1;
    i = is[e_ij] - 1;
    j = js[e_ij] - 1;

    if (rnorm(1)[1] < 0) {
      int t = i;
      i = j;
      j = t;
    }
    y_i = coords.row(i);
    y_j = coords.row(j);
    // Calculate gradients
    if (alpha != 0)   grads =   - w * 2 * (y_i - y_j) * alpha                 / ((alpha * dist(y_i, y_j)) + 1);
    else              grads =   - w * 2 * (y_i - y_j) * exp(dist(y_i, y_j))   / (exp( dist(y_i, y_j)) + 1);
    // Update parameters
    y_i = y_i + (grads * rho * 2);
    coords(j,_) = y_j - (grads * rho);
    // Negative samples
    NumericVector samples = runif(M).sort();
    double maxweight = 1;
    int jidxb = ps[i], jidxe = ps[i + 1]; // position in the i, j * x vectors of the first element in column i
    int nidx = 0; // count of negative elements found so far, index to next negative element to fill
    // subtract from maxweight the weights of each present edge, so the non-present ones
    // can be rescaled without having to loop through twice
    for (int idx = ps[i]; idx < ps[i + 1]; idx++) maxweight = maxweight - negativeSampleWeights[idx];
    double runningCount = 0;
    int kidx = 0;
    int k;
    do {
      // If this is one of the negative indices, skip it
      if (kidx >= jidxb && kidx < jidxe ) continue;
      runningCount += negativeSampleWeights[kidx] / maxweight;
      if (runningCount > samples[nidx]) {
        k = js[kidx] - 1;
        y_j = coords.row(k);
        if (alpha != 0 )  grads =     w * gamma * (y_i - y_j) * 2 / (dist(y_i, y_j) * ( 1 + (alpha * dist(y_i, y_j))));
        else              grads =     w * gamma * (y_i - y_j)     / (exp(dist(y_i, y_j)) + 1);
        y_i = y_i + (grads * rho / negjs.length());
        coords(k,_) = y_j - (grads * rho);
        nidx++;
      }
    } while (nidx < M && ++kidx < ps.length());

    coords(i,_) = y_i;

    rho = rho - ((rho - minRho) / (positiveEdges.length() + 1));
    if (eIdx > 0 && eIdx % 1000 == 0) callback(5000);
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
void distMatrixTowij(
  NumericVector is,
  NumericVector js,
  NumericVector xs,
  NumericVector sigmas,
  NumericVector outVector,
  int N,
  Function callback
) {
  int n_e = is.length();
  NumericVector rowSums = NumericVector(N);
  NumericVector pjis = NumericVector(n_e* 2);
  for (int idx=0; idx < N; idx++) rowSums[idx] = 0;
  int i, j;
  double pji;
  // We have to compute pji twice, once for the lower and once for the upper triangle
  // We can accumulate rowSums at the same time
#pragma omp parallel for shared(pjis, rowSums) private(pji, i)
  for (int e=0; e < n_e; e++) {
    i = is[e];
    pji = exp(- pow(xs[e], 2)) / sigmas[i - 1];
    pjis[e] = pji;
#pragma omp atomic
    rowSums[i - 1] = rowSums[i - 1] + pji;
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
#pragma omp parallel for shared(pjis, rowSums) private(pji, i)
  for (int e=0; e < n_e; e++) {
    i = js[e];
    pji = exp(-pow(xs[e], 2)) / sigmas[i - 1];
    pjis[e + n_e] = pji;
#pragma omp atomic
    rowSums[i - 1] = rowSums[i - 1] + pji;
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
  // Now convert p(j|i) to w_{ij} by symmetrizing.
#pragma omp parallel for shared(pjis, rowSums, outVector)
  for (int e=0; e < n_e; e++) {
    outVector[e] = ((pjis[e] / rowSums[is[e] - 1]) + (pjis[e] / rowSums[js[e] - 1])) / ( 2 * N );
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
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
