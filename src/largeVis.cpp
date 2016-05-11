#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <math.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <algorithm>
#include <queue>
#include <vector>
#include <set>
using namespace Rcpp;
using namespace std;

struct heapObject {
  double d;
  int n;

  heapObject(double d, int n) : d(d), n(n) {}

  bool operator<(const struct heapObject& other) const {
    return d < other.d;
  }
};

double dist(NumericVector i, NumericVector j) {
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

arma::sp_mat convertSparse(S4 mat) {
  IntegerVector dims = mat.slot("Dim");
  arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
  arma::urowvec j  = Rcpp::as<arma::urowvec>(mat.slot("j"));
  arma::umat locs;
  locs.insert_rows(0, i);
  locs.insert_rows(0,j);
  arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));
  int nrow = dims[0], ncol = dims[1];
  arma::sp_mat res(locs, x, nrow, ncol, true, true);
  return res;
};

// [[Rcpp::export]]
void sgd(NumericMatrix coords,
         NumericVector positiveEdges,
         NumericVector is,
         NumericVector js,
         NumericVector ws,
         NumericVector negativeSampleWeights,
         int gamma,
         int rho,
         int minRho,
         bool useWeights,
         int M,
         int alpha,
         Function callback) {

  int i,j,e_ij,k;
  double w = 1;
  IntegerVector availjs, negjs;
  NumericVector y_i, y_j, grads, j_row;
  // Make sparse matrix
  arma::umat locs; // = join_cols(Rcpp::as<arma::urowvec> (is), Rcpp::as<arma::urowvec> (js)).t();
  arma::urowvec isvec = Rcpp::as<arma::urowvec> (is);
  arma::urowvec jsvec = Rcpp::as<arma::urowvec> (js);
  locs.insert_rows(0,isvec - 1);
  locs.insert_rows(1,jsvec - 1);
  int N = max(is);

  arma::sp_mat wij(locs, Rcpp::as<arma::vec> (ws), N, N, true, true);

  #pragma omp parallel for shared(coords, rho) private(i,j,e_ij,y_i,y_j,grads,j_row,availjs,negjs,k)
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
    if (alpha != 0)   grads =   - w * 2 * (y_i - y_j) * alpha                 / ((alpha * dist(y_i, y_j)) + 1);
    else              grads =   - w * 2 * (y_i - y_j) * exp(dist(y_i, y_j))   / (exp( dist(y_i, y_j)) + 1);
    y_i = y_i + (grads * rho);
    coords(j,_) = y_j - (grads * rho);
    j_row = wij.row(i);
    availjs = NumericVector::create();
    for (int jidx = 0; jidx < j_row.length(); jidx++) {
      if (jidx != i && j_row[jidx] == 0) availjs.push_front(jidx);
    }
    negjs = RcppArmadillo::sample(availjs, (M < availjs.length()) ? M : availjs.length(), false, negativeSampleWeights[availjs]);

    coords(i,_) = y_i + (grads * rho);
    for (int jidx = 0; jidx < negjs.length(); jidx++) {
      k = negjs[jidx] - 1;
      y_j = coords.row(k);
      if (alpha != 0 )  grads =     w * gamma * (y_i - y_j) * 2 / (dist(y_i, y_j) * ( 1 + (alpha * dist(y_i, y_j))));
      else              grads =     w * gamma * (y_i - y_j)     / (exp(dist(y_i, y_j)) + 1);
//      coords(i,_) = y_i + (grads * rho / negjs.length());
      coords(k,_) = y_j - (grads * rho);
    }
    rho = rho - ((rho - minRho) / (positiveEdges.length() + 1));
    if (eIdx > 0 && eIdx % 1000 == 0) callback(1000);
  }
};

// [[Rcpp::export]]
void distance(NumericVector is, NumericVector js, NumericVector xs, NumericMatrix data) {
  for (int i=0; i < is.length(); i++) {
    xs[i] = sqrt(sum(pow(data.row(is[i] - 1) - data.row(js[i] - 1), 2)));
  }
};

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
#pragma omp parallel for shared(pjis, rowSums) private(pji, i)
  for (int e=0; e < n_e; e++) {
    i = is[e];
    pji = exp(- pow(xs[e], 2)) / sigmas[i - 1];
    pjis[e] = pji;
    rowSums[i - 1] = rowSums[i - 1] + pji;
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
#pragma omp parallel for shared(pjis, rowSums) private(pji, i)
  for (int e=0; e < n_e; e++) {
    i = js[e];
    pji = exp(-pow(xs[e], 2)) / sigmas[i - 1];
    pjis[e + n_e] = pji;
    rowSums[i - 1] = rowSums[i - 1] + pji;
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
#pragma omp parallel for shared(pjis, rowSums) private(outVector)
  for (int e=0; e < n_e; e++) {
    outVector[e] = ((pjis[e] / rowSums[is[e] - 1]) + (pjis[e] / rowSums[js[e] - 1])) / ( 2 * N );
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
};

// [[Rcpp::export]]
void searchTree(int threshold, NumericVector indices,
                NumericMatrix data, NumericMatrix output, Function callback) {
  if (indices.length() == 1) {
    output(0, indices[1] - 1) = indices[1];
  }
  if (indices.length() <= threshold) {
    for (int i = 0; i < indices.length(); i++) {
      for (int j = i + 1; j < indices.length(); j++) {
          output(i, indices[i] - 1) = indices[j];
          output(i, indices[j] - 1) = indices[i];
      }
    }
    callback(indices.length());
    return;
  }

  NumericVector selections = RcppArmadillo::sample(indices, 2, false);
  NumericVector x1, x2, v, m;
  x1 = data.row(selections[0] - 1);
  x2 = data.row(selections[1] - 1);
  v = x2 - x1;
  m = (x1 + x2) / 2;
  double mv = sum(m * v);
  NumericVector direction = NumericVector(indices.length());
  for (int idx = 0; idx < indices.length();idx++) direction[idx] = sum(v * data.row(indices[idx] - 1)) - mv;
  double branch = sum(direction > 0);
  if (branch < 3 || branch > indices.length() - 3) {
    searchTree(threshold, indices, data, output, callback);
    return;
  }
  searchTree(threshold, indices[direction > 0], data, output, callback);
  searchTree(threshold, indices[direction <= 0], data, output, callback);
};
