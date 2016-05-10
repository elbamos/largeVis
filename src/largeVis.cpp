#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <algorithm>
#include <queue>
#include <vector>
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

void pushHeap(std::priority_queue<heapObject> heap, int j, NumericVector x_i, NumericMatrix data, int K) {
  double d = sqrt(dist(x_i, data.row(j-1)));
  heap.push(heapObject(d, j));
  if (heap.size() > K) heap.pop();
}

// [[Rcpp::export('neighbors_inner')]]
NumericMatrix neighbors_inner(int K, NumericMatrix old_knns, NumericMatrix data, Function callback) {
  int N = old_knns.ncol();
  int oldK = old_knns.nrow();
  int dataD = data.ncol();
  double d;
  int j, k;

  NumericMatrix output = NumericMatrix(K, N);
  //for (int i = 0; i < K; i++) for (int j = 0; j < N; j++) output(i,j) = 0;
  std::fill(output.begin(), output.end(), 0);

  NumericVector x_i;

  for (int i = 0; i < N; i++) {
    if (i % 1000 == 999) callback(1000);
    x_i = data.row(i);
    std::priority_queue<heapObject> heap;
    #pragma omp parallel for
    for (int jidx = 0; jidx < oldK; jidx++) {
      j = old_knns.column(i)[jidx];
      if (j == 0) break;
      pushHeap(heap, j, x_i, data, K);

      for (int kidx = 0; kidx < oldK; kidx++) {
        k = old_knns.column(j - 1)[kidx];
        if (k == 0) break;
        pushHeap(heap, k, x_i, data, K);
      }
      int j = 0;
      while (j < K && ! heap.empty()) {
        output(j, i) = heap.top().n;
        heap.pop();
        j++;
      }
    }
  }
  return (output);
};

arma::sp_mat convertSparse(S4 mat) {
  IntegerVector dims = mat.slot("Dim");
  arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
  arma::urowvec j  = Rcpp::as<arma::urowvec>(mat.slot("j"));
  arma::umat locs;
  Rcout << locs.n_cols << " " << locs.n_rows << "\n";
  arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));
  int nrow = dims[0], ncol = dims[1];
  arma::sp_mat res(locs, x, nrow, ncol, true, true);
  return res;
};

// [[Rcpp::export('sgd')]]
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
         S4 wij,
         int M,
         int alpha,
         Function callback) {

    bool useAlpha = true;

  int i,j,e_ij;
  int w = 1;
  IntegerVector availjs, negjs;
  NumericVector y_i, y_j, grads, j_row;
  arma::sp_mat mat = convertSparse(wij);

  for (int eIdx=0; eIdx < positiveEdges.length(); eIdx++) {
    e_ij = positiveEdges[eIdx];
    i = is[e_ij];
    j = js[e_ij];
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

    j_row = mat.row(i);
    availjs = NumericVector::create();
    for (int jidx = 0; jidx < j_row.length(); jidx++) {
      if (jidx != i && j_row[jidx] == 0) availjs.push_front(jidx);
    }
    negjs = RcppArmadillo::sample(availjs, (M < availjs.length()) ? M : availjs.length(), false, negativeSampleWeights[availjs]);
    for (int jidx = 0; jidx < negjs.length(); jidx++) {
      y_j = coords.row(negjs[jidx]);
      if (alpha != 0 )  grads =     w * gamma * (y_i - y_j) * 2                   / (dist(y_i, y_j) * ( 1 + (alpha * dist(y_i, y_j))));
      else              grads =     w * gamma * (y_i - y_j)                       / (exp(dist(y_i, y_j)) + 1);
      coords(i,_) = y_i + (grads * rho / negjs.length());
      coords(j,_) = y_j - (grads * rho);
    }
    rho = rho - ((rho - minRho) / (positiveEdges.length() + 1));
    if (eIdx % 1000 == 0) callback(1000);
  }
};
