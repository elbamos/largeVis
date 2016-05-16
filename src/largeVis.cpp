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

inline double dist(arma::vec i, arma::vec j) {
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
        if (j == 0) continue;
        if (j - 1 != i && seen.insert(j - 1).second) {
          d = dist(x_i, data.row(j-1));
          heap.push(heapObject(d, j));
          if (heap.size() > K) heap.pop();
        }

        for (int kidx = 0; kidx < oldK; kidx++) {
          k = old_knns.column(j - 1)[kidx];
          if (k == 0) continue;
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

void checkVector(const arma::vec& x,
                 const std::string& label) {
  if (x.has_nan() || x.has_inf())
    Rcout << "\n Failure at " << label;
};

double objective(const arma::mat& inputs, double gamma, double alpha) {
  double objective = log(1 / (1 + (alpha * dist(inputs.col(0), inputs.col(1)))));
  for (int i = 2; i < 7; i++)
    objective += gamma * log(1 - (1 / (1 + alpha * dist(inputs.col(0), inputs.col(i)))));
  return objective;
}

void checkGrad(const arma::vec& x,
               const arma::vec& y,
               const arma::vec& grad,
               bool together,
               const string& label) {
  double oldDist = dist(x,y);
  double newDist = dist(x + grad, y - grad);
  if (together && newDist > oldDist) Rcout << "\nGrad " << label << " yi " << x << " other " << y << " grad " << grad << " moved further apart.";
  else if (! together && newDist < oldDist) Rcout << "\nGrad " << label << " yi " << x << " other " << y << " grad " << grad << "moved closer together.";
};
/*
 * w * (log( 1 / (1 + (a * sqrt((b - j)^2 + (c - f)^2)^2))) + g(log(1 - (1 / (1 + (a * sqrt((b - g)^2 + (c - k)^2)^2)))) + log(1 - (1 / (1 + (a * sqrt((b - m)^2 + (c - n)^2)^2))))))
 */
// [[Rcpp::export]]
arma::mat sgd(NumericMatrix coords,
         NumericVector positiveEdges,
         NumericVector is, // vary randomly
         NumericVector js, // ordered
         NumericVector ps, // N+1 length vector of indices to start of each row j in vector is
         NumericVector ws, // w{ij}
         double gamma,
         double rho,
         double minRho,
         bool useWeights,
         int M,
         double alpha,
         Function callback) {

  arma::vec p = as<arma::vec>(ps);
  arma::vec i_idx = as<arma::vec>(is);
  arma::vec j_idx = as<arma::vec>(js);
  arma::vec pEdges = as<arma::vec>(positiveEdges);
  arma::mat coordinates = as<arma::mat>(coords);
  int N = p.size() - 1;
  // Calculate negative sample weights, d_{i}^0.75.
  // Stored as a vector of cumulative sums, normalized, so it can
  // be readily searched using binary searches.
  // TODO:  CREATE A FORM FOR WHEN USEWEIGHTS = TRUE
  arma::vec negativeSampleWeights = pow(diff(p), 0.75);
  double scale = sum(negativeSampleWeights);
  negativeSampleWeights = negativeSampleWeights / scale;
  negativeSampleWeights = cumsum(negativeSampleWeights);

  // arma::vec losses = arma::vec(pEdges.size()).zero();

  // Iterate through the edges in the positiveEdges vector
#pragma omp parallel for
  for (int eIdx=0; eIdx < pEdges.size(); eIdx++) {
    const int e_ij = pEdges[eIdx];
    int i = i_idx[e_ij];
    int j = j_idx[e_ij];

    const double localRho = rho - ((rho - minRho) * eIdx / (pEdges.size()));

    if ((arma::randn<arma::vec>(1))[0] < 0) swap(i, j);

    arma::vec y_i = coordinates.col(i);
    arma::vec y_j = coordinates.col(j);

    // arma::mat before = arma::mat(2,M + 2);
    // arma::mat after = arma::mat(2,M + 2);
    // before.col(0) = y_i;
    // before.col(1) = y_j;
    // wij
    const double w = (useWeights) ? ws[e_ij] : 1;
// should be negative below
    const double dist_ij = sqrt(dist(y_i, y_j));
    const arma::vec d_dist_ij = (y_i - y_j) / dist_ij;
    const double p_ij = 1 / (1 + (alpha * pow(dist_ij,2)));
    const arma::vec d_p_ij = d_dist_ij * -pow(dist_ij, 2) / pow((pow(dist_ij,2) * alpha) + 1,2);
 //   const double o_ij = log(p_ij);
    const arma::vec d_ij = (1 / p_ij) * d_p_ij;
   // arma::vec igrad = -2 * alpha * (y_i - y_j) / (1 + (alpha * dist(y_i,y_j)));
  //  checkGrad(y_i, y_j, igrad, true, "pos");
#pragma omp critical
{
    coordinates.col(j) -= (w * localRho * d_ij);
    // after.col(1) = coordinates.col(j);
}
    // checkVector(igrad, "ijgrad");

    arma::vec samples = arma::randu<arma::vec>(M * 2);
    arma::vec::iterator targetIt = samples.begin();
    int sampleIdx = 1;
    int m = 0;
    // int ms[5];
    // The indices of the nodes with edges to i
    arma::vec searchVector = i_idx.subvec(p[i], p[i + 1] - 1);
    arma::vec d_i = d_ij;
    while (m < M) {
      if (sampleIdx % (M * 2) == 0) samples.randu();
      // binary search for lowest number greater than the sampled number
      const double target = targetIt[sampleIdx++ % (M * 2)];
      int k;
      if (useWeights) k = target * (N - 1);
      else {
        arma::vec::iterator loc = std::upper_bound(negativeSampleWeights.begin(),
                                                 negativeSampleWeights.end(), target);
        k = std::distance(negativeSampleWeights.begin(), loc);
      }
      if (k == i ||
          k == j ||
          sum(searchVector == k) > 0) continue;
      // ms[m] = k;
      // Calculate gradient for a single negative sample
      arma::vec y_k = coordinates.col(k);
      const double dist_ik = sqrt(dist(y_i, y_k));
      const arma::vec d_dist_ik = (y_i - y_k) / dist_ik;
      const double p_ik = 1 - (1 / (1 + (alpha * pow(dist_ik,2))));
      const arma::vec d_p_ik = d_dist_ik * (1 / ((pow(dist_ik,2) * pow(alpha, 2)) + alpha));
    //  const double o_ik = log(p_ik);
      const arma::vec d_ik = (gamma / p_ik) * d_p_ik;


      //
      // const double alphadk1 = (alpha * dist(y_i, y_k)) + 1;
      // arma::vec gradk = 2 * gamma * (y_i - y_k) / (dist(y_i, y_k) * alphadk1);
  //    checkGrad(y_i, y_k, igrad, false, "neg");
      d_i = d_i + d_ik;
      // before.col(m + 2) = y_k;
#pragma omp critical
{
      coordinates.col(k) -=  (d_ik * localRho * w);
}
// after.col(m + 2) = coordinates.col(k);
m++;


    }
    // checkVector(igrad, "jgrad");

#pragma omp critical
{
    coordinates.col(i) += (d_i * w * localRho);
}
   // after.col(0) = coordinates.col(i);
    // double bef = objective(before, gamma, alpha);
    // double aft = objective(after, gamma, alpha);
    // if (bef >= aft) Rcout << eIdx << " BAD " << bef << " to " << aft << " i " << i << " j " << j << " k " << ms << "\n";
    //
    // if (bef >= aft && eIdx % 10 == 0) {
    //   for (int idx = 0; idx < M + 2; idx++) {
    //     Rcout << before.col(idx)[0] << " " << before.col(idx)[1] << " aft " << after.col(idx)[0] << " " <<
    //       after.col(idx)[1] << "\n";
    //   }
    // }
    if (eIdx > 0 && eIdx % 10000 == 0) callback(10000);
  }
  return coordinates;
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

void searchTree(int threshold,
                const arma::uvec& indices,
                const arma::mat& data,
                arma::imat& output,
                Function callback) {
  if (indices.size() < threshold) {
  #pragma omp critical
  {
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
  }
  callback(indices.size());
  return;
  }
  // Get hyperplane
  arma::uvec selections = indices.elem(arma::randi<arma::uvec>(2, arma::distr_param(0, indices.size() - 1)));
  arma::vec v =  data.col(selections[1]) - data.col(selections[0]);
  arma::vec m = sum(data.cols(selections), 1) / 2;
  double mv = dot(m,v); // This is the hyperplane
  arma::vec direction = arma::vec(indices.size());
  for (int i = 0; i < indices.size(); i++) direction[i] = sum(data.col(indices[i]) % v) - mv;
  int branch = sum(direction > 0);
  // Don't create branches that have only 2 nodes; if the split is that lopsided, recurse and try again
  if (branch < threshold / 3 || branch > indices.size() - (threshold / 3)) {
    searchTree(threshold, indices, data, output, callback);
    return;
  }
//  #pragma omp parallel sections
  {
//    #pragma omp section
    {
      searchTree(threshold, indices.elem(arma::find(direction > 0)), data, output, callback);
    }
//    #pragma omp section
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
  #pragma omp parallel for
  for (int t = 0; t < n_trees; t++) {
    arma::imat output = arma::imat(threshold + 1, inputData.n_cols);
    output.fill(-1);
    searchTree(threshold,
               arma::regspace<arma::uvec>(0, inputData.n_cols - 1),
               inputData,
               output,
               callback);
    #pragma omp critical
    {
        for (int i = 0; i < inputData.n_cols; i++) for (int j = 0; j <= threshold; j++)  heap[i].insert(output(j,i));
    }
  }
  int sz = 0;
  for (int i = 0; i < heap.size(); i++) if (heap[i].size() > sz) sz = heap[i].size();
  arma::mat output = arma::mat(sz, inputData.n_cols).zeros();
  for (int i = 0; i < heap.size(); i++) {
    std::set<int> neighbors = heap[i];
    int j = 0;
    for (std::set<int>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
      if (*it == 0) continue;
      output(j, i) = *it;
      j++;
    }
  }
  return output;
}


// [[Rcpp::export]]
double sigFunc(double sigma, NumericVector x_i, double perplexity) {
  NumericVector xs = exp(- pow(x_i, 2) / sigma);
  NumericVector softxs = xs / sum(xs);
  double p2 = - sum(log(softxs) / log(2)) / xs.length();
  return perplexity - p2;
};


