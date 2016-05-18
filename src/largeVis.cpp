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

// The Euclidean distance between two vectors
inline double dist(NumericVector i, NumericVector j) {
  return sum(pow(i - j, 2));
}

inline double dist(arma::vec i, arma::vec j) {
  return sum(pow(i - j, 2));
}

/*
 * Some helper functions useful in debugging.
 */
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
 * The stochastic gradient descent function. Asynchronicity is enabled by openmp.
 */

// [[Rcpp::export]]
void sgd(NumericMatrix coords,
         NumericVector is, // vary randomly
         NumericVector js, // ordered
         NumericVector ps, // N+1 length vector of indices to start of each row j in vector is
         NumericVector ws, // w{ij}
         double gamma,
         double rho,
         double minRho,
         bool useWeights,
         int nBatches,
         int M,
         double alpha,
         Function callback) {

  int N = ps.size() - 1;
  const arma::vec i_idx = as<arma::vec>(is);
  // Calculate negative sample weights, d_{i}^0.75.
  // Stored as a vector of cumulative sums, normalized, so it can
  // be readily searched using binary searches.
  arma::vec negativeSampleWeights = pow(diff(ps), 0.75);
  const double scale = sum(negativeSampleWeights);
  negativeSampleWeights = negativeSampleWeights / scale;
  negativeSampleWeights = cumsum(negativeSampleWeights);

  // positive edges for sampling
  arma::vec positiveEdgeWeights;
  arma::vec positiveSamples;
  const int posSampleLength = min(1000000, nBatches);
  positiveSamples = arma::randu<arma::vec>(posSampleLength);
  if (! useWeights) {
    const double posScale = sum(ws);
    positiveEdgeWeights = as<arma::vec>(ws);
    positiveEdgeWeights = cumsum(positiveEdgeWeights / posScale);
  }

  // Iterate through the edges in the positiveEdges vector
#pragma omp parallel for shared(coords, positiveSamples)
  for (int eIdx=0; eIdx < nBatches; eIdx++) {
    arma::vec::iterator posIt = positiveSamples.begin();
    const double posTarget = posIt[eIdx % posSampleLength];
    int k;
    int e_ij;
    if (useWeights) {
      e_ij = posTarget * posSampleLength;
    } else {
      e_ij = std::distance(positiveEdgeWeights.begin(),
                           std::upper_bound(positiveEdgeWeights.begin(),
                                            positiveEdgeWeights.end(),
                                            posTarget));
    }
    int i = is[e_ij];
    int j = js[e_ij];

    const double localRho = rho - ((rho - minRho) * eIdx / nBatches);

    if ((arma::randn<arma::vec>(1))[0] < 0) swap(i, j);

    NumericVector y_i = coords.row(i);
    NumericVector y_j = coords.row(j);

    // wij
    const double w = (useWeights) ? ws[e_ij] : 1;
    // TODO: RE-ADD TO USE EXP IF ALPHA = 0
    const double dist_ij = sqrt(dist(y_i, y_j));

    const NumericVector d_dist_ij = (y_i - y_j) / sqrt(dist_ij);
    double p_ij;
    if (alpha == 0)   p_ij =   1 / (1 +      exp(pow(dist_ij,2)));
    else              p_ij =   1 / (1 + (alpha * pow(dist_ij,2)));
    NumericVector d_p_ij;
    if (alpha == 0) d_p_ij =   2 *  d_dist_ij * -exp(pow(dist_ij,2)) /    (1 + exp(pow(dist_ij,2)));
    else            d_p_ij =        d_dist_ij *     -pow(dist_ij,2)  / pow(1 +    (pow(dist_ij,2) * alpha),2);
 //   const double o_ij = log(p_ij);
    const NumericVector d_ij = (1 / p_ij) * d_p_ij;

    arma::vec samples = arma::randu<arma::vec>(M * 2);
    arma::vec::iterator targetIt = samples.begin();
    int sampleIdx = 1;
    // The indices of the nodes with edges to i
    arma::vec searchVector = i_idx.subvec(ps[i], ps[i + 1] - 1);
    NumericVector d_i = d_ij;
    int m = 0;
    while (m < M) {
      if (sampleIdx % (M * 2) == 0) samples.randu();
      // binary search implementing weighted sampling
      const double target = targetIt[sampleIdx++ % (M * 2)];
      int k;
      if (useWeights) k = target * (N - 1);
      else k = std::distance(negativeSampleWeights.begin(),
                             std::upper_bound(negativeSampleWeights.begin(),
                                              negativeSampleWeights.end(),
                                              target)
                               );

      if (k == i ||
          k == j ||
          sum(searchVector == k) > 0) continue;

      const NumericVector y_k = coords.row(k);
      const double dist_ik = sqrt(dist(y_i, y_k));
      if (dist_ik == 0) continue; // Duplicates

      const NumericVector d_dist_ik = (y_i - y_k) / sqrt(dist_ik);
      const double p_ik = 1 - (1 / (1 + (alpha * pow(dist_ik,2))));
      const NumericVector d_p_ik = d_dist_ik * (1 / ((pow(dist_ik,2) * pow(alpha, 2)) + alpha));
    //  const double o_ik = log(p_ik);
      const NumericVector d_ik = (gamma / p_ik) * d_p_ik;

      d_i = d_i + d_ik;
      coords(k,_) = coords(k,_) - (d_ik * localRho * w);
      m++;
    }
    coords(j,_) = coords(j,_) - (d_ij * w * localRho);
    coords(i,_) = coords(i,_) + (d_i  * w * localRho);
    if (eIdx > 0 && eIdx % 10000 == 0) callback(10000);
    if (eIdx >0 && eIdx % posSampleLength == 0) positiveSamples.randu();
  }
};




