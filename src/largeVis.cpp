#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#include "progress.hpp"
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
arma::mat sgd(arma::mat coords,
         const NumericVector is, // vary randomly
         const NumericVector js, // ordered
         const NumericVector ps, // N+1 length vector of indices to start of each row j in vector is
         const NumericVector ws, // w{ij}
         const double gamma,
         const double rho,
         const double minRho,
         const bool useWeights,
         const int nBatches,
         const int M,
         const double alpha,
         bool verbose) {

  Progress progress(nBatches, verbose);

  const int updateRate = 1e6;
  const int D = coords.n_rows;
  const int N = ps.size() - 1;
  const int E = ws.length();
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
  if (! useWeights) {
    const double posScale = sum(ws);
    positiveEdgeWeights = arma::vec(E);
    positiveEdgeWeights[0] = ws[0] / posScale;
    for (int idx = 1; idx < E; idx++)
      positiveEdgeWeights[idx] = positiveEdgeWeights[idx - 1] + (ws[idx] / posScale);
  }

  Rcpp::List progressParams = Rcpp::List::create();
  double trailingLoss = 0;

  const int posSampleLength = min(1000000, nBatches);
  arma::vec positiveSamples = arma::randu<arma::vec>(posSampleLength);

  // Iterate through the edges in the positiveEdges vector
#pragma omp parallel for shared(coords, positiveSamples, trailingLoss) schedule(static)
  for (int eIdx=0; eIdx < nBatches; eIdx++) {
    if (progress.increment()) {
      const double posTarget = *(positiveSamples.begin() + (eIdx % posSampleLength));
      int k;
      int e_ij;
      if (useWeights) {
        e_ij = posTarget * (E - 1);
      } else {
        e_ij = std::distance(positiveEdgeWeights.begin(),
                             std::upper_bound(positiveEdgeWeights.begin(),
                                              positiveEdgeWeights.end(),
                                              posTarget));
      }
      const int i = is[e_ij];
      const int j = js[e_ij];

      const double localRho = rho - ((rho - minRho) * eIdx / nBatches);

      //if ((arma::randn<arma::vec>(1))[0] < 0) swap(i, j);

      // const arma::vec y_i = coords.row(i);
      // const arma::vec y_j = coords.row(j);
      arma::vec y_i = arma::vec(D);
      arma::vec y_j = arma::vec(D);
      for (int idx = 0; idx < D; idx++) {
        y_i[idx] = coords(idx,i);
        y_j[idx] = coords(idx,j);
      }

      // wij
      const double w = (useWeights) ? ws[e_ij] : 1;
      // TODO: RE-ADD TO USE EXP IF ALPHA = 0
      const double dist_ij = (alpha == 0) ? std::min(16.0, sqrt(dist(y_i, y_j))) : sqrt(dist(y_i, y_j));

      const arma::vec d_dist_ij = (y_i - y_j) / sqrt(dist_ij);
      double p_ij;
      if (alpha == 0)   p_ij =   1 / (1 +      exp(pow(dist_ij,2)));
      else              p_ij =   1 / (1 + (alpha * pow(dist_ij,2)));

      arma::vec d_p_ij;
      if (alpha == 0) d_p_ij =  - 2 * d_dist_ij * dist_ij * exp(pow(dist_ij,2)) / pow(1 + exp(pow(dist_ij,2)), 2);
      else            d_p_ij =  - 2 * d_dist_ij * dist_ij * alpha               / pow(1 +    (pow(dist_ij,2) * alpha),2);

      //double o = log(p_ij);
      const arma::vec d_j = (1 / p_ij) * d_p_ij;

      arma::vec samples = arma::randu<arma::vec>(M * 2);
      arma::vec::iterator targetIt = samples.begin();
      int sampleIdx = 1;
      // The indices of the nodes with edges to i
      arma::vec searchVector = i_idx.subvec(ps[i], ps[i + 1] - 1);
      arma::vec d_i = d_j;
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
        arma::vec y_k = arma::vec(D);
        for (int idx = 0; idx < D; idx++) y_k[idx] = coords(idx,k);

        double dist_ik = sqrt(dist(y_i, y_k));
        if (dist_ik == 0) continue; // Duplicates
        if (alpha == 0 && dist_ik > 16) dist_ik = 16; // Avoids numeric overflow

        const arma::vec d_dist_ik = (y_i - y_k) / sqrt(dist_ik);

        double p_ik;
        if (alpha == 0) p_ik  =  1 - (1 / (1 +      exp(pow(dist_ik,2))));
        else            p_ik  =  1 - (1 / (1 + (alpha * pow(dist_ik,2))));

        arma::vec d_p_ik;
        if (alpha == 0) d_p_ik = 2 * d_dist_ik * p_ik * exp(pow(dist_ik,2)) / pow(1 +      exp(pow(dist_ik,2)),2);
        else            d_p_ik = 2 * d_dist_ik * p_ik * alpha               / pow(1 + (alpha * pow(dist_ik,2)),2);

        //o += (gamma * log(p_ik));

        const arma::vec d_k = (gamma / p_ik) * d_p_ik;

        d_i += d_k;
        for (int idx = 0; idx < D; idx++) coords(idx,k) -= d_k[idx] * localRho * w;

        m++;
      }

      for (int idx = 0; idx < D; idx++) {
        coords(idx,j) -=  d_j[idx] * w * localRho;
        coords(idx,i) +=  d_i[idx] * w * localRho;
      }

      if (eIdx >0 && eIdx % posSampleLength == 0) positiveSamples.randu();
    }
  }
  return coords;
};




