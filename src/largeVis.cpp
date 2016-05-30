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

inline double dist(arma::vec i, arma::vec j) {
  return sum(square(i - j));
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
              const arma::vec& is, // vary randomly
              const NumericVector js, // ordered
              const NumericVector ps, // N+1 length vector of indices to start of each row j in vector is
              const NumericVector ws, // w{ij}
              const double gamma,
              const double rho,
              const double minRho,
              const bool useWeights,
              const long nBatches,
              const int M,
              const double alpha,
              bool verbose) {

  Progress progress(nBatches, verbose);

  const int D = coords.n_rows;
  const int N = ps.size() - 1;
  const int E = ws.length();
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

  const int posSampleLength = ((nBatches > 1000000) ? 1000000 : (int) nBatches);
  arma::vec positiveSamples = arma::randu<arma::vec>(posSampleLength);

  // Iterate through the edges in the positiveEdges vector
#pragma omp parallel for shared(coords, positiveSamples) schedule(static)
  for (long eIdx=0; eIdx < nBatches; eIdx++) {
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

      const arma::vec y_i = coords.col(i);
      const arma::vec y_j = coords.col(j);

      // wij
      const double w = (useWeights) ? ws[e_ij] : 1;

      const double dist_ij = dist(y_i, y_j);

      const arma::vec d_dist_ij = (y_i - y_j) / sqrt(dist_ij);
      double p_ij;
      if (alpha == 0)   p_ij =   1 / (1 +      exp(dist_ij));
      else              p_ij =   1 / (1 + (alpha * dist_ij));

      arma::vec d_p_ij;
      if (alpha == 0) d_p_ij =  d_dist_ij * -2 * dist_ij * exp(dist_ij) / pow(1 + exp(dist_ij), 2);
      else            d_p_ij =  d_dist_ij * -2 * dist_ij * alpha        / pow(1 +    (dist_ij * alpha),2);

      //double o = log(p_ij);
      const arma::vec d_j = (1 / p_ij) * d_p_ij;
      // alternative: d_i - 2 * alpha * (y_i - y_j) / (alpha * sum(square(y_i - y_j)))

      arma::vec samples = arma::randu<arma::vec>(M * 2);
      arma::vec::iterator targetIt = samples.begin();
      int sampleIdx = 1;
      // The indices of the nodes with edges to i
      arma::vec searchVector = is.subvec(ps[i], ps[i + 1] - 1);
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
        const arma::vec y_k = coords.col(k);

        const double dist_ik = dist(y_i, y_k);
        if (dist_ik == 0) continue; // Duplicates

        const arma::vec d_dist_ik = (y_i - y_k) / sqrt(dist_ik);

        double p_ik;
        if (alpha == 0) p_ik  =  1 - (1 / (1 +      exp(dist_ik)));
        else            p_ik  =  1 - (1 / (1 + (alpha * dist_ik)));

        arma::vec d_p_ik;
        if (alpha == 0) d_p_ik =  d_dist_ik * 2 * dist_ik * exp(dist_ik) / pow(1 +      exp(dist_ik),2);
        else            d_p_ik =  d_dist_ik * 2 * dist_ik * alpha        / pow(1 + (alpha * dist_ik),2);
        //o += (gamma * log(p_ik));

        const arma::vec d_k = (gamma / p_ik) * d_p_ik;
        // alternative:  d_k = 2 * alpha * (y_i - y_k) / (square(1 + (alpha * sum(square(y_i - y_k)))) * (1 - (1 / (alpha * sum(square(y_i - y_k))))))

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
