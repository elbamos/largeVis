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

    // wij
    const double w = (useWeights) ? ws[e_ij] : 1;
    // TODO: RE-ADD TO USE EXP IF ALPHA = 0
    const double dist_ij = sqrt(dist(y_i, y_j));
    const arma::vec d_dist_ij = (y_i - y_j) / dist_ij;
    const double p_ij = //(alpha == 0) ?  (1 / (1 + exp(pow(dist_ij,2)))):
                                        1 / (1 + (alpha * pow(dist_ij,2)));
    arma::vec d_p_ij;
    if (alpha == 0) d_p_ij = - 2 *  d_dist_ij * dist_ij * exp(pow(dist_ij,2)) / (1 + exp(pow(dist_ij,2)));
    else            d_p_ij =        d_dist_ij * -pow(dist_ij, 2) / pow((pow(dist_ij,2) * alpha) + 1,2);
 //   const double o_ij = log(p_ij);
    const arma::vec d_ij = (1 / p_ij) * d_p_ij;
    #pragma omp critical
    {
      coordinates.col(j) -= (w * localRho * d_ij);
    }
    arma::vec samples = arma::randu<arma::vec>(M * 2);
    arma::vec::iterator targetIt = samples.begin();
    int sampleIdx = 1;
    int m = 0;
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

      arma::vec y_k = coordinates.col(k);
      const double dist_ik = sqrt(dist(y_i, y_k));
      if (dist_ik == 0) continue;
      const arma::vec d_dist_ik = (y_i - y_k) / dist_ik;
      const double p_ik = 1 - (1 / (1 + (alpha * pow(dist_ik,2))));
      const arma::vec d_p_ik = d_dist_ik * (1 / ((pow(dist_ik,2) * pow(alpha, 2)) + alpha));
    //  const double o_ik = log(p_ik);
      const arma::vec d_ik = (gamma / p_ik) * d_p_ik;

      d_i = d_i + d_ik;

      #pragma omp critical
      {
        coordinates.col(k) -=  (d_ik * localRho * w);
      }
      m++;
    }
    #pragma omp critical
    {
        coordinates.col(i) += (d_i * w * localRho);
    }
    if (eIdx > 0 && eIdx % 10000 == 0) callback(10000);
  }
  return coordinates;
};




