// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.hpp"

double relDist(const arma::vec& i, const arma::vec& j) {
  const int lim = i.n_elem;
  double cnt = 0;
  for (int idx = 0; idx < lim; idx++) cnt += ((i[idx] - j[idx]) * (i[idx] - j[idx]));
  return cnt;
}
// Vanilla euclidean
double dist(const arma::vec& i, const arma::vec& j) {
  return sqrt(relDist(i,j));
}

// Copies the vector sums into a vector while it computes distance^2 -
// useful in calculating the gradients during SGD
double distAndVector(double *x_i,
                     double *x_j,
                     double *output,
                     const int& D) {
  double cnt = 0;
  for (int d = 0; d < D; d++) {
    double t = x_i[d] - x_j[d];
    output[d] = t;
    cnt += t * t;
  }
  return cnt;
}

// Vanilla cosine distance calculation
double cosDist(const arma::vec& i, const arma::vec& j) {
  int D = i.n_elem;
  double pp = 0, qq = 0, pq = 0;
  for (int d = 0; d < D; d++) {
    pp += (i[d]) * (i[d]);
    qq += (j[d]) * (j[d]);
    pq += (i[d]) * (j[d]);
  }
  double ppqq = pp * qq;
  if (ppqq > 0) return 2.0 - 2.0 * pq / sqrt(ppqq);
  else return 2.0; // cos is 0
}
// Versions of the distance functions for finding the neighbors
// of sparse matrices.  Not optimized.
double sparseDist(const sp_mat& i, const sp_mat& j) {
  return as_scalar(sqrt(sum(square(i - j))));
}
double sparseCosDist(const sp_mat& i, const sp_mat& j) {
  return 2.0 - 2.0 * (as_scalar((dot(i,j)) / as_scalar(norm(i,2) * norm(j,2))));
}
double sparseRelDist(const sp_mat& i, const sp_mat& j) {
  return as_scalar(sum(square(i - j)));
}

/*
 * Fast calculation of pairwise distances with the result stored in a pre-allocated vector.
 */
// [[Rcpp::export]]
arma::vec fastDistance(const NumericVector is,
                       const NumericVector js,
                       const arma::mat& data,
                       const std::string& distMethod,
                       bool verbose) {

  Progress p(is.size(), verbose);
  vec xs = vec(is.size());
  double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = dist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = cosDist;
#ifdef _OPENMP
#pragma omp parallel for shared (xs)
#endif
  for (int i=0; i < is.length(); i++) if (p.increment()) xs[i] =
    distanceFunction(data.col(is[i]), data.col(js[i]));
  return xs;
};

arma::vec fastSparseDistance(const arma::vec& is,
                             const arma::vec& js,
                             const sp_mat& data,
                             const std::string& distMethod,
                             bool verbose) {

  Progress p(is.size(), verbose);
  vec xs = vec(is.size());
  double (*distanceFunction)(
      const sp_mat& x_i,
      const sp_mat& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = sparseDist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = sparseCosDist;
#ifdef _OPENMP
#pragma omp parallel for shared (xs)
#endif
  for (int i=0; i < is.size(); i++) if (p.increment()) xs[i] =
    distanceFunction(data.col(is[i]), data.col(js[i]));
  return xs;
};

// [[Rcpp::export]]
arma::vec fastCDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& p_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose) {
  const int N = p_locations.size() - 1;
  const sp_mat data = sp_mat(i_locations, p_locations, x, N, N);
  return fastSparseDistance(is,js,data,distMethod,verbose);
}

// [[Rcpp::export]]
arma::vec fastSDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& j_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose) {
  const umat locations = join_cols(i_locations, j_locations);
  const sp_mat data = sp_mat(locations, x);
  return fastSparseDistance(is,js,data,distMethod,verbose);
}
