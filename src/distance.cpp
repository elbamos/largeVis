// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

distancetype relDist(const arma::vec& i, const arma::vec& j) {
  const dimidxtype D = i.n_elem;
  distancetype cnt = 0;
  for (dimidxtype idx = 0; idx < D; idx++) cnt += ((i[idx] - j[idx]) * (i[idx] - j[idx]));
  return cnt;
}
// Vanilla euclidean
distancetype dist(const arma::vec& i, const arma::vec& j) {
  return sqrt(relDist(i,j));
}

// Vanilla cosine distance calculation
distancetype cosDist(const arma::vec& i, const arma::vec& j) {
  const dimidxtype D = i.n_elem;
	distancetype pp = 0, qq = 0, pq = 0;
  for (dimidxtype d = 0; d < D; d++) {
    pp += (i[d]) * (i[d]);
    qq += (j[d]) * (j[d]);
    pq += (i[d]) * (j[d]);
  }
  distancetype ppqq = pp * qq;
  if (ppqq > 0) return 2.0 - 2.0 * pq / sqrt(ppqq);
  else return 2.0; // cos is 0
}
// Versions of the distance functions for finding the neighbors
// of sparse matrices.  Not optimized.
distancetype sparseDist(const sp_mat& i, const sp_mat& j) {
  return as_scalar(sqrt(sum(square(i - j))));
}
distancetype sparseCosDist(const sp_mat& i, const sp_mat& j) {
  return 2.0 - 2.0 * (as_scalar((dot(i,j)) / as_scalar(norm(i,2) * norm(j,2))));
}
distancetype sparseRelDist(const sp_mat& i, const sp_mat& j) {
  return as_scalar(sum(square(i - j)));
}

/*
 * Fast calculation of pairwise distances with the result stored in a pre-allocated vector.
 */
// [[Rcpp::export]]
arma::vec fastDistance(const IntegerVector is,
                       const IntegerVector js,
                       const arma::mat& data,
                       const std::string& distMethod,
                       Rcpp::Nullable<Rcpp::NumericVector> threads,
                       bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  Progress p(is.size(), verbose);
  vec xs = vec(is.size());
  distancetype (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = dist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = cosDist;
#ifdef _OPENMP
#pragma omp parallel for shared (xs)
#endif
  for (vertexidxtype i=0; i < is.length(); i++) if (p.increment()) xs[i] =
    distanceFunction(data.col(is[i]), data.col(js[i]));
  return xs;
};

arma::vec fastSparseDistance(const arma::ivec& is,
                             const arma::ivec& js,
                             const sp_mat& data,
                             const std::string& distMethod,
                             bool verbose) {

  Progress p(is.size(), verbose);
  vec xs = vec(is.size());
  distancetype (*distanceFunction)(
      const sp_mat& x_i,
      const sp_mat& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = sparseDist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = sparseCosDist;
#ifdef _OPENMP
#pragma omp parallel for shared (xs)
#endif
  for (vertexidxtype i=0; i < is.n_elem; i++) if (p.increment()) xs[i] =
    distanceFunction(data.col(is[i]), data.col(js[i]));
  return xs;
};

// [[Rcpp::export]]
arma::vec fastCDistance(const arma::ivec& is,
                        const arma::ivec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& p_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        Rcpp::Nullable<Rcpp::NumericVector> threads,
                        bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  const vertexidxtype N = p_locations.n_elem - 1;
  const sp_mat data = sp_mat(i_locations, p_locations, x, N, N);
  return fastSparseDistance(is,js,data,distMethod,verbose);
}

// [[Rcpp::export]]
arma::vec fastSDistance(const arma::ivec& is,
                        const arma::ivec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& j_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        Rcpp::Nullable<Rcpp::NumericVector> threads,
                        bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  const umat locations = join_cols(i_locations, j_locations);
  const sp_mat data = sp_mat(locations, x);
  return fastSparseDistance(is,js,data,distMethod,verbose);
}
