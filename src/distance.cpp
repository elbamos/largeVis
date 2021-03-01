#include "largeVis.h"
#include "distance.h"
#include <progress.hpp>

distancetype relDist(const arma::vec& i, const arma::vec& j) {
  const dimidxtype D = i.n_elem;
  distancetype cnt = 0;
  for (dimidxtype idx = 0; idx < D; ++idx) cnt += ((i[idx] - j[idx]) * (i[idx] - j[idx]));
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
  for (dimidxtype d = 0; d < D; ++d) {
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

class DistanceWorker : public RcppParallel::Worker {
private:
	const arma::mat  data;
	const arma::ivec is;
	const arma::ivec js;
	distancetype (*distFunction)(const vec&, const vec&) = nullptr;
	Progress p;
	vec* out;

public:


	DistanceWorker(
		const arma::mat &data,
		const arma::ivec &is,
		const arma::ivec &js,
		distancetype (*distFunction)(const vec&, const vec&),
		vec *out,
		Progress &p) :
	data{data},
	is{is},
	js{js},
	distFunction{distFunction},
	out{out},
	p{p} {
	}

	void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t pos=begin; pos < end; ++pos) {
			(*out)[pos] = distFunction(data.col(is[pos]), data.col(js[pos]));
		}
	}
};

/*
 * Fast calculation of pairwise distances with the result stored in a pre-allocated vector.
 */
// [[Rcpp::export]]
arma::vec fastDistance(const arma::ivec& is,
                       const arma::ivec& js,
                       const arma::mat& data,
                       const std::string& distMethod,
                       bool verbose) {
	Progress p(is.size(), verbose);
	vec xs = vec(is.size());
	xs.fill(0);

	distancetype (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j);
	if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = dist;
	else distanceFunction = cosDist;

	DistanceWorker worker(data, is, js, distanceFunction, &xs, p);
	parallelFor(0, is.n_elem, worker);
	return xs;
}

class SparseDistanceWorker : public RcppParallel::Worker {
private:
	const arma::sp_mat  data;
	const arma::ivec is;
	const arma::ivec js;
	distancetype (*distFunction)(const sp_mat&, const sp_mat&) = nullptr;
	Progress p;
	vec* out;

public:


	SparseDistanceWorker(
		const arma::sp_mat &data,
		const arma::ivec &is,
		const arma::ivec &js,
		distancetype (*distFunction)(const sp_mat&, const sp_mat&),
		vec *out,
		Progress &p) :
	data{data},
	is{is},
	js{js},
	distFunction{distFunction},
	out{out},
	p{p} {
	}

	void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t pos=begin; pos < end; ++pos) {
			(*out)[pos] = distFunction(data.col(is[pos]), data.col(js[pos]));
		}
	}
};

vec fastSparseDistance(const ivec& is,
                       const ivec& js,
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

  SparseDistanceWorker worker(data, is, js, distanceFunction, &xs, p);
  parallelFor(0, is.n_elem, worker);
  return xs;
}

// [[Rcpp::export]]
arma::vec fastCDistance(const arma::ivec& is,
                        const arma::ivec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& p_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose) {
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
                        bool verbose) {
  const umat locations = join_cols(i_locations, j_locations);
  const sp_mat data = sp_mat(locations, x);
  return fastSparseDistance(is,js,data,distMethod,verbose);
}
