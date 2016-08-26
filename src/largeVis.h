#ifndef _LARGEVIS
#define _LARGEVIS
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "progress.hpp"
#include <math.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <set>
#include <memory>

using namespace Rcpp;
using namespace std;
using namespace arma;
/*
 * Global types
 */
typedef double distancetype;
typedef double coordinatetype;
#ifdef __MINGW32__
typedef long long int vertexidxtype;
#else
typedef long long vertexidxtype;
#endif
typedef long long edgeidxtype;
typedef long long iterationtype;
typedef int dimidxtype;
typedef int kidxtype;

#ifdef _OPENMP
void checkCRAN(Rcpp::Nullable<Rcpp::NumericVector> threads);
#endif

/*
 * Distance Functions
 */
distancetype dist(const arma::vec& i, const arma::vec& j);
distancetype relDist(const arma::vec& i, const arma::vec& j);
distancetype cosDist(const arma::vec& i, const arma::vec& j);
distancetype sparseDist(const sp_mat& i, const sp_mat& j);
distancetype sparseRelDist(const sp_mat& i, const sp_mat& j);
distancetype sparseCosDist(const sp_mat& i, const sp_mat& j);

// Exported distance functions for high dimensional space
arma::vec fastDistance(const NumericVector is,
                       const NumericVector js,
                       const arma::mat& data,
                       const std::string& distMethod,
                       bool verbose);
arma::vec fastSparseDistance(const arma::vec& is,
                             const arma::vec& js,
                             const sp_mat& data,
                             const std::string& distMethod,
                             bool verbose);
arma::vec fastCDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& p_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose);
arma::vec fastSDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& j_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose);

/*
 * Functions related to the alias algorithm
 */

template <class T>
class AliasTable {
private:
	std::unique_ptr< coordinatetype[] > probs;
	std::unique_ptr< T[] > aliases;
  std::uniform_real_distribution< coordinatetype > rnd = std::uniform_real_distribution< coordinatetype >();
  std::mt19937_64 mt;
  T N;

public:
  AliasTable() {
  }

	void initialize(const distancetype* weights, T N) {
		this -> N = N;
		probs = std::unique_ptr< coordinatetype[] >( new coordinatetype[N] );
		aliases = std::unique_ptr< T[] >(new T[N]);
		distancetype sm = 0;
		for (T i = 0; i != N; i++) sm += weights[i];
  	for (T i = 0; i != N; i++) probs[i] = weights[i] * N / sm;
  	queue<T> small = queue<T>();
  	queue<T> large = queue<T>();
  	for (T i = 0; i < N; i++) ((probs[i] < 1) ?
                                small :
                                large).push(i);
  	while (! large.empty() & ! small.empty()) {
  		T big = large.front();
  	  large.pop();
  		T little = small.front();
  		small.pop();
  		aliases[little] = big;
  		probs[big] = probs[big] + probs[little] - 1;
  		(probs[big] < 1 ? small : large).push(big);
  	}
  	long double accu = 0;
  	while (! large.empty()) {
  	  accu += 1 - large.front();
  		probs[large.front()] = 1;
  		large.pop();
  	}
  	while (! small.empty()) {
  	  accu += 1 - small.front();
  	  probs[small.front()] = 1;
  	  small.pop();
  	}
  	if (accu > 1e-5) warning("Numerical instability in alias table " + to_string(accu));
  };

  T search(coordinatetype random, coordinatetype random2) const {
    T candidate = random * N;
    return (random2 >= probs[candidate]) ? aliases[candidate] : candidate;
  };

  long long initRandom(long long seed) {
  	mt = mt19937_64(seed);
  	return mt();
  }
  void initRandom() {
  	std::random_device seed;
  	initRandom(seed());
  }
  coordinatetype getRand() {
  	return rnd(mt);
  }
  T operator()() {
  	return search(rnd(mt), rnd(mt));
  }
};

/*
 * Gradients
 */
class Gradient {
protected:
  const distancetype gamma;
  distancetype cap;
  const dimidxtype D;
  Gradient(const distancetype g,
           const dimidxtype d) : gamma{g}, cap(5), D{d} {};
  virtual void _positiveGradient(const distancetype dist_squared,
                                coordinatetype* holder) const = 0;
  virtual void _negativeGradient(const distancetype dist_squared,
                                 coordinatetype* holder) const = 0;
  inline void multModify(coordinatetype *col, coordinatetype adj) const;
  inline coordinatetype clamp(coordinatetype val) const;

public:
  virtual void positiveGradient(const coordinatetype* i,
                                const coordinatetype* j,
                                coordinatetype* holder) const;
  virtual void negativeGradient(const coordinatetype* i,
                                const coordinatetype* k,
                                coordinatetype* holder) const;
  inline distancetype distAndVector(const coordinatetype *x_i,
			                        const coordinatetype *x_j,
			                        coordinatetype *output) const;
};

class AlphaGradient: public Gradient {
  const coordinatetype alpha;
  const coordinatetype twoalpha;
protected:
  const coordinatetype alphagamma;
  virtual void _positiveGradient(const double dist_squared,
                                 coordinatetype* holder) const;
  virtual void _negativeGradient(const double dist_squared,
                                 coordinatetype* holder) const;
public:
  AlphaGradient(const distancetype a,
                const distancetype g,
                const dimidxtype D) : Gradient(g, D),
                               alpha{a},
                               twoalpha(alpha * -2),
                               alphagamma(alpha * gamma * 2) { } ;
};

class AlphaOneGradient: public AlphaGradient {
public:
  AlphaOneGradient(const distancetype g,
                   const dimidxtype D);
protected:
  virtual void _positiveGradient(const distancetype dist_squared,
                                 coordinatetype* holder) const;
  virtual void _negativeGradient(const distancetype dist_squared,
                                 coordinatetype* holder) const;
};

class ExpGradient: public Gradient {
public:
  const coordinatetype gammagamma;
  ExpGradient(const distancetype g, const dimidxtype d) : Gradient(g, d),
                                          	 gammagamma(gamma * gamma) {
    cap = gamma;
  };
protected:
  virtual void _positiveGradient(const distancetype dist_squared,
                                 coordinatetype* holder) const;
  virtual void _negativeGradient(const distancetype dist_squared,
                                 coordinatetype* holder) const;
};
#endif
