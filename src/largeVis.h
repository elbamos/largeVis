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

using namespace Rcpp;
using namespace std;
using namespace arma;

/*
 * Neighbor search
 */
struct HeapObject {
  double d;
  int n;
  HeapObject(double d, int n) : d(d), n(n) {}
  bool operator<(const struct HeapObject& other) const {
    return d < other.d;
  }
};
typedef priority_queue<HeapObject> MaxHeap;
typedef vector< imat::col_iterator > PositionVector;
typedef vector<int> Neighborhood;
Neighborhood** createNeighborhood(int N);
void copyHeapToMatrix(set<int>* tree,
                      const int K,
                      const int i,
                      arma::imat& knns);
void addDistance(const arma::vec& x_i,
                 const arma::mat& data,
                 const int j,
                 MaxHeap& heap,
                 const int K,
                 double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j));
void heapToSet(MaxHeap& heap, set<int>* set);
arma::imat annoy(const int n_trees,
                 const int threshold,
                 const arma::mat& data,
                 const int N,
                 const int K,
                 double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j),
                 Progress& p);
void addNeighbors(const arma::ivec& indices,
                  Neighborhood* heap[],
                  const int I);
/*
 * Distance Functions
 */
double dist(const arma::vec& i, const arma::vec& j);
double relDist(const arma::vec& i, const arma::vec& j);
double cosDist(const arma::vec& i, const arma::vec& j);
double sparseDist(const sp_mat& i, const sp_mat& j);
double sparseRelDist(const sp_mat& i, const sp_mat& j);
double sparseCosDist(const sp_mat& i, const sp_mat& j);

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
typedef double realsies;

template <class T>
class AliasTable {
private:
	realsies* probs;
  T* aliases;
  T N;

public:
  AliasTable(const NumericVector& weights) : N(weights.length()) {
    aliases = new T[N];
    probs = new realsies[N];

  	const long double sm = sum(weights);
  	for (T i = 0; i < N; i++) probs[i] = weights[i] * N / sm;
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

  T search(realsies *random) const {
    return search(random[0], random[1]);
  };

  T search(realsies random, realsies random2) const {
    T candidate = random * N;
    return (random2 >= probs[candidate]) ? aliases[candidate] : candidate;
  };

  // const gsl_rng_type *gsl_T = NULL;
  // gsl_rng *gsl_r = NULL;
  //
  // void initRandom() {
  //   initRandom(314159265);
  // }
  // void initRandom(long seed) {
  //   gsl_T = gsl_rng_rand48;
  //   gsl_r = gsl_rng_alloc(gsl_T);
  //   gsl_rng_set(gsl_r, seed);
  // }
  std::uniform_real_distribution<realsies> rnd;
  std::mt19937_64 mt;

  void initRandom(long seed) {
  	mt = mt19937_64(seed);
  	rnd = uniform_real_distribution<realsies>();
  }
  void initRandom() {
  	std::random_device seed;
  	initRandom(seed());
  }

  T sample() {
  	realsies dub1 = rnd(mt);
  	realsies dub2 = rnd(mt);
  	return search(dub1, dub2);
    // return search(gsl_rng_uniform(gsl_r), gsl_rng_uniform(gsl_r));
  }
};

/*
 * Gradients
 */
class Gradient {
protected:
  const double gamma;
  double cap;
  const int D;
  Gradient(const double g,
           const int d) : gamma{g}, cap(5), D{d} {};
  virtual void _positiveGradient(const double dist_squared,
                                double* holder) const = 0;
  virtual void _negativeGradient(const double dist_squared,
                                double* holder) const = 0;
  inline void multModify(double *col, const double adj) const;
  inline double clamp(double val) const;

public:
  virtual void positiveGradient(const double* i,
                                const double* j,
                                double* holder) const;
  virtual void negativeGradient(const double* i,
                                const double* k,
                                double* holder) const;
  inline double distAndVector(const double *x_i,
                       const double *x_j,
                       double *output) const;
};

class AlphaGradient: public Gradient {
  const double alpha;
  const double twoalpha;
protected:
  const double alphagamma;
  virtual void _positiveGradient(const double dist_squared,
                                 double* holder) const;
  virtual void _negativeGradient(const double dist_squared,
                                 double* holder) const;
public:
  AlphaGradient(const double a,
                const double g,
                const int d) : Gradient(g, d),
                               alpha{a},
                               twoalpha(alpha * -2),
                               alphagamma(alpha * gamma * 2) { } ;
};

class AlphaOneGradient: public AlphaGradient {
public:
  AlphaOneGradient(const double g,
                   const int d);
protected:
  virtual void _positiveGradient(const double dist_squared,
                                 double* holder) const;
  virtual void _negativeGradient(const double dist_squared,
                                 double* holder) const;
};

class ExpGradient: public Gradient {
public:
  const double gammagamma;
  ExpGradient(const double g, const int d) : Gradient(g, d),
                                          	 gammagamma(gamma * gamma) {
    cap = gamma;
  };
protected:
  virtual void _positiveGradient(const double dist_squared,
                        				 double* holder) const;
  virtual void _negativeGradient(const double dist_squared,
                        				 double* holder) const;
};
#endif
