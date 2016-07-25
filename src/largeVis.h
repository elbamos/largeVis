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
#define _GSLRANDOM
#ifdef _GSLRANDOM
#include <gsl/gsl_rng.h>
#endif
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

template <class T>
class AliasTable {
private:
  long double* probs;
  T* aliases;
  T N;

public:
  AliasTable(const T n,
             const NumericVector& weights) {
  	N = n;
  	// norm_probs = new double[n];
  	probs = new long double[n];
  	aliases = new T[n];
  	
  	const double sm = sum(weights); //accu(weights);
  	for (T i = 0; i < n; i++) probs[i] = weights[i] * n / sm;
  	queue<T> small = queue<T>();
  	queue<T> large = queue<T>();
  	for (T i = 0; i < n; i++) ((probs[i] < 1) ?
                                small :
                                large).push(i);
  	while (! large.empty() & ! small.empty()) {
  		T big = large.front();
  		large.pop();
  		T little = small.front();
  		small.pop();
  		// probs[little] = norm_probs[little];
  		aliases[little] = big;
  		probs[big] = probs[big] + probs[little] - 1;
  		(probs[big] < 1 ? small : large).push(big);
  		// norm_probs[big] = norm_probs[big] + norm_probs[little] - 1;
  		// ((norm_probs[big] < 1) ? small : large).push(big);
  	}
  	while (! large.empty()) {
  		probs[large.front()] = 1;
  		large.pop();
  	}
  	while (! small.empty()) {
  	  warning("Numeric instability in alias table.");
  	  probs[small.front()] = 1;
  	  small.pop();
  	}
  };

  T search(double *random) const {
    return search(random[0], random[1]);
  };
  
  T search(double random, double random2) const {
    T candidate = random * N;
    return (random2 >= probs[candidate]) ? aliases[candidate] : candidate;
  };
  
#ifdef _GSLRANDOM
  const gsl_rng_type *gsl_T = NULL;
  gsl_rng *gsl_r = NULL;
  
  void initRandom() {
    initRandom(314159265);
  }
  void initRandom(long seed) {
    gsl_T = gsl_rng_rand48;
    gsl_r = gsl_rng_alloc(gsl_T);
    gsl_rng_set(gsl_r, seed);
  }
  T sample() const {
    return search(gsl_rng_uniform(gsl_r), gsl_rng_uniform(gsl_r));
  }
#endif
  

};

/*
 * Gradients
 */
class Gradient {
protected:
  double gamma;
  double cap;
  int D;
  Gradient(const double g,
           const int d);
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
  double alpha;
  double twoalpha;
public:
  AlphaGradient(const double a,
                const double g,
                const int d);
protected:
  double alphagamma;
  virtual void _positiveGradient(const double dist_squared,
                        double* holder) const;
  virtual void _negativeGradient(const double dist_squared,
                        double* holder) const;
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
  double gammagamma;
  ExpGradient(const double g,
              const int d);
protected:
  virtual void _positiveGradient(const double dist_squared,
                        double* holder) const;
  virtual void _negativeGradient(const double dist_squared,
                        double* holder) const;
};

arma::vec testNegativeGradient(arma::vec i, arma::vec j,
                               NumericVector alpha, NumericVector gamma, NumericVector f);
arma::vec testPositiveGradient(arma::vec i, arma::vec j,
                               NumericVector alpha, NumericVector f);

/*
 * The SGD loop
 */
arma::mat sgd(arma::mat coords,
              arma::ivec& is,
              const IntegerVector js,
              const IntegerVector ps,
              const arma::vec ws,
              const double gamma,
              const double rho,
              const double minRho,
              const bool useWeights,
              const long nBatches,
              const int M,
              const double alpha,
              bool verbose);
#endif
