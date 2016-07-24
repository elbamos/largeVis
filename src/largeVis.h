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
                 const int max_recursion_degree,
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
  double* probs;
  T* aliases;
  T N;

public:
  AliasTable(const T n,
             const NumericVector& weights);
  T search(double *random) const;
  T search(double random, double random2) const;
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
  inline void multModify(double *col, double adj) const;
  // inline double max(double val) const;
  // inline double min(double val) const;
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
// class LookupGradient: public Gradient {
// private:
//   double Lookup(const double dist_squared, double* table) const;
// protected:
//   double bound;
//   double alpha;
//   int steps;
//   double boundsteps;
//   double* negativeLookup;
//   virtual void _positiveGradient(const double dist_squared,
//                                  double* holder) const;
//   virtual bool _negativeGradient(const double dist_squared,
//                                  double* holder) const;
// public:
//   LookupGradient(double alpha,
//                  double gamma,
//                  int d,
//                  double bound,
//                  int steps);
// };
class AlphaGradient: public Gradient {
  double alphagamma;
  double alpha;
  double twoalpha;
public:
  AlphaGradient(const double a,
                const double g,
                const int d);
protected:
  double cap;
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
