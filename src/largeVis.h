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

struct heapObject {
  /*
   * Unite node index and distance in a single object, facilitates
   * creation of the max heap
   */
  double d;
  int n;
  heapObject(double d, int n) : d(d), n(n) {}
  bool operator<(const struct heapObject& other) const {
    return d < other.d;
  }
};
typedef priority_queue<heapObject> maxHeap;
typedef vector< imat::col_iterator > positionVector;
typedef vector<int> neighborhood;

/*
 * Distance Functions
 */
double relDist(const arma::vec& i, const arma::vec& j);
double sparseDist(const sp_mat& i, const sp_mat& j);
double dist(const arma::vec& i, const arma::vec& j);
double distAndVector(double *x_i,
                     double *x_j,
                     double *output,
                     const int& d);
double cosDist(const arma::vec& i, const arma::vec& j);
double sparseCosDist(const sp_mat& i, const sp_mat& j);
double sparseRelDist(const sp_mat& i, const sp_mat& j);
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
void makeAliasTable(int n, arma::vec weights, double *probs, int *alias);
int searchAliasTable(const int& n,
                     double *random,
                     double *probs,
                     int *alias,
                     const int& sample_n);

/*
 * Gradient Functions
 *
 *  *Gradient -- should be the correct analytical gradients
 *  test*Gradient -- exports results of positive and negativeGradient to R for comparison with
 *                    autodifferentiation
 */
void positiveGradient(double* i, double* j,
                      double* holder,
                      const double alpha,
                      const int D);
bool negativeGradient(double* i,
                      double* k,
                      double* holder,
                      const double alpha, const double gamma,
                      const int D);
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
