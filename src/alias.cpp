// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

/*
 * Functions for implementing the alias algorithm to get weighted samples
 * in O(1)
 */

void makeAliasTable(int n, arma::vec weights, double *probs, int *alias) {
  const double sm = accu(weights);
  for (int i = 0; i < n; i++) probs[i] = weights[i] * n / sm;
  queue<int> small = queue<int>();
  queue<int> large = queue<int>();
  for (int i = 0; i < n; i++) ((probs[i] < 1) ?
                                 small :
                                 large).push(i);
  while (! large.empty() & ! small.empty()) {
    int big = large.front();
    large.pop();
    int little = small.front();
    small.pop();
    alias[little] = big;
    ((probs[big] + probs[little] - 1 < 1) ? small : large).push(big);
  }
  while (! large.empty()) {
    probs[large.front()] = 1;
    large.pop();
  }
  if (! small.empty()) stop("Numeric instability in alias table.");
};
int searchAliasTable(double *random, // Two random numbers
                     double *probs, // Probs table for alias
                     int *alias, // Alias table
                     const int& N // Total number of entries in alias table
                       ) {
  int candidate = random[0] * (N - 1);
  return (random[1] >= probs[candidate]) ? alias[candidate] : candidate;
};
