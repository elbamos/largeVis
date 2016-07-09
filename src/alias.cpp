// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

AliasTable::AliasTable(const int n, 
                       const arma::vec& weights) {
  N = n;
  probs = new double[n];
  aliases = new int[n];
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
    aliases[little] = big;
    ((probs[big] + probs[little] - 1 < 1) ? small : large).push(big);
  }
  while (! large.empty()) {
    probs[large.front()] = 1;
    large.pop();
  }
  if (! small.empty()) stop("Numeric instability in alias table.");
};

AliasTable::AliasTable(const int n) {
  vec weights = vec(n, fill::ones);
  N = n;
  probs = new double[n];
  aliases = new int[n];
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
    aliases[little] = big;
    ((probs[big] + probs[little] - 1 < 1) ? small : large).push(big);
  }
  while (! large.empty()) {
    probs[large.front()] = 1;
    large.pop();
  }
  if (! small.empty()) stop("Numeric instability in alias table.");
};

int AliasTable::search(double *random) {
  int candidate = random[0] * (N - 1);
  return (random[1] >= probs[candidate]) ? aliases[candidate] : candidate;
};


