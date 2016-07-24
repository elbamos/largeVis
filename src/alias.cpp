// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

template <class T>
AliasTable<T>::AliasTable(const T n,
                       		const NumericVector& weights) {
  N = n;
  probs = new double[n];
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
    aliases[little] = big;
    ((probs[big] + probs[little] - 1 < 1) ? small : large).push(big);
  }
  while (! large.empty()) {
    probs[large.front()] = 1;
    large.pop();
  }
  if (! small.empty()) stop("Numeric instability in alias table.");
};

template <class T>
T AliasTable<T>::search(double *random) const {
  return search(random[0], random[1]);
};

template <class T>
T AliasTable<T>::search(double random, double random2) const {
  T candidate = random * (N - 1);
  return (random2 >= probs[candidate]) ? aliases[candidate] : candidate;
};

