// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]   
#include "largeVis.h"

// Copies the vector sums into a vector while it computes distance^2 -
// useful in calculating the gradients during SGD
inline double Gradient::distAndVector(const double *x_i,
                                      const double *x_j,
                                      double *output) const {
  double cnt = 0;
  for (int d = 0; d < D; d++) {
    double t = x_i[d] - x_j[d];
    output[d] = t;
    cnt += t * t;
  }
  return cnt;
}

inline void Gradient::multModify(double *col, int D, double adj) const {
  for (int i = 0; i != D; i++) {
    col[i] *= adj;
    if (col[i] > cap) col[i] = cap;
    else if (col[i] < - cap) col[i] = -cap;
  }
}

// Parent class
Gradient::Gradient(const double g, 
                   const int d) {
  gamma = g;
  D = d;
  cap = 5;
}
void Gradient::positiveGradient(const double* i, 
                                const double* j, 
                                double* holder) const { 
  const double dist_squared = distAndVector(i, j, holder);
  _positiveGradient(dist_squared, holder);
}
bool Gradient::negativeGradient(const double* i, 
                                const double* k,
                                double* holder) const {
  const double dist_squared = distAndVector(i, k, holder);
  bool val = _negativeGradient(dist_squared, holder);
  return val;
}

/*
 * Class for a gradient found using a pre-calculated lookup table
 */
LookupGradient::LookupGradient(double alpha, 
                               double gamma, 
                               int d,
                               double bound, 
                               int steps) : Gradient(gamma, d) {
  this -> alpha = alpha;
  this -> bound = bound;
  this -> steps = steps;
  this -> boundsteps = (steps - 1) / bound;
  negativeLookup = new double[steps];
  double step = bound / steps;
  double distsquared = 0;
  for (int i = 1; i != steps; i++) {
    distsquared +=step;
    const double adk = alpha * distsquared;
    negativeLookup[i] = alpha * gamma / (distsquared * (adk + 1));
  }
  negativeLookup[0] = negativeLookup[1];
}

inline int fastround(double r) {
  return r + 0.5;
}

double LookupGradient::Lookup(const double dist_squared, 
                            double* table) const {
  return (dist_squared > bound) ? table[steps - 1] :
                            table[fastround(dist_squared * boundsteps)];
}

void LookupGradient::_positiveGradient(const double dist_squared, 
                               double* holder) const {
  const double grad = - alpha / (1 + alpha * dist_squared);
  multModify(holder, D, grad);
}
bool LookupGradient::_negativeGradient(const double dist_squared, 
                               double* holder) const {
  if (dist_squared == 0) return true;
  double grad = Lookup(dist_squared, negativeLookup);
  multModify(holder, D, grad);
  return false;
}

/*
 * Generalized gradient with an alpha parameter
 */

AlphaGradient::AlphaGradient(const double a, 
                             const double g, 
                             const int d) : Gradient(g, d) {
  alpha = a;
  alphagamma = a * g;
}
void AlphaGradient::_positiveGradient(const double dist_squared, 
                                      double* holder) const {
  const double grad = - alpha / (1 + alpha * dist_squared);
  multModify(holder, D, grad);
}
bool AlphaGradient::_negativeGradient(const double dist_squared, 
                                      double* holder) const {
  if (dist_squared == 0) return true; // If the two points are in the same place, skip
  const double adk = alpha * dist_squared;
  double grad = alphagamma / (dist_squared * (adk + 1));
  multModify(holder, D, grad);
  return false;
}

/*
 * Optimized gradient for alpha == 1
 */
AlphaOneGradient::AlphaOneGradient(const double g, 
                                   const int d) : AlphaGradient(1, g, d) {
  this -> gamma = 2 * g;
}
void AlphaOneGradient::_positiveGradient(const double dist_squared,  
                                         double* holder) const {
  const double grad = - 2 / (1 + dist_squared);
  multModify(holder, D, grad);
}
bool AlphaOneGradient::_negativeGradient(const double dist_squared, 
                                         double* holder) const {
  if (dist_squared == 0) return true; // If the two points are in the same place, skip
  double grad = gamma / (1 + dist_squared) / (0.1 + dist_squared); //(dist_squared * (dist_squared + 1));
  multModify(holder, D, grad);
  return false;
}

/*
 * Alternative probabilistic function (sigmoid)
 */

ExpGradient::ExpGradient(const double g,
                         const int d) : Gradient(g, d) {
  gammagamma = g * g;
}
void ExpGradient::_positiveGradient(const double dist_squared, 
                                   double* holder) const {
  const double expsq = exp(dist_squared);
  const double grad = (dist_squared > 4) ? -1 : 
                                           -(expsq / (expsq + 1));
  multModify(holder, D, grad);
}
bool ExpGradient::_negativeGradient(const double dist_squared, 
                                   double* holder) const {
  if (dist_squared == 0) return true; 
  const double grad = (dist_squared > gammagamma) ? 0 : 
                                                    gamma / (1 + exp(dist_squared));
  multModify(holder, D, grad);
  return false;
}

// [[Rcpp::export]]
arma::vec testPositiveGradient(arma::vec i, arma::vec j,
                               NumericVector alpha, NumericVector f) {
  double a = alpha[0];
  vec holder = vec(i.size());
  Gradient* grad;
  if (a == 0) grad = new ExpGradient(1, i.size());
  else grad = new LookupGradient(a,
                            1,
                            i.size(),
                            7, //bound
                            1000);
  grad -> positiveGradient(i.memptr(), j.memptr(), holder.memptr());
  return holder;
};
// [[Rcpp::export]]
arma::vec testNegativeGradient(arma::vec i, arma::vec j,
                               NumericVector alpha, NumericVector gamma, NumericVector f) {
  double a = alpha[0];
  double g = gamma[0];
  vec holder = vec(i.size());
  Gradient* grad;
  if (a == 0) grad = new ExpGradient(g, i.size());
  else grad = new LookupGradient(a,
                                 g,
                                 i.size(),
                                 49, //bound
                                 1000);  grad -> negativeGradient(i.memptr(), j.memptr(), holder.memptr());
  return holder;
};
