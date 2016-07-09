// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

// utility function
inline void multModify(double *col, int D, double adj) {
  for (int i = 0; i != D; i++) col[i] *= adj;
}

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

// Parent class
Gradient::Gradient(const double g, 
                   const int d) {
  gamma = g;
  D = d;
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
  return _negativeGradient(dist_squared, holder);
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
  positiveLookup = new double[steps];
  negativeLookup = new double[steps];
  double step = bound / steps;
  double distsquared = 0;
  for (int i = 1; i != steps; i++) {
    distsquared +=step;
    const double adk = alpha * distsquared;
    negativeLookup[i] = alpha * gamma / (sqrt(distsquared) * (adk + 1));
  //  positiveLookup[i] = - alpha / (1 + alpha * distsquared);
  }
  negativeLookup[0] = negativeLookup[1];
  //positiveLookup[0] = positiveLookup[1];
}

inline int fastround(double r) {
  return (r > 0.0) ? (r + 0.5) : (r - 0.5);
}

double LookupGradient::Lookup(const double dist_squared, 
                            double* table) const {
  double grad;
  if (dist_squared > bound) grad = table[steps - 1];
  else grad = table[fastround(dist_squared * boundsteps)];
  return grad;
}

void LookupGradient::_positiveGradient(const double dist_squared, 
                               double* holder) const {
  // double grad = Lookup(dist_squared, positiveLookup);
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
  cap = 2;
}
void AlphaGradient::_positiveGradient(const double dist_squared, 
                                      double* holder) const {
  const double grad = - alpha / (1 + alpha * dist_squared);
  multModify(holder, D, grad);
}
bool AlphaGradient::_negativeGradient(const double dist_squared, 
                                      double* holder) const {
  if (dist_squared == 0) return true; // If the two points are in the same place, skip
  const double dist = sqrt(dist_squared);
  const double adk = alpha * dist_squared;
  double grad = alphagamma / (dist * (adk + 1));
 // if (grad > cap) grad = cap;
  // multModify(holder, D, (grad > gamma / 4) ? gamma / 4 : grad);
  // return false;
  
  multModify(holder, D, grad);
  for (int d = 0; d != D; d++) {
    if (holder[d] > 2 * alphagamma) stop(std::to_string(holder[d]) + " " + 
        std::to_string(grad) + " " + std::to_string(dist_squared) + " " + 
        "negative > 2");
    if (holder[d] < -2 * alphagamma) stop(std::to_string(holder[d]) + " " + 
        std::to_string(grad) + " " + std::to_string(dist_squared) + " " + 
        "negative < -2");
    // if (holder[d] > cap) holder[d] = cap;
    // if (holder[d] < - cap) holder[d] = - cap;
  }
  return false;
}

/*
 * Optimized gradient for alpha == 1
 */
AlphaOneGradient::AlphaOneGradient(const double g, 
                                   const int d) : AlphaGradient(1, g, d) {
  cap = g;
}
void AlphaOneGradient::_positiveGradient(const double dist_squared,  
                                         double* holder) const {
  const double grad = - 1 / (1 + dist_squared);
  multModify(holder, D, grad);
}
bool AlphaOneGradient::_negativeGradient(const double dist_squared, 
                                         double* holder) const {
  if (dist_squared == 0) return true; // If the two points are in the same place, skip
  double grad = gamma / (sqrt(dist_squared) * (dist_squared + 1));
  multModify(holder, D, grad);
  for (int d = 0; d != D; d++) {
    if (holder[d] > 2 * gamma) stop(std::to_string(holder[d]) + " " + 
        std::to_string(grad) + " " + std::to_string(dist_squared) + " " + 
        "negative1 > 2");
    if (holder[d] < -2 * gamma) stop(std::to_string(holder[d]) + " " + 
        std::to_string(grad) + " " + std::to_string(dist_squared) + " " + 
        "negative` < -2");
    // if (holder[d] > cap) holder[d] = cap;
    // if (holder[d] < - cap) holder[d] = - cap;
  }
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
  else grad = new AlphaGradient(a, 1, i.size() );
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
  else grad = new AlphaGradient(a, g, i.size() );
  grad -> negativeGradient(i.memptr(), j.memptr(), holder.memptr());
  return holder;
};
