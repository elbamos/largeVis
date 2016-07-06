#include <testthat.h>
#include "largeVis.h"

context("distances") {

  test_that("distAndVector sets the vector") {
    arma::vec left = arma::randu<arma::vec>(5);
    arma::vec right = arma::randu<arma::vec>(5);
    double tmpOut[5];
    double d = distAndVector(left.memptr(), right.memptr(), tmpOut, left.size());
    for (int d = 0; d < 5; d++) expect_true(tmpOut[d] == left[d] - right[d]);
  }
  test_that("distAndVector does the squared distance") {
    arma::vec left = arma::randu<arma::vec>(5);
    arma::vec right = arma::randu<arma::vec>(5);
    double tmpOut[5];
    double d = distAndVector(left.memptr(), right.memptr(), tmpOut, left.size());
    double comp = relDist(left, right);
    expect_true(comp == d);
  }
}

void positiveTrigGradient(double* i, double* j,
                          double* holder,
                          const double alpha, const double f,
                          const int D) {
  const double dist_ij = sqrt(distAndVector(j, i, holder, D));
  const double grad = (-1) * ((alpha == 0) ?
                                exp(dist_ij * dist_ij) * dist_ij / (exp(dist_ij * dist_ij) + 1) :
                                2 * alpha * dist_ij / ((alpha * dist_ij * dist_ij) + 1));
  const double angle = atan2(holder[1],  holder[0]);
  holder[1] = - sin(angle) * grad;
  holder[0] = - cos(angle) * grad;
};

bool negativeTrigGradient(double* i,
                          double* k,
                          double* holder,
                          const double alpha, const double gamma, const double f,
                          const int D) {
  const double dist_ik = sqrt(distAndVector(k, i, holder, D));
  if (dist_ik == 0) return true; // If the two points are in the same place, skip
  double grad = gamma * ((alpha == 0) ?
                                     dist_ik / (exp(dist_ik * dist_ik) + 1) :
                                     1 / ( alpha * dist_ik * dist_ik * dist_ik + dist_ik));
  const double angle = atan2(holder[1], holder[0]);
  if (grad > 1) grad = 1;
  holder[1] = - sin(angle) * grad;
  holder[0] = - cos(angle) * grad;
  return false;
};

void checkMatch(double *a, double *b, double *ina, double *inb) {
  if (a[0] != b[0] || a[1] != b[1]) {
    Rcout << "\n" << ina[0] << " " << ina[1] << "\t" << inb[0] << " " << inb[1] <<
      "\t" << a[0] << " " << a[1] << "\t" << b[0] << " " << b[1];
  }
}

context("Gradients") {
  test_that("Positive matches trigonometric, alpha == 0") {
    arma::vec rand = arma::randu<arma::vec>(4);
    double holder[2];
    double holder2[2];
    double alpha = 0;
    double f = 2;
    int D = 2;
    positiveGradient(rand.memptr(), rand.memptr() + D, holder, alpha, D);
    positiveTrigGradient(rand.memptr(), rand.memptr() + D, holder2, alpha, f, D);
    checkMatch(holder, holder2, rand.memptr(), rand.memptr() + D);
    for (int d = 0; d < D; d++) expect_true(pow(holder[d] - holder2[d], 2) < 1e-4);
  }
  test_that("Negative matches trigonometric, alpha == 0") {
    arma::vec rand = arma::vec(4);
    rand[0] = -1; rand[1] = -1; rand[2] = 1; rand[3] = 1;
    double holder[5];
    double holder2[5];
    double alpha = 0;
    double f = 2;
    double gamma = 7;
    int D = 2;
    negativeGradient(rand.memptr(), rand.memptr() + D, holder, 
                     alpha, gamma, gamma / 4, D);
    negativeTrigGradient(rand.memptr(), rand.memptr() + D, holder2, alpha, gamma, f, D);
    checkMatch(holder, holder2, rand.memptr(), rand.memptr() + D);
    for (int d = 0; d < D; d++) expect_true(pow(holder[d] - holder2[d], 2) < 1e-4);
  }
  test_that("Positive matches trigonometric, alpha == 1") {
    arma::vec rand = arma::randu<arma::vec>(4);
    double holder[2];
    double holder2[2];
    double alpha = 1;
    double f = 2;
    int D = 2;
    positiveGradient(rand.memptr(), rand.memptr() + D, holder, alpha, D);
    positiveTrigGradient(rand.memptr(), rand.memptr() + D, holder2, alpha, f, D);
    checkMatch(holder, holder2, rand.memptr(), rand.memptr() + D);
    for (int d = 0; d < D; d++) expect_true(pow(holder[d] - holder2[d], 2) < 1e-4);
  }
  test_that("Negative matches trigonometric, alpha == 1") {
    arma::vec rand = arma::vec(4);
    rand[0] = -1; rand[1] = -1; rand[2] = 1; rand[3] = 1;
    double holder[5];
    double holder2[5];
    double alpha = 1;
    double f = 2;
    double gamma = 7;
    int D = 2;
    negativeGradient(rand.memptr(), rand.memptr() + D, holder, 
                     alpha, gamma, gamma / 4, D);
    negativeTrigGradient(rand.memptr(), rand.memptr() + D, holder2, alpha, gamma, f, D);
    checkMatch(holder, holder2, rand.memptr(), rand.memptr() + D);
    for (int d = 0; d < D; d++) expect_true(pow(holder[d] - holder2[d], 2) < 1e-3);
  }
}
