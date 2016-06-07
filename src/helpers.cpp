#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
using namespace Rcpp;

double relDist(const arma::vec& i, const arma::vec& j) {
  const int lim = i.n_elem;
  double cnt = 0;
  for (int idx = 0; idx < lim; idx++) cnt += ((i[idx] - j[idx]) * (i[idx] - j[idx]));
  return cnt;
}

double dist(const arma::vec& i, const arma::vec& j) {
  return sqrt(relDist(i,j));
}

// Stolen directly from Erik B's annoylib.h
double cosDist(const arma::vec& i, const arma::vec& j) {
  int lim = i.n_elem;
  double pp = 0, qq = 0, pq = 0;
  for (int z = 0; z < lim; z++) {
    pp += (i[z]) * (i[z]);
    qq += (j[z]) * (j[z]);
    pq += (i[z]) * (j[z]);
  }
  double ppqq = pp * qq;
  if (ppqq > 0) return 2.0 - 2.0 * pq / sqrt(ppqq);
  else return 2.0; // cos is 0
}

double sparseDist(const arma::sp_mat& i, const arma::sp_mat& j) {
  return arma::as_scalar(sqrt(sum(square(i - j))));
}

double sparseCosDist(const arma::sp_mat& i, const arma::sp_mat& j) {
  return 1 - (arma::as_scalar((dot(i,j)) / arma::as_scalar(arma::norm(i,2) * arma::norm(j,2))));
}

/*
 * Function for reading file in the LINE format used by the paper authors for two of their examples.
 */
// [[Rcpp::export]]
NumericMatrix readLINE(std::string path) {
  FILE * infile = fopen(path.c_str(), "r");
  int N, D;
  if (fscanf(infile, "%d %d", &N, &D) != 2) stop("File error");
  NumericMatrix val = NumericMatrix(N, D);
  std::vector< std::string > nodeNames;
  double * buffer = (double *) malloc(sizeof(double) * D);
  for (int n = 0; n < N; n++) {
    char nodeName [80];
    if (fscanf(infile, "\n%s ", &nodeName) != 1) stop("Name error.");
    nodeNames.push_back(nodeName);
    // double r;
    // for (int d = 0; d < D; d++) {
    //   int tim = fscanf(infile, "%lf ", &r);
    //   Rcout << tim << " ";
    //   val(n,d) = r;
    // }
    // float f;
    // for (int d = 0; d < D; d++) {
    //   if (fread(&f, sizeof(float), 1, infile) != 1) stop("Read error.");
    //   val(n,d) = static_cast<double>(f);
    // }
    if (fread((void *)buffer, sizeof(double), D, infile) != D) stop("Read error.");
    for (int d = 0; d < D; d++) if (buffer[d] < 1e36) val(n,d) = buffer[d];
    else stop("Bad range");
  }
  rownames(val) = Rcpp::wrap(nodeNames);
  return val;
}
