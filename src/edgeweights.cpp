// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

// Take four vectors (i indices, j indices, edge distances, and sigmas), and calculate
// p(j|i) and then w_{ij}.
// [[Rcpp::export]]
arma::sp_mat distMatrixTowij( const NumericVector sources,
                              const NumericVector targets,
                              const NumericVector weights,
                              const NumericVector sigmas,
                              const int N,
                              bool verbose) {
  Progress p(weights.size() * 2, verbose);
  vec rowSums = vec(N);
  vec pjis = vec(sources.length());
  for (int idx=0; idx != N; idx++) rowSums[idx] = 0;
  // Compute pji, accumulate rowSums at the same time
#ifdef _OPENMP
#pragma omp parallel for shared(pjis, rowSums)
#endif
  for (int e=0; e < pjis.size(); e++) if (p.increment()){
    const int i = sources[e];
    const double pji = exp(- pow(weights[e], 2)) / sigmas[i];
    pjis[e] = pji;
#ifdef _OPENMP
#pragma omp atomic
#endif
    rowSums[i] += pji;
  }
  if (p.check_abort()) return sp_mat(0);
  // Now convert p(j|i) to w_{ij} by symmetrizing.
  // Loop through the edges, and populate a location matrix and value vector for
  // the sp_mat batch insertion constructor.  Put all coordinates in the
  // lower triangle.  The constructor will automatically add duplicates.
  vec values = vec(pjis.size());
  umat locations = umat(2, pjis.size());
#ifdef _OPENMP
#pragma omp parallel for shared(locations, values)
#endif
  for (int e = 0; e < pjis.size(); e++) if (p.increment()) {
    int newi = sources[e], newj = targets[e];
    if (newi < newj) swap(newi, newj);
    values[e] =  ((pjis[e] / rowSums[sources[e]]) / (2 * N));
    locations(1,e) = newi;
    locations(0,e) = newj;
  }
  sp_mat wij = sp_mat(
    true, // add_values
    locations,
    values,
    N, N // n_col and n_row
  );
  wij = wij + wij.t();
  return wij;
};


// [[Rcpp::export]]
double sigFunc(const double& twosigmasquared,
               const NumericVector& x_i,
               const double& perplexity) {
  const NumericVector xs = exp(- x_i / twosigmasquared);
  const double sum_weight = sum(xs);
  double H = sum(x_i * xs) / twosigmasquared;
  H = (H / sum_weight) + log(sum_weight);
  return pow(perplexity - H, 2);
};
// double sigFunc(const double& sigma,
//                const NumericVector& x_i,
//                const double& perplexity) {
//   const NumericVector xs = exp(- x_i * x_i / sigma);
//   const NumericVector softxs = xs / sum(xs);
//   const double p2 = - sum(log(softxs) / log(2)) / xs.length();
//   return pow(perplexity - p2, 2);
// };



class ReferenceEdges {
protected:
  arma::vec sigmas;
  int n_vertices;
  int n_edges;
  std::vector<int> edge_from, edge_to, head, next, reverse;
  std::vector<double> edge_weight;
  double perplexity;
public:
  void similarityOne(int id) {
    double beta, lo_beta, hi_beta, sum_weight, H, tmp;
    int p;
    beta = 1;
    lo_beta = hi_beta = -1;

    for (int iter = 0; iter < 200; ++iter) {
      H = sum_weight = 0;
      for (p = head[id]; p >= 0; p = next[p]) {
        sum_weight += tmp = exp(-beta * edge_weight[p]);
        H += beta * (edge_weight[p] * tmp);
      }
      H = (H / sum_weight) + log(sum_weight);
      if (fabs(H - log(perplexity)) < 1e-5) break;
      if (H > log(perplexity)) {
        lo_beta = beta;
        if (hi_beta < 0) beta *= 2; else beta = (beta + hi_beta) / 2;
      } else {
        hi_beta = beta;
        if (lo_beta < 0) beta /= 2; else beta = (lo_beta + beta) / 2;
      }
    }
    for (p = head[id], sum_weight = 0; p >= 0; p = next[p]) {
      sum_weight += edge_weight[p] = exp(-beta * edge_weight[p]);
    }
    for (p = head[id]; p >= 0; p = next[p]){
      edge_weight[p] /= sum_weight;
    }
    sigmas[id] = beta;
  }

  void searchReverse(int id) {
    int x, y, p, q;
    for (p = head[x]; p >= 0; p = next[p]) {
      y = edge_to[p];
      for (q = head[x]; q >= 0; q = next[q]) {
        if (edge_to[q] == x) break;
      }
      reverse[p] = q;
    }
  }

  ReferenceEdges(double perplexity) {
    edge_from = std::vector<int>();
    edge_to = std::vector<int>();
    edge_weight = std::vector<double>();
    this -> perplexity = perplexity;
    head = std::vector<int>();
    next = std::vector<int>();
    reverse = std::vector<int>();
  }

  void build(const arma::ivec& from,
             const arma::ivec& to,
             const arma::vec& weights) {
    n_edges = from.size();
    n_vertices = from[n_edges - 1] + 1;
    sigmas = vec(n_vertices);
    int n_edge = 0;
    for (int i = 0; i < n_vertices; i++) head.push_back(-1);
    for (int x = 0; x < n_vertices; x++) {
      while (from[n_edge] == x) {
        edge_from.push_back(x);
        edge_to.push_back(to[n_edge]);
        edge_weight.push_back(weights[n_edge] * weights[n_edge]);
        next.push_back(head[x]);
        reverse.push_back(-1);
        head[x] = n_edge++;
      }
    }
  }

  arma::sp_mat getWIJ() {
   #pragma omp parallel for
    for (int id = 0; id < n_vertices; id++) {
      similarityOne(id);
    }
   #pragma omp parallel for
    for (int id = 0; id < n_vertices; id++) {
      searchReverse(id);
    }
    int n_edge = edge_to.size();
    double sum_weight = 0;
    for (int id = 0; id != n_vertices; id++) {
      for (int p = head[id]; p >= 0; p = next[p]) {
        int y = edge_to[p];
        int q = reverse[p];
        if (q == -1) {
          edge_from.push_back(y);
          edge_to.push_back(id);
          edge_weight.push_back(0);
          next.push_back(head[y]);
          reverse.push_back(p);
          q = reverse[p] = head[y] = n_edge++;
        }
        if (id > y){
          sum_weight += edge_weight[p] + edge_weight[q];
          edge_weight[p] = edge_weight[q] = (edge_weight[p] + edge_weight[q]) / 2;
        }
      }
    }
    // return sigmas;
    umat locations = umat(2, edge_from.size());
    vec values = vec(edge_weight.size());
    for (int i = 0; i < edge_from.size(); i++) {
      locations(0, i) = edge_from[i];
      locations(1, i) = edge_to[i];
      values[i] = edge_weight[i];
    }
    sp_mat wij = sp_mat(
      true, // add_values
      locations,
      values,
      n_vertices, n_vertices // n_col and n_row
    );
    return wij;
  }

  arma::vec getSigmas() {
    return sigmas;
  }
};

// [[Rcpp::export]]
arma::sp_mat referenceWij(const arma::ivec& i,
                  const arma::ivec& j,
                  arma::vec& d,
                  double perplexity) {
  ReferenceEdges ref = ReferenceEdges(perplexity);
  ref.build(i, j, d);
  vec sigmas = ref.getSigmas();
  sp_mat wij = ref.getWIJ();
  return wij;
  // return List::create(Named("wij", wij),
  //                     Named("sigmas", sigmas));
}
