// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

class ReferenceEdges {
protected:
  // arma::vec sigmas;
  const double perplexity;
	const long long n_edges;
	const int n_vertices;
  std::vector<long long> edge_from, edge_to, head, next, reverse;
  std::vector<double> edge_weight;

public:
	ReferenceEdges(double perplexity,
                 const arma::ivec& from,
                 const arma::ivec& to,
                 const arma::vec& weights) : perplexity{perplexity},
                 														 n_edges(from.size()),
                                             n_vertices(from[(long) n_edges - 1] + 1),
																						 edge_from(std::vector<long long>()),
																						 edge_to(std::vector<long long>()),
																						 head(std::vector<long long>()),
																						 next(std::vector<long long>()),
																						 reverse(std::vector<long long>()),
																						 edge_weight(std::vector<double>()) {
		// sigmas = vec(n_vertices);
		long n_edge = 0;
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
  void similarityOne(long id) {
    double beta, lo_beta, hi_beta, sum_weight, H, tmp;
  	long p;
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
    // sigmas[id] = beta;
  }

  void searchReverse(int id) {
    long long y, p, q;
    for (p = head[id]; p >= 0; p = next[p]) {
      y = edge_to[p];
      for (q = head[id]; q >= 0; q = next[q]) {
        if (edge_to[q] == id) break;
      }
      reverse[p] = q;
    }
  }

  void run() {
   #pragma omp parallel for
    for (int id = 0; id < n_vertices; id++) {
      similarityOne(id);
    }
   #pragma omp parallel for
    for (int id = 0; id < n_vertices; id++) {
      searchReverse(id);
    }
    long long n_edge = edge_to.size();
    double sum_weight = 0;
    for (int id = 0; id != n_vertices; id++) {
      for (long long p = head[id]; p >= 0; p = next[p]) {
      	long long y = edge_to[p];
      	long long q = reverse[p];
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
  }

  arma::sp_mat getWIJ() {
    umat locations = umat(2, edge_from.size());
    vec values = vec(edge_weight.size());
    for (long long i = 0; i < edge_from.size(); i++) {
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

  // arma::vec getSigmas() {
  //   return sigmas;
  // }
};

// [[Rcpp::export]]
arma::sp_mat referenceWij(const arma::ivec& i,
				                  const arma::ivec& j,
				                  arma::vec& d,
				                  double perplexity) {
  ReferenceEdges ref = ReferenceEdges(perplexity, i, j, d);
  // vec sigmas = ref.getSigmas();
  ref.run();
  return ref.getWIJ();
  // return wij;
}
