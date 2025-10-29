#include "edgeweights.h"

//#define DEBUG

using namespace Rcpp;
using namespace std;
using namespace arma;

void ReferenceEdges::run() {
	SimilarityOneWorker sim_worker(this);
	parallelFor(0, n_vertices, sim_worker);

	SearchReverseWorker rev_worker(this);
	parallelFor(0, n_vertices, rev_worker);

  for (vertexidxtype id = 0; id < n_vertices; ++id) {
    searchReverse(id);
  }
  edgeidxtype n_edge = edge_to.size();
  double sum_weight = 0;
  for (vertexidxtype id = 0; id != n_vertices; ++id) {
    for (edgeidxtype p = head[id]; p >= 0; p = next[p]) {
    	vertexidxtype y = edge_to[p];
    	edgeidxtype q = reverse[p];
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

// [[Rcpp::export]]
arma::sp_mat referenceWij(const arma::ivec& i,
				                  const arma::ivec& j,
				                  arma::vec& d,
				                  double perplexity) {
  ReferenceEdges ref = ReferenceEdges(perplexity, i, j, d);
  // vec sigmas = ref.getSigmas();
  ref.run();
  sp_mat wij = ref.getWIJ();
  return wij;
}
