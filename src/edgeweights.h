#include "largeVis.h"
#include <vector>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

class ReferenceEdges {
protected:
	// arma::vec sigmas;
	const double perplexity;
	const edgeidxtype n_edges;
	const vertexidxtype n_vertices;
	vector< vertexidxtype > edge_from, edge_to;
	vector< edgeidxtype> head, next, reverse;
	vector< double > edge_weight;

public:
	ReferenceEdges(double perplexity,
                const arma::ivec& from,
                const arma::ivec& to,
                const arma::vec& weights) : perplexity{perplexity},
                n_edges(from.size()),
                n_vertices(from[(long) n_edges - 1] + 1),
                edge_from(vector< vertexidxtype >()),
                edge_to(vector< vertexidxtype >()),
                head(vector< edgeidxtype >(n_vertices, -1)),
                next(vector< edgeidxtype >()),
                reverse(vector< edgeidxtype >()),
                edge_weight(vector<double>()) {
                	// sigmas = vec(n_vertices);
                	edgeidxtype n_edge = 0;
                	edge_from.reserve(n_edges);
                	edge_to.reserve(n_edges);
                	edge_weight.reserve(n_edges);
                	next.reserve(n_edges);
                	reverse.reserve(n_edges);
                	for (vertexidxtype x = 0; x < n_vertices; ++x) {
                		while (n_edge < n_edges && from[n_edge] == x) {
                			edge_from.push_back(x);
                			edge_to.push_back(to[n_edge]);
                			edge_weight.push_back(weights[n_edge] * weights[n_edge]);
                			next.push_back(head[x]);
                			reverse.push_back(-1);
                			head[x] = n_edge++;
                		}
                	}
                }

	void similarityOne(vertexidxtype id) {
		double beta, lo_beta, hi_beta, sum_weight, tmp;
		vertexidxtype p;
		beta = 1;
		lo_beta = hi_beta = -1;

		for (int iter = 0; iter < 200; ++iter) {
			double H = sum_weight = 0;
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

	void searchReverse(vertexidxtype id) {
		edgeidxtype p, q;
		for (p = head[id]; p >= 0; p = next[p]) {
			for (q = head[edge_to[p]]; q >= 0; q = next[q]) {
				if (edge_to[q] == id) break;
			}
			reverse[p] = q;
		}
	}

	void run();

	arma::sp_mat getWIJ() {
		umat locations = umat(2, edge_from.size());
		for (vector< vertexidxtype >::size_type i = 0; i < edge_from.size(); ++i) {
			locations(0, i) = edge_from[i];
			locations(1, i) = edge_to[i];
		}
		sp_mat wij = sp_mat(
			true, // add_values
			locations,
			vec(edge_weight),
			n_vertices, n_vertices, true, false
		);
		return wij;
	}
};

class SimilarityOneWorker : public RcppParallel::Worker {
public:
	ReferenceEdges *edges;

	SimilarityOneWorker(ReferenceEdges *edges) : edges {edges} {};

	void operator()(std::size_t begin, std::size_t end) {
		for (vertexidxtype i = begin; i < end; ++i) {
			edges->similarityOne(i);
		}
	}
};

class SearchReverseWorker : public RcppParallel::Worker {
public:
	ReferenceEdges *edges;

	SearchReverseWorker(ReferenceEdges *edges) : edges {edges} {};

	void operator()(std::size_t begin, std::size_t end) {
		for (vertexidxtype i = begin; i < end; ++i) {
			edges->searchReverse(i);
		}
	}
};
