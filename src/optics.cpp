// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "minindexedpq.h"
#include <queue>
#include <Rmath.h>
#include <progress.hpp>

using namespace Rcpp;
using namespace std;
using namespace arma;

typedef pair<long long, double> iddist;

class CompareDist {
public:
  bool operator()(iddist n1, iddist n2) {
    return n1.second > n2.second;
  }
};

typedef std::priority_queue<iddist,
                            vector<iddist>,
                            CompareDist> NNheap;
typedef std::vector<iddist> NNlist;

class OPTICS {
protected:
	arma::sp_mat* edges;
	arma::imat* neighbors;
	const double eps;
	const long long N;
	vector< bool > visited;
	vector< long long > orderedPoints;
	vector< double > reachdist, coredist;
	vector< long long > predecessor;
	Progress progress;

	long double reachabilityDistance(long long& p,
                                  long long& q) const {
		double dist = max((*edges)(p, q), (*edges)(q, p));
		return max(coredist[p], dist);
	}

	NNlist getNeighbors(long long& p,
                      PairingHeap< long long, double >& seeds) {
		NNlist ret = NNlist(neighbors -> n_rows);
		if (coredist[p] == INFINITY) return ret;
		bool exceeded = false;
		for (auto it = neighbors -> begin_col(p); it != neighbors -> end_col(p) && *it != -1; it++) {
			if (visited[*it]) continue;
			double d = (*edges)(p, *it);
			if (d < eps) addNeighbor(p, *it, seeds);
			else {
				exceeded = true;
				break;
			}
		}
		if (! exceeded) for (auto it = edges -> begin_col(p); it != edges -> end_col(p); it++) {
			if (! visited[it.row()] && *it < eps) {
				addNeighbor(p, it.row(), seeds);
			}
		}
		return ret;
	}

	void addNeighbor(long long& p,
                   long long q,
                   PairingHeap< long long, double>& seeds) {
		if (visited[q]) return;
		double newReachabilityDistance = reachabilityDistance(p, q);

		if(reachdist[q] == INFINITY) {
			seeds.insert(q, newReachabilityDistance);
			predecessor[q] = p;
		} else if (seeds.decreaseIf(q, newReachabilityDistance))  predecessor[q] = p;
	}

public:
	OPTICS(arma::sp_mat& edges,
         arma::imat& neighbors,
         double eps,
         int minPts,
         bool verbose) : edges{&edges}, neighbors{&neighbors},
         								 eps{eps}, N(neighbors.n_cols),
         								 visited(vector< bool >(N, false)),
								         orderedPoints(vector<long long>()),
								         reachdist(vector< double >(N, INFINITY)),
								         coredist(vector< double >(N)),
								         predecessor(vector< long long >(N, NA_INTEGER)),
								         progress(Progress(N, verbose)) {
         	if (neighbors.n_rows < minPts) stop("Insufficient neighbors.");
         	orderedPoints.reserve(N);
         	for (long long n = 0; n != N; n++) {
         		double nthDistance = edges(n, neighbors(minPts - 1, n));
         		coredist[n] = (nthDistance < eps) ? nthDistance : INFINITY;
         	}
         }

	List run() {
		PairingHeap<long long, double> seeds(N);
		for (long long p = 0; p < N; p++) if (progress.increment() && ! visited[p]) {
			visited[p] = true;
			orderedPoints.push_back(p);
			if (p > 0) reachdist[p] = seeds.keyOf(p);
			if (coredist[p] == INFINITY) continue; // core-dist is undefined
			getNeighbors(p, seeds);
			while (!seeds.isEmpty()) {
				long long q = seeds.pop();
				visited[q] = true;
				orderedPoints.push_back(q);
				reachdist[q] = seeds.keyOf(q);
				if (coredist[q] == INFINITY) continue;
				getNeighbors(q, seeds);
			}
		}
		List ret;
		ret["order"] = IntegerVector(orderedPoints.begin(), orderedPoints.end()) +1;
		ret["reachdist"] = NumericVector(reachdist.begin(), reachdist.end());
		ret["coredist"] = NumericVector(coredist.begin(), coredist.end());
		ret["predecessor"] = IntegerVector(predecessor.begin(), predecessor.end()) + 1;
		return ret;
	}
};

// [[Rcpp::export]]
List optics_cpp(arma::sp_mat& edges,
                arma::imat& neighbors,
                double eps,
                int minPts,
                bool verbose) {
	OPTICS opt = OPTICS(edges, neighbors, eps, minPts, verbose);
	Rcout << "\nrun\n";
	return opt.run();
}