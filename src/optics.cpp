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

//#define DEBUG

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

#ifdef DEBUG
set< long long > debugPts = set< long long >();

void startDebug() {
	debugPts.insert(4);
}

bool testDebug(long long& p) {
	if (debugPts.find(p) != debugPts.end()) return true;
	else return false;
}
#endif

	long double reachabilityDistance(long long& p,
                                   long long& q) const {
		double dist = max((*edges)(p, q), (*edges)(q, p));
		return max(coredist[p], dist);
	}

	void getNeighbors(long long& p,
                      PairingHeap< long long, double >& seeds) {
		bool exceeded = false;
		arma::sp_colvec pEdges = edges->col(p);
		for (auto it = neighbors -> begin_col(p);
       	 (it != neighbors -> end_col(p)) && (*it != -1);
       	 it++) if (! visited[*it]) {
			long long q = *it;
			if (pEdges[q] < eps) addNeighbor(p, q, seeds);
			else {
				exceeded = true;
				break;
			}
		}
		if (! exceeded) for (auto it = edges -> begin_col(p);
                         it != edges -> end_col(p);
                         it++) {
			if (! visited[it.row()] && *it < eps) {
				addNeighbor(p, it.row(), seeds);
			}
		}
	}

	void addNeighbor(long long& p,
                   long long q,
                   PairingHeap< long long, double>& seeds) {
		if (visited[q]) return;
		double newReachabilityDistance = reachabilityDistance(p, q);

		if (! seeds.contains(q)) {
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
         	if (minPts < 2) stop("minPts must be >= 2");
         	orderedPoints.reserve(N);
         	for (long long n = 0; n != N; n++) {
         		double nthDistance = edges(n, neighbors(minPts - 2, n));
         		coredist[n] = (nthDistance < eps) ? nthDistance : INFINITY;
         	}
        }

	List run() {
		PairingHeap<long long, double> seeds(N);
		for (long long p = 0; p < N; p++) if (progress.increment() && ! visited[p]) {
			visited[p] = true;
			orderedPoints.push_back(p);
			if (coredist[p] == INFINITY) continue; // core-dist is undefined
			getNeighbors(p, seeds);
			while (!seeds.isEmpty()) {
				long long q = seeds.pop();
				double key = seeds.keyOf(q);
				if (key == seeds.topKey()) {
					long long r = seeds.pop();
					if (r > q) swap(q, r);
					seeds.insert(r, key);
				}
				visited[q] = true;
				orderedPoints.push_back(q);
				reachdist[q] = key;
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
	return opt.run();
}