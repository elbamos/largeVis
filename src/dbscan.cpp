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

class DBSCAN {
protected:
	arma::sp_mat* edges;
	long double eps;
	int minPts;
	long long N;
	std::vector< bool > visited;
	std::vector< int > clusterAssignments;
	Progress progress;

	list< long long > regionQuery(long long& p) const {
		set<long long> trackSet = set<long long>();
		trackSet.insert(p);
		bool exceeded = false;
		for (auto it = edges -> begin_row(p); it != edges -> end_row(p); it++) {
			if (*it < eps) trackSet.insert(it.col());
			else exceeded = true;
		}
		if (! exceeded) {
			for (auto it = edges -> begin_col(p); it != edges -> end_col(p); it++) {
				if (*it < eps) trackSet.insert(it.row());
				else exceeded = true;
			}
		}
		if (! exceeded) Rf_warning("eps not exceeded. Cannot be certain that the full e-neighborhood has been found. Consider reducing eps or providing more neighbors.");
		list<long long> ret = list<long long>(trackSet.begin(), trackSet.end());
		return ret;
	}

	void expandCluster(long long& P, list< long long >& pNeighbors, int C) {
		clusterAssignments[P] = C;
		for (auto pprime = pNeighbors.begin(); pprime != pNeighbors.end(); pprime++) {
			if (! visited[*pprime]) {
				visited[*pprime] = true;
				list< long long > pprimeNeighbors = regionQuery(*pprime);
				if (pprimeNeighbors.size() >= minPts) {
					pNeighbors.insert(pNeighbors.end(), pprimeNeighbors.begin(), pprimeNeighbors.end());
				}
			}
			if (clusterAssignments[*pprime] == -1) clusterAssignments[*pprime] = C;
		}
	}

public:

	DBSCAN( arma::sp_mat& edges,
         const arma::imat& neighbors,
         double eps,
         int minPts,
         bool verbose) : edges{&edges}, eps{eps}, minPts{minPts}, N(neighbors.n_cols), visited(vector< bool >(N, false)),
								         clusterAssignments(vector<int>(N, -1)),
								         progress(Progress(N, verbose)) {
         	if (neighbors.n_rows < minPts) stop("Insufficient neighbors.");
         }


	IntegerVector run() {
		int C = -1;
		for (long long p = 0; p < N; p++) if (progress.increment() && ! visited[p]) {
			visited[p] = true;
			list< long long > pNeighbors = regionQuery(p);
			if (pNeighbors.size() >= minPts) {
				++C;
				expandCluster(p, pNeighbors, C);
			}
		}
		return IntegerVector(clusterAssignments.begin(), clusterAssignments.end()) + 1;
	}
};

// [[Rcpp::export]]
List optics_cpp(arma::sp_mat& edges,
              arma::imat& neighbors,
                double eps,
                int minPts,
                bool verbose) {
	if (neighbors.n_rows < minPts) stop("Insufficient neighbors.");
	long long N = neighbors.n_cols;
	Progress progress(N, verbose);
	std::vector< bool > visited(N, false);
	std::vector<long long> orderedPoints = vector<long long>(), predecessor(N, -1);
	std::vector<double> reachdist(N, INFINITY), coredist(N, INFINITY);
	orderedPoints.reserve(N);

	for (long long n = 0; n != N; n++) {
		double nthDistance = edges(n, neighbors(minPts - 1, n));
		coredist[n] = (nthDistance < eps) ? nthDistance : INFINITY;
	}

	PairingHeap<long long, double> seeds = PairingHeap<long long, double>(N);
	seeds.batchInsert(N, 0);
	while (! seeds.isEmpty()) {
		long long p = seeds.pop();
		orderedPoints.push_back(p);
		visited[p] = true;
		reachdist[p] = seeds.keyOf(p);
		progress.increment();
		if (coredist[p] == INFINITY) continue;
		double cdp = coredist[p];

		for (auto it = edges.begin_row(p); it != edges.end_row(p); it++) {
			double d = *it;
			long long q = it.col();
			if (visited[q]) continue;
			double newReachabilityDistance = max(cdp, d);

			if (seeds.decreaseIf(q, newReachabilityDistance))  predecessor[q] = p;
		}
	}
	List ret;
	ret["order"] = IntegerVector(orderedPoints.begin(), orderedPoints.end()) +1;
	ret["reachdist"] = NumericVector(reachdist.begin(), reachdist.end());
	ret["coredist"] = NumericVector(coredist.begin(), coredist.end());
	ret["predecessor"] = IntegerVector(predecessor.begin(), predecessor.end()) + 1;
	return ret;
}

// [[Rcpp::export]]
IntegerVector dbscan_cpp(arma::sp_mat& edges,
                     arma::imat& neighbors,
                     double eps,
                     int minPts,
                     bool verbose) {
	DBSCAN db = DBSCAN(edges, neighbors, eps, minPts, verbose);
	return db.run();
}

