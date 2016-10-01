// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "minindexedpq.h"
#include <queue>
#include <Rmath.h>
#include <progress.hpp>

//#define DEBUG

using namespace Rcpp;
using namespace std;
using namespace arma;

class DBSCAN {
protected:
	const sp_mat* edges;
	const imat* neighbors;
	long double eps;
	int minPts;
	long long N;
	vector< bool > visited;
	vector< int > clusterAssignments;
	Progress progress;

	list< long long > regionQuery(long long& p) const {
		set< long long > holder = set< long long >();
		bool exceeded = false;
		for (auto it = neighbors -> begin_col(p);
         it != neighbors -> end_col(p);
         it ++) {
			if (*it == -1 || (*edges)(p, *it) > eps) {
				exceeded = true;
				break;
			}
			holder.insert(*it);
		}
		if (! exceeded) {
			for (auto it = edges -> begin_col(p);
        	 it != edges -> end_col(p);
        	 it++) {
				if (*it < eps) holder.insert(it.row());
			}
		}
		list< long long > ret = list< long long >(holder.begin(), holder.end());
		return ret;
	}

	void expandCluster(long long& P, list< long long >& pNeighbors, int& C) {
		clusterAssignments[P] = C;
		for (auto pprime = pNeighbors.begin(); pprime != pNeighbors.end(); pprime++) {
			if (! visited[*pprime]) {
				visited[*pprime] = true;
				list< long long > pprimeNeighbors = regionQuery(*pprime);
				if (pprimeNeighbors.size() >= minPts - 1) {
					pNeighbors.insert(pNeighbors.end(), pprimeNeighbors.begin(), pprimeNeighbors.end());
				}
			}
			if (clusterAssignments[*pprime] == -1) clusterAssignments[*pprime] = C;
		}
	}

public:

	DBSCAN(const sp_mat& edges,
         const imat& neighbors,
         double eps,
         int minPts,
         bool verbose) : edges{&edges}, neighbors{&neighbors},
         								 eps{eps}, minPts{minPts}, N(neighbors.n_cols),
         								 visited(vector< bool >(N, false)),
								         clusterAssignments(vector<int>(N, -1)),
								         progress(Progress(N, verbose)) {
         	if (neighbors.n_rows < minPts) stop("Insufficient neighbors.");

       }


	IntegerVector run() {
		int C = -1;
		for (long long p = 0; p < N; p++) if (progress.increment() && ! visited[p]) {
			visited[p] = true;
			list< long long > pNeighbors = regionQuery(p);
			if (pNeighbors.size() >= minPts - 1) {
				++C;
				expandCluster(p, pNeighbors, C);
			}
		}
		return IntegerVector(clusterAssignments.begin(), clusterAssignments.end()) + 1;
	}
};

// [[Rcpp::export]]
IntegerVector dbscan_cpp(const arma::sp_mat& edges,
                         const arma::imat& neighbors,
                         double eps,
                         int minPts,
                         bool verbose) {
	DBSCAN db = DBSCAN(edges, neighbors, eps, minPts, verbose);
	return db.run();
}