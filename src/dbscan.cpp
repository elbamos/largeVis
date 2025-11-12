#include <RcppArmadillo.h>
#include <Rmath.h>
#include <progress.hpp>

//#define DEBUG

using namespace Rcpp;
using namespace std;
using namespace arma;

typedef sword vertexidx;

class DBSCAN {
protected:
	const sp_mat* edges;
	const imat* neighbors;
	const long double eps;
	const unsigned int minPts;
	const vertexidx N;
	vector< bool > visited;
	vector< int > clusterAssignments;
	Progress progress;

	list< vertexidx > regionQuery(vertexidx& p) const {
		set< vertexidx > holder = set< vertexidx >();
		bool exceeded = false;
		for (auto it = neighbors -> begin_col(p);
         it != neighbors -> end_col(p);
         ++it) {
			if (*it == -1 || (*edges)(p, *it) >= eps) {
				exceeded = true;
				break;
			}
			holder.insert(*it);
		}
		if (! exceeded) {
			for (auto it = edges -> begin_col(p);
        	 it != edges -> end_col(p);
        	 ++it) {
				if (*it < eps) holder.insert(it.row());
			}
		}
		list< vertexidx > ret = list< vertexidx >(holder.begin(), holder.end());
		return ret;
	}

	void expandCluster(vertexidx& P, list< vertexidx >& pNeighbors, int& C) {
		clusterAssignments[P] = C;
		for (auto pprime = pNeighbors.begin(); pprime != pNeighbors.end(); ++pprime) {
			if (! visited[*pprime]) {
				visited[*pprime] = true;
				list< vertexidx > pprimeNeighbors = regionQuery(*pprime);
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
         const double& eps,
         const unsigned int& minPts,
         bool verbose) : edges{&edges}, neighbors{&neighbors},
         								 eps{eps}, minPts{minPts}, N(neighbors.n_cols),
         								 visited(vector< bool >(N, false)),
								         clusterAssignments(vector<int>(N, -1)),
								         progress(Progress(N, verbose)) {
         	if (neighbors.n_rows < minPts) throw Rcpp::exception("Insufficient Neighbors.");
       }


	IntegerVector run() {
		int C = -1;
		for (vertexidx p = 0; p < N; ++p) if (progress.increment() && ! visited[p]) {
			visited[p] = true;
			list< vertexidx > pNeighbors = regionQuery(p);
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
