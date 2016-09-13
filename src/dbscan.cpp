// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "minindexedpq.h"
#include <Rmath.h>

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
		list<long long> ret = list<long long>();
		ret.push_back(p);
		bool exceeded = false;
		for (auto it = edges -> begin_row(p); it != edges -> end_row(p); it++) {
			if (*it < eps) ret.push_back(it.col());
			else exceeded = true;
		}
		// If exceeded is false, then the minNeighborhood is bigger than we have.
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

class OPTICS {
protected:
	arma::sp_mat* edges;
	long double eps;
	long long N;
	std::vector< bool > visited;
	std::vector<long long> orderedPoints;
  std::vector<long double> reachdist, coredist;
	Progress progress;

	long double reachabilityDistance(long long& p,
                                   long long& q) const {
		long double dist = max((*edges)(p, q), (*edges)(q, p));
		return max(coredist[p], dist);
	}

	NNlist getNeighbors(long long& p) const {
		NNlist ret = NNlist();
		NNheap heap = NNheap();
		if (coredist[p] == INFINITY) return ret;
		bool exceeded = false;
		for (auto it = edges -> begin_row(p); it != edges -> end_row(p); it++) {
			if (*it < eps) heap.emplace(it.col(), *it);
			else exceeded = true;
		}
		// If exceeded is false, then the minNeighborhood is bigger than we have.
		ret.reserve(heap.size());
		while (! heap.empty()) {
			ret.push_back(heap.top());
			heap.pop();
		}
		return ret;
	}

	void update(NNlist& pNeighbors,
             PairingHeap<long long, long double>& seeds,
             long long& p) {

		while(!pNeighbors.empty()) {
			iddist o = pNeighbors.back();
			pNeighbors.pop_back();

			if (visited[o.first]) continue; // already processed

			long double newReachabilityDistance = reachabilityDistance(p, o.first);

			if(reachdist[o.first] == INFINITY) {
				reachdist[o.first] = newReachabilityDistance;
				seeds.insert(o.first, newReachabilityDistance);
			} else  if(newReachabilityDistance < reachdist[o.first]) {
				reachdist[o.first] = newReachabilityDistance;
				seeds.decreaseIf(o.first, newReachabilityDistance);
			}
		}
	}


public:

  OPTICS( arma::sp_mat& edges,
          const arma::imat& neighbors,
          double eps,
          int minPts,
          bool verbose) : edges{&edges}, eps{eps}, N(neighbors.n_cols), visited(vector< bool >(N, false)),
          								orderedPoints(vector<long long>()),
          								reachdist(vector<long double>(N, INFINITY)),
								          coredist(vector<long double>(N, INFINITY)), progress(Progress(N, verbose)) {
		if (neighbors.n_rows < minPts) stop("Insufficient neighbors.");
    orderedPoints.reserve(N);
  	for (long long n = 0; n != N; n++) {
  		long long nthNeighbor = neighbors(minPts - 1, n);
  		coredist[n] = (edges(n, nthNeighbor) < eps) ? edges(n, nthNeighbor) : INFINITY;
  	}
  }

  List run() {
    for (long long p = 0; p < N; p++) if (progress.increment() && ! visited[p]) {

      NNlist pNeighbors = getNeighbors(p);
    	visited[p] = true;
			orderedPoints.push_back(p);

      if (coredist[p] == INFINITY) continue; // core-dist is undefined
			PairingHeap<long long, long double> seeds = PairingHeap<long long, long double>(N);
			update(pNeighbors, seeds, p);

      while (!seeds.isEmpty()) {
      	long long q = seeds.pop();
      	NNlist qNeighbors = getNeighbors(q);
      	visited[q] = true;
      	orderedPoints.push_back(q);
      	if (coredist[q] != INFINITY) {
      		update(qNeighbors, seeds, q);
      	}
      }
    }
    List ret;
    ret["order"] = IntegerVector(orderedPoints.begin(), orderedPoints.end()) +1;
    ret["reachdist"] = NumericVector(reachdist.begin(), reachdist.end());
    ret["coredist"] = NumericVector(coredist.begin(), coredist.end());
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

// [[Rcpp::export]]
IntegerVector dbscan_cpp(arma::sp_mat& edges,
                     arma::imat& neighbors,
                     double eps,
                     int minPts,
                     bool verbose) {
	DBSCAN db = DBSCAN(edges, neighbors, eps, minPts, verbose);
	return db.run();
}

