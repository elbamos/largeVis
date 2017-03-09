// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
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
	const sp_mat* edges;
	const imat* neighbors;
	const double eps;
	const long long N;

	bool* visited;
	vector< long long > orderedPoints;
	vector< double > reachdist, coredist;
	priority_queue< pair<double, long> > seedQueue;
	vector< long long > predecessor;

	Progress progress;

	const long double reachabilityDistance(const long long& p,
                                  			 const long long& q) const {
		double dist = max((*edges)(p, q), (*edges)(q, p));
		return max(coredist[p], dist);
	}

	void getNeighbors(const long long& p,
                    PairingHeap< long long, double >& seeds) {
		bool exceeded = false;
		const sp_colvec pEdges = edges->col(p);
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

	void addNeighbor(const long long& p,
                   const long long& q,
                   PairingHeap< long long, double>& seeds) {
		if (visited[q]) return;
		const double newReachabilityDistance = reachabilityDistance(p, q);

		if (! seeds.contains(q)) {
			seeds.insert(q, newReachabilityDistance);
			predecessor[q] = p;
		} else if (seeds.decreaseIf(q, newReachabilityDistance))  predecessor[q] = p;
	}

public:
	OPTICS(const sp_mat& edges,
         const imat& neighbors,
         const double& eps,
         const unsigned int& minPts,
         const bool& verbose) : edges{&edges}, neighbors{&neighbors},
         								 eps{eps}, N(neighbors.n_cols),
         								 visited(new bool[N]),
								         orderedPoints(vector<long long>()),
								         reachdist(vector< double >(N, INFINITY)),
								         coredist(vector< double >(N)),
								         predecessor(vector< long long >(N, NA_INTEGER)),
								         progress(Progress(N, verbose)) {
         	if (neighbors.n_rows < minPts) throw Rcpp::exception("Insufficient neighbors.");
         	if (minPts < 2) throw Rcpp::exception("minPts must be >= 2");
         	orderedPoints.reserve(N);
         	for (long long n = 0; n != N; n++) {
         		double nthDistance = edges(n, neighbors(minPts - 2, n));
         		visited[n] = FALSE;
         		coredist[n] = (nthDistance < eps) ? nthDistance : INFINITY;
         	}
        }

	~OPTICS() {
		delete[] visited;
	}

	void queue() {
		for (long long n = 0; n != N; n++) {
			seedQueue.emplace(coredist[n], n);
		}
	}

	inline void runOne(const long long &p, PairingHeap<long long, double> &seeds) {
		visited[p] = true;
		orderedPoints.push_back(p);
		if (coredist[p] == INFINITY) return; // core-dist is undefined
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

	void runAll() {
		PairingHeap<long long, double> seeds(N);
		for (long long p = 0; p != N && progress.increment(); p++) {
			if (! visited[p]) runOne(p, seeds);
		}
	}

	void runQueue() {
		PairingHeap<long long, double> seeds(N);
		while (! seedQueue.empty() && progress.increment()) {
			const long long p = seedQueue.top().second;
			seedQueue.pop();
			if (visited[p]) continue;
			runOne(p, seeds);
		}
	}

	List run() {
		if (seedQueue.empty()) runAll();
		else runQueue();
		List ret;
		ret["order"] = IntegerVector(orderedPoints.begin(), orderedPoints.end()) +1;
		ret["reachdist"] = NumericVector(reachdist.begin(), reachdist.end());
		ret["coredist"] = NumericVector(coredist.begin(), coredist.end());
		ret["predecessor"] = IntegerVector(predecessor.begin(), predecessor.end()) + 1;
		return ret;
	}
};

// [[Rcpp::export]]
List optics_cpp(const arma::sp_mat& edges,
                const arma::imat& neighbors,
                const double& eps,
                const int& minPts,
                const bool& useQueue,
                const bool& verbose) {
	OPTICS opt = OPTICS(edges, neighbors, eps, minPts, verbose);
	if (useQueue) opt.queue();
	return opt.run();
}