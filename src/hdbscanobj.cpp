// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "hdbscan.h"
#include "primsalgorithm.h"
//#define DEBUG

void HDBSCAN::condense(const unsigned int& minPts) {
#ifdef _OPENMP
#pragma omp parallel shared(p)
{
	const int level = std::log2(omp_get_max_threads()) + 1;
#pragma omp master
#else
	const int level = 0;
#endif
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		HDCluster* thisone = it->second;
		thisone->condense(minPts, level);
		p.increment(thisone->sz);
	}
#ifdef _OPENMP
}
#endif
}

void HDBSCAN::condense(const unsigned int& minPts, vector<HDCluster*>& points) {
	//	set<HDCluster*> theseroots;
	//	std::for_each(points.begin(), points.end(), [&theseroots](HDCluster* it) {theseroots.insert(it->getRoot());});
#ifdef _OPENMP
#pragma omp parallel
{
	const int level = std::log2(omp_get_max_threads()) + 1;
#else
	const int level = 0;
#endif
	{
#ifdef _OPENMP
#pragma omp master
#endif
		for (auto it = roots.begin(); it != roots.end(); ++it) {
			HDCluster* thisone = it->second;
			if (thisone->left != nullptr)
#ifdef _OPENMP
#pragma omp task
#endif
			{
				thisone->condense(minPts, level);
				thisone->newparent(points, thisone);
			}
		}
	}
#ifdef _OPENMP
}
#endif
}


void HDBSCAN::determineStability(const unsigned int& minPts) {
	if (roots.size() == 1) {
		HDCluster& root = *(roots.begin()->second);
		root.determineSubStability(minPts, p);
	} else {
		for (auto it = roots.begin(); it != roots.end(); ++it) {
			HDCluster& thisone = *(it->second);
			thisone.determineStability(minPts, p);
			p.increment(thisone.sz);
		}
	}
}






void HDBSCAN::extractClusters(double* ret) {
	arma::uword selectedClusterCnt = 1; //NA_INTEGER;
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		HDCluster& thisone = *(it->second);
		thisone.extract(ret, selectedClusterCnt, p);
		p.increment(thisone.sz);
	}
}




HDBSCAN::HDBSCAN(const arma::uword& N, const bool& verbose) :
	N{N},
	p(Progress(6 * N, verbose)) {
		coreDistances = new double[N];
		roots.reserve(N);
	}

void HDBSCAN::buildHierarchy(const vector<pair<double, arma::uword>>& mergeSequence,
                             const unsigned int& minPts,
                             const arma::uword* minimum_spanning_tree) {
	arma::uword cnt = 0;
	std::vector<HDCluster*> points;
	points.reserve(N);
	std::generate_n(points.end(), N, [&cnt](){return new HDCluster(cnt++);});
	roots.max_load_factor(3);
	for (auto it = mergeSequence.begin(); it != mergeSequence.end();  ++it) if (p.increment()) {
		const arma::uword& n = it -> second;
		if (minimum_spanning_tree[n] == -1) continue;
#ifdef DEBUG
		if (it->first == 0) stop("Zero distance");
		if (it->first == NA_INTEGER) stop("NA distance");
		if (it->first == INFINITY) stop("infinite distance");
#endif
		HDCluster  *a, *b;
		a = points[n]->getRoot();
		b = points[minimum_spanning_tree[n]]->getRoot();
		HDCluster* newparent = new HDCluster(a, b, cnt++, it->first);
		roots.erase(a->id);
		roots.erase(b->id);
		roots.emplace(newparent->id, newparent);
		points[n] = points[minimum_spanning_tree[n]] = newparent;
		if (newparent->rank % 4096 == 0) condense(minPts, points);
	}
	roots.rehash(roots.size());
	}

HDBSCAN::~HDBSCAN() {
	delete[] coreDistances;
}

void HDBSCAN::makeCoreDistances(const sp_mat& edges,
                                const IntegerMatrix& neighbors,
                                const int& K) {
	if (neighbors.nrow() < K) stop("Specified K bigger than the number of neighbors in the adjacency matrix.");
	const IntegerVector kthNeighbors = neighbors.row(K - 1);
	for (arma::uword n = 0; n < N; n++) if (p.increment()) {
		const arma::uword q = kthNeighbors[n];
		if (q == -1 || q == NA_INTEGER) stop("Insufficient neighbors.");
		coreDistances[n] = edges(n, q);
		if (coreDistances[n] == 0) coreDistances[n] = max(edges(q, n), 1e-5);
	}
}

IntegerVector HDBSCAN::build( const unsigned int& K,
                              const sp_mat& edges,
                              const unsigned int& minPts,
                              const IntegerMatrix& neighbors) {
	makeCoreDistances(edges, neighbors, K); // 1 N
	PrimsAlgorithm<arma::uword, double> prim = PrimsAlgorithm<arma::uword, double>(N, coreDistances);
	const arma::uword* minimum_spanning_tree = prim.run(edges, neighbors, p, 0); // 1N
	vector< pair<double, arma::uword> > mergeSequence = prim.getMergeSequence();
	buildHierarchy(mergeSequence, minPts, minimum_spanning_tree); // 1 N
	vector<arma::uword> treevector(minimum_spanning_tree, minimum_spanning_tree + N);
	return IntegerVector(treevector.begin(), treevector.end());
}

void HDBSCAN::condenseAndExtract(const unsigned int& minPts, double* clusters) {
	condense(minPts); // 1 N
	determineStability(minPts); // 1 N
	extractClusters(clusters); // 1 N
};

Rcpp::List HDBSCAN::getHierarchy() const {
	vector<arma::uword> nodemembership(N);
	vector<double> lambdas(N);
	vector<arma::uword> clusterParent;
	vector<bool> clusterSelected;
	vector<double> clusterStability;

	arma::uword clusterCnt = 0;
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		HDCluster& thisone = *(it->second);
		thisone.reportHierarchy(clusterCnt, nodemembership, lambdas, clusterParent, clusterSelected, clusterStability);
		delete &thisone;
	}

	return  List::create(Named("nodemembership") = IntegerVector(nodemembership.begin(), nodemembership.end()),
                      Named("lambda") = NumericVector(lambdas.begin(), lambdas.end()),
                      Named("parent") = IntegerVector(clusterParent.begin(), clusterParent.end()),
                      Named("stability") = NumericVector(clusterStability.begin(), clusterStability.end()),
                      Named("selected") = LogicalVector(clusterSelected.begin(), clusterSelected.end()),
                      Named("coredistances") = wrap(vector<double>(coreDistances, coreDistances + N)));
}