// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "hdbscan.h"
#include "primsalgorithm.h"
#define DEBUG

void HDBSCAN::condense(const unsigned int& minPts) {
#ifdef _OPENMP
	const int level = std::log2(omp_get_max_threads()) + 1;
#pragma omp parallel shared(p)
#else
	const int level = 0;
#endif
	{
#ifdef _OPENMP
#pragma omp master
#endif
		for (auto it = roots.begin(); it != roots.end(); ++it) {
			HDCluster* thisone = it->second;
			thisone->condense(minPts, level, p);
		}
	}
}

void HDCluster::condense(const unsigned int minPts, unsigned int level, Progress& p) {
	if (left == nullptr) return;
	if (level != 0) {
#ifdef _OPENMP
#pragma omp task shared(p)
#endif
		{
			left->condense(minPts, level - 1, p);
		}
		right->condense(minPts, level - 1, p);
#ifdef _OPENMP
#pragma omp taskwait
#endif
	} else {
		left->condense(minPts, level, p);
		right->condense(minPts, level, p);
	}
	// We don't have to check the right side because we always put the bigger child on the right
	if (left->sz < minPts) {
		condenseTooSmall();
		swap(left, right); // right definitely null, left not null but could be big or small
		if (left->sz < minPts) condenseTooSmall();
		else condenseSingleton();
	}
#ifdef DEBUG
	if (left != nullptr && left->sz < minPts) stop("bad left");
	if (right != nullptr && right->sz < minPts) stop("bad right");
	if ((left == nullptr) != (right == nullptr)) stop("Singleton!");
#endif
	p.increment();
}

void HDCluster::mergeUp() {
	if (left->sz == 1) {
		sum_lambda_p += left->lambda_birth;
		fallenPoints.emplace_back(left->id, left->lambda_birth);
	}
	else sum_lambda_p += left->sum_lambda_p;
	fallenPoints.splice(fallenPoints.end(), left->fallenPoints);
}

// Entering function we know that child has a split of its own
void HDCluster::condenseSingleton() {
	mergeUp();
	lambda_death = max(lambda_death, left->lambda_death);
#ifdef DEBUG
	if (lambda_death == INFINITY) stop("max infinity");
#endif

	right = left->right;
	left->right = nullptr;
	HDCluster* keep = left;
	left = keep->left;
	keep->left = nullptr;
	delete keep;
}

// Entering function we know that child has no split of its own
void HDCluster::condenseTooSmall() {
	mergeUp();
#ifdef DEBUG
	if (lambda_death == INFINITY) stop("infinity is too small");
#endif

	left->left = nullptr;
	left->right = nullptr;

	delete left;
	left = nullptr;
}




double HDCluster::determineStability(const unsigned int& minPts, Progress& p) {
#ifdef DEBUG
	if (sz < minPts && parent != nullptr) stop("Condense failed");
#endif
	stability = sum_lambda_p - (lambda_birth * fallenPoints.size());
	if (left == nullptr) { // leaf node
		if (sz >= minPts) selected = true; // Otherwise, this is a parent singleton smaller than minPts.
		p.increment();
	} else {
		const double childStabilities = left->determineStability(minPts, p) +
																		right->determineStability(minPts, p);
		stability += lambda_death * (left->sz + right->sz);

		if (stability > childStabilities) {
			selected = true;
			left->deselect();
			right->deselect();
		} else stability = childStabilities;
	}
	return stability;
}

// Prevents agglomeration in a single cluster.
void HDCluster::determineSubStability(const unsigned int& minPts, Progress& p) {
#ifdef DEBUG
	if (sz < minPts && parent != nullptr) stop("Condense failed");
#endif
	stability = sum_lambda_p - (lambda_birth * fallenPoints.size());
	if (left != nullptr) {
		left->determineStability(minPts, p);
		right->determineStability(minPts, p);
	} else {
		p.increment(sz);
	}
}

void HDBSCAN::determineStability(const unsigned int& minPts) {
	if (roots.size() == 1) {
		HDCluster& root = *(roots.begin()->second);
		root.determineSubStability(minPts, p);
	} else {
		for (auto it = roots.begin(); it != roots.end(); ++it) {
			HDCluster& thisone = *(it->second);
			thisone.determineStability(minPts, p);
		}
	}
}






void HDBSCAN::extractClusters(double* ret) {
	arma::uword selectedClusterCnt = 1; //NA_INTEGER;
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		HDCluster& thisone = *(it->second);
		thisone.extract(ret, selectedClusterCnt, p);
	}
}

void HDCluster::extract(
		double* ret, // An N * 2 array where for each point n * 2 is the cluster id for the point and n * 2 + 1 is lambda_p.
		arma::uword& selectedClusterCnt,
		Progress& p
) const { // tracker of the clusterx
	extract(ret, selectedClusterCnt, NA_REAL, p);
}

void HDCluster::extract( double* ret,
                         arma::uword& selectedClusterCnt,
                         arma::uword currentSelectedCluster,
                         Progress& p) const {
	if (selected) currentSelectedCluster = selectedClusterCnt++;
	std::for_each(fallenPoints.begin(), fallenPoints.end(),
               [&ret, &currentSelectedCluster](const std::pair<arma::uword, double>& it) {
               	ret[it.first * 2] = (currentSelectedCluster == 0) ? NA_REAL : currentSelectedCluster;
               	ret[it.first * 2 + 1] = it.second;
               });
	if (left != nullptr) {
		left->extract(ret, selectedClusterCnt, currentSelectedCluster, p);
		right->extract(ret, selectedClusterCnt, currentSelectedCluster, p);
	} else {
		p.increment(sz);
	}
}




void HDCluster::reportHierarchy(
		arma::uword& clusterCnt,
		vector<arma::uword>& nodeMembership, // The clusterid of the immediate parent for each point
		vector<double>& lambdas,
		vector<arma::uword>& clusterParent,
		vector<bool>& clusterSelected,
		vector<double>& clusterStability) {
	reportHierarchy(clusterCnt, nodeMembership, lambdas, clusterParent, clusterSelected, clusterStability, NA_REAL);
}

void HDCluster::reportHierarchy(
		arma::uword& clusterCnt,
		vector<arma::uword>& nodeMembership, // The clusterid of the immediate parent for each point
		vector<double>& lambdas,
		vector<arma::uword>& clusterParent,
		vector<bool>& clusterSelected,
		vector<double>& clusterStability,
		const arma::uword parentCluster) const {
	arma::uword thisCluster = clusterCnt++;
	std::for_each(fallenPoints.begin(), fallenPoints.end(),
               [&nodeMembership, &lambdas, &thisCluster](const std::pair<arma::uword, double>& it) {
	               	nodeMembership[it.first] = thisCluster;
	               	lambdas[it.first] = it.second;
               });
	clusterParent.emplace_back(parentCluster);
	clusterSelected.push_back(selected);
	clusterStability.emplace_back(stability);
	if (left != nullptr) left->reportHierarchy(clusterCnt, nodeMembership, lambdas, clusterParent, clusterSelected, clusterStability, thisCluster);
	if (right != nullptr) right->reportHierarchy(clusterCnt, nodeMembership, lambdas, clusterParent, clusterSelected, clusterStability, thisCluster);
}






HDCluster::~HDCluster() {
	if (left != nullptr) delete left;
	if (right != nullptr) delete right;
}

HDCluster::HDCluster(const arma::uword& id) : sz(1), id(id) { }

HDCluster::HDCluster(HDCluster* a, HDCluster* b, const arma::uword& id, const double& d) :
	sz(a->sz + b->sz), id(id), lambda_birth(0), lambda_death(1/d) {
#ifdef DEBUG
	if (lambda_death == INFINITY) stop("death is infinity");
#endif
	a->parent = b->parent = this;
	a->lambda_birth = b->lambda_birth = lambda_death;

	if (a->sz < b->sz) {
		left = a;
		right = b;
	} else {
		right = a;
		left = b;
	}
}





HDCluster* HDCluster::getRoot() {
	if (parent == nullptr) return this;
	return parent->getRoot();
}

void HDCluster::deselect() {
	if (selected) selected = false;
	else if (left != nullptr) {
		left->deselect();
		right->deselect();
	}
}






HDBSCAN::HDBSCAN(const arma::uword& N, const bool& verbose) :
	N{N},
	p(Progress(10 * N, verbose)) {
		coreDistances = new double[N];
		roots.reserve(N);
}

void HDBSCAN::buildHierarchy(const vector<pair<double,
                             arma::uword>>& mergeSequence,
                             const arma::uword* minimum_spanning_tree) {
	arma::uword cnt = 0;
	std::vector<HDCluster*> points;
	points.reserve(N);
	for (cnt = 0; cnt != N; ++cnt) {
		HDCluster* ret = new HDCluster(cnt);
		points.push_back(ret);
		roots.insert({cnt, ret});
	}
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
                              const IntegerMatrix& neighbors) {
	makeCoreDistances(edges, neighbors, K); // 1 N
	PrimsAlgorithm<arma::uword, double> prim = PrimsAlgorithm<arma::uword, double>(N, coreDistances);
	const arma::uword* minimum_spanning_tree = prim.run(edges, neighbors, p, 0); // 1N
	vector< pair<double, arma::uword> > mergeSequence = prim.getMergeSequence();
	buildHierarchy(mergeSequence, minimum_spanning_tree); // 2 N
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

// [[Rcpp::export]]
List hdbscanc(const arma::sp_mat& edges,
              const IntegerMatrix& neighbors,
              const int& K,
              const int& minPts,
              const Rcpp::Nullable<Rcpp::NumericVector> threads,
              const bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
	HDBSCAN object = HDBSCAN(edges.n_cols, verbose);
	// 1 N
	IntegerVector tree = object.build(K, edges, neighbors); // 4N
	NumericMatrix clusters = NumericMatrix(2, edges.n_cols);
	object.condenseAndExtract(minPts, REAL(clusters)); // 3N
	List hierarchy = object.getHierarchy();
	return List::create(Named("clusters") = clusters,
                      Named("tree") = IntegerVector(tree),
                      Named("hierarchy") = hierarchy);
}
