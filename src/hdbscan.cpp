// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "hdbscan.h"
#include "primsalgorithm.h"
//#define DEBUG

void HDBSCAN::condense(const unsigned int& minPts) const {
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		HDCluster& thisone = **it;
		thisone.condense(minPts);
	}
}

void HDCluster::condense(const unsigned int& minPts) {
	if (left != nullptr) left->condense(minPts);
	if (right != nullptr) right->condense(minPts);
	if (left != nullptr && left->sz < minPts) condenseTooSmall(true);
	if (right != nullptr && right->sz < minPts)	condenseTooSmall(false);
	if (left == nullptr && right != nullptr) {
		HDCluster* p = right;
		condenseSingleton(p);
		delete p;
	}
	if (right == nullptr && left != nullptr) {
		HDCluster* p = left;
		condenseSingleton(p);
//		delete p;
	}
#ifdef DEBUG
	if ((left == nullptr) != (right == nullptr)) stop("Singleton!");
#endif
}

// Entering function we know that child has a split of its own
void HDCluster::condenseSingleton(HDCluster* p) {
	if (p->sz == 1) {
		sum_lambda_p += p->lambda_birth;
		fallenPoints.emplace(p->id, p->lambda_birth);
	}
	else sum_lambda_p += p->sum_lambda_p;
	fallenPoints.insert(p->fallenPoints.begin(), p->fallenPoints.end());
	lambda_death = max(lambda_death, p->lambda_death);
#ifdef DEBUG
	if (lambda_death == INFINITY) stop("max infinity");
#endif
	left = p->left;
	right = p->right;
	if (left != nullptr) left->parent = this;
	if (right != nullptr) right->parent = this;
}

// Entering function we know that child has no split of its own
void HDCluster::condenseTooSmall(const bool& l) {
	HDCluster* p = (l) ? left : right;
	if (p->sz == 1) {
		sum_lambda_p += p->lambda_birth;
		fallenPoints.emplace(p->id, p->lambda_birth);
	}
	else sum_lambda_p += p->sum_lambda_p;
	fallenPoints.insert(p->fallenPoints.begin(), p->fallenPoints.end());
//	lambda_death = max(lambda_death, p->lambda_death);
#ifdef DEBUG
	if (lambda_death == INFINITY) stop("infinity is too small");
#endif
	if (l) left = nullptr;
	else right = nullptr;
	delete p;
}




void HDCluster::determineStability(const unsigned int& minPts) {
	if (sz < minPts && parent != nullptr) stop("Condense failed");
#ifdef DEBUG
	if (stability == INFINITY) stop("Unstable");
#endif
	stability = sum_lambda_p - (lambda_birth * fallenPoints.size());
#ifdef DEBUG
	if (stability == INFINITY) stop("Unstable 2");
#endif
	if (left == nullptr) { // leaf node
		if (sz >= minPts) selected = true; // Otherwise, this is a parent singleton smaller than minPts.
		else stop("Root merge failure - provide more neighbors or reduce minPts.");
	} else {
		left->determineStability(minPts);
		right->determineStability(minPts);
		stability += lambda_death * (left->sz + right->sz);
#ifdef DEBUG
		if (stability == INFINITY) stop("Unstable 3");
#endif
		double childStabilities = left->stability + right->stability;
		if (stability > childStabilities) {
			selected = true;
			left->deselect();
			right->deselect();
		} else stability = childStabilities;
	}
}

void HDBSCAN::determineStability(const unsigned int& minPts) const {
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		HDCluster& thisone = **it;
		thisone.determineStability(minPts);
	}
}






void HDBSCAN::extractClusters(double* ret) const {
	arma::uword selectedClusterCnt = 0;
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		HDCluster& thisone = **it;
		thisone.extract(ret, selectedClusterCnt);
	}
}

void HDCluster::extract(
		double* ret, // An N * 2 array where for each point n * 2 is the cluster id for the point and n * 2 + 1 is lambda_p.
		arma::uword& selectedClusterCnt
) const { // tracker of the clusterx
	extract(ret, selectedClusterCnt, NA_REAL);
}

void HDCluster::extract( double* ret, arma::uword& selectedClusterCnt, arma::uword currentSelectedCluster) const {
	if (selected) currentSelectedCluster = selectedClusterCnt++;
	std::for_each(fallenPoints.begin(), fallenPoints.end(),
               [&ret, &currentSelectedCluster](const std::pair<arma::uword, double>& it) {
               	ret[it.first * 2] = currentSelectedCluster;
               	ret[it.first * 2 + 1] = it.second;
               });
	if (left != nullptr) left->extract(ret, selectedClusterCnt, currentSelectedCluster);
	if (right != nullptr) right->extract(ret, selectedClusterCnt, currentSelectedCluster);
}




void HDCluster::reportHierarchy(
		arma::uword& clusterCnt,
		vector<arma::uword>& nodeMembership, // The clusterid of the immediate parent for each point
		vector<double>& lambdas,
		vector<arma::uword>& clusterParent,
		vector<bool>& clusterSelected,
		vector<double>& clusterStability) {
	reportHierarchy(clusterCnt, nodeMembership, lambdas, clusterParent, clusterSelected, clusterStability, NA_INTEGER);
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
               	nodeMembership[it.first] = (thisCluster == NA_INTEGER) ? NA_REAL : thisCluster;
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

HDCluster::HDCluster(const arma::uword& id) : sz(1), id(id) {}

HDCluster::HDCluster(HDCluster* point1, HDCluster* point2, set<HDCluster*>& roots, const double& d) :
	id(0), lambda_birth(1/d), lambda_death(1/d) {
#ifdef DEBUG
	if (lambda_death == INFINITY) stop("death is infinity");
#endif
	HDCluster* a = point1->getRoot();
	HDCluster* b = point2->getRoot();
	sz = a->sz + b->sz;
	a->parent = this;
	b->parent = this;
	if (d == 0 || (1 / d) == INFINITY) stop("Infinite lambda");
	a->lambda_birth = b->lambda_birth = lambda_death;
	if (a->sz > b->sz) {
		left = a;
		right = b;
	} else {
		right = a;
		left = b;
	}
	findAndErase(roots, left);
	findAndErase(roots, right);
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






HDBSCAN::HDBSCAN(const unsigned long long& N, const bool& verbose) : N{N}, p(Progress(10 * N, verbose)) {
	coreDistances = new double[N];
}

void HDBSCAN::buildHierarchy(const vector<pair<double, arma::uword>>& container, const arma::uword* minimum_spanning_tree) {
	vector<HDCluster*> points = vector<HDCluster*>();
	points.reserve(N);
	arma::uword cnt = 0;
	std::generate_n(points.begin(), N, [&cnt](){return new HDCluster(cnt++);});
	roots.insert(points.begin(), points.end());

	for (auto it = container.begin(); it != container.end();  ++it) {
		arma::uword n = it -> second;
		if (minimum_spanning_tree[n] == -1) continue;
#ifdef DEBUG
		if (it->first == 0) stop("Zero distance");
		if (it->first == NA_INTEGER) stop("NA distance");
		if (it->first == INFINITY) stop("infinite distance");
#endif
		HDCluster* newparent = new HDCluster(points[n], points[minimum_spanning_tree[n]], roots, it->first);
		roots.insert(newparent);
	}
}

HDBSCAN::~HDBSCAN() {
	delete[] coreDistances;
}

void HDBSCAN::makeCoreDistances(const sp_mat& edges,
                                const IntegerMatrix& neighbors,
                                const int& K) {
	if (neighbors.nrow() < K) stop("Specified K bigger than the number of neighbors in the adjacency matrix.");
	//if (K < 4) stop("K must be >= 4 when used with neighbors.");
	IntegerVector kthNeighbors = neighbors.row(K - 1);
	for (arma::uword n = 0; n < N; n++) if (p.increment()) {
		arma::uword q = kthNeighbors[n];
		if (q == -1 || q == NA_INTEGER) stop("Insufficient neighbors.");
		coreDistances[n] = edges(n, q);
		if (coreDistances[n] == 0) coreDistances[n] = max(edges(q, n), 1e-5);
	}
}

IntegerVector HDBSCAN::build( const unsigned int& K,
                              const sp_mat& edges,
                              const IntegerMatrix& neighbors) {
	makeCoreDistances(edges, neighbors, K);
	PrimsAlgorithm<arma::uword, double> prim = PrimsAlgorithm<arma::uword, double>(N, coreDistances);
	const arma::uword* minimum_spanning_tree = prim.run(edges, neighbors, p, 0);
	vector< pair<double, arma::uword> > mergeSequence = prim.getMergeSequence();
	buildHierarchy(mergeSequence, minimum_spanning_tree); // 2 N
	vector<arma::uword> treevector(minimum_spanning_tree, minimum_spanning_tree + N);
	return IntegerVector(treevector.begin(), treevector.end());
}

void HDBSCAN::condenseAndExtract(const unsigned int& minPts, double* clusters) const {
	condense(minPts); // 2 N
	determineStability(minPts); // N
	extractClusters(clusters );
};

Rcpp::List HDBSCAN::getHierarchy() const {
	vector<arma::uword> nodemembership(N);
	vector<double> lambdas(N);

	vector<arma::uword> clusterParent;
	vector<bool> clusterSelected;
	vector<double> clusterStability;

	arma::uword clusterCnt = 0;
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		HDCluster& thisone = **it;
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
              Rcpp::Nullable<Rcpp::NumericVector> threads,
              const bool& verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
	HDBSCAN object = HDBSCAN(edges.n_cols, verbose);
	// 1 N
	IntegerVector tree = object.build(K, edges, neighbors);
	NumericMatrix clusters = NumericMatrix(2, edges.n_cols);
	object.condenseAndExtract(minPts, REAL(clusters));
	List hierarchy = object.getHierarchy();
	return List::create(Named("clusters") = clusters,
                      Named("tree") = IntegerVector(tree),
                      Named("hierarchy") = hierarchy);
}
