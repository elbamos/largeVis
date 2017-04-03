// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "hdbscan.h"
#include "primsalgorithm.h"
//#define DEBUG

void HDCluster::newparent(vector<HDCluster*>& points, HDCluster* newparent) {
	for (auto it = fallenPoints.begin(); it != fallenPoints.end(); ++it) {
		points[it->first] = newparent;
	}
	if (left != nullptr) {
		left->newparent(points, newparent);
		right->newparent(points, newparent);
	}
}


void HDCluster::condense(const unsigned int minPts, unsigned int level) {
	if (left != nullptr) {
		const unsigned int newlevel = (level == 0) ? level : level - 1;
#ifdef _OPENMP
#pragma omp task
#endif
{
	left->condense(minPts, newlevel);
}
right->condense(minPts, newlevel);
#ifdef _OPENMP
#pragma omp taskwait
#endif
innerCondense(minPts);
	}
}

void HDCluster::innerCondense(const unsigned int minPts) {
	if (left->sz < minPts) {
		condenseTooSmall();
		swap(left, right); // right definitely null, left not null but could be big or small
		if (left->sz < minPts) {
			condenseTooSmall();
			rank = 0;
		} else {
			condenseSingleton();
			rank = (left == nullptr) ? 0 : max(left->rank, right->rank) + 1;
		}
	}
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
	if (lambda_death == INFINITY) throw Rcpp::exception("max infinity");
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
	if (lambda_death == INFINITY) throw Rcpp::exception("Infinity is too small.");
#endif

	left->left = nullptr;
	left->right = nullptr;

	delete left;
	left = nullptr;
}




double HDCluster::determineStability(const unsigned int& minPts, Progress& p) {
#ifdef DEBUG
	if (sz < minPts && parent != nullptr) throw Rcpp::exception("Condense failed.");
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
	if (sz < minPts && parent != nullptr) throw Rcpp::exception("Condense failed.");
#endif
	stability = sum_lambda_p - (lambda_birth * fallenPoints.size());
	if (left != nullptr) {
			const double childStabilities = left->determineStability(minPts, p) +
																			right->determineStability(minPts, p);
			if (childStabilities > stability) stability = childStabilities;
	} else {
		p.increment(sz);
	}
}


void HDCluster::extract(
		int* clusters,
		double* lambdas, // An N * 2 array where for each point n * 2 is the cluster id for the point and n * 2 + 1 is lambda_p.
		int& selectedClusterCnt,
		Progress& p
) const { // tracker of the clusterx
	extract(clusters, lambdas, selectedClusterCnt, NA_INTEGER, p);
}

void HDCluster::extract( int* clusters,
                         double* lambdas,
                         int& selectedClusterCnt,
                         int currentSelectedCluster,
                         Progress& p) const {
	if (selected) currentSelectedCluster = selectedClusterCnt++;
	std::for_each(fallenPoints.begin(), fallenPoints.end(),
               [&clusters, &lambdas, &currentSelectedCluster](const std::pair<arma::uword, double>& it) {
               	clusters[it.first] = (currentSelectedCluster == 0) ? NA_INTEGER : currentSelectedCluster;
               	lambdas[it.first] = it.second;
               });
	if (left != nullptr) {
		left->extract(clusters, lambdas, selectedClusterCnt, currentSelectedCluster, p);
		right->extract(clusters, lambdas, selectedClusterCnt, currentSelectedCluster, p);
	} else {
		p.increment(sz);
	}
}




void HDCluster::reportHierarchy(
		int& clusterCnt,
		vector<int>& nodeMembership, // The clusterid of the immediate parent for each point
		vector<double>& lambdas,
		vector<int>& clusterParent,
		vector<bool>& clusterSelected,
		vector<double>& clusterStability,
		vector<double>& lambdaBirth,
		vector<double>& lambdaDeath) {
	reportHierarchy(clusterCnt, nodeMembership, lambdas, clusterParent, clusterSelected, clusterStability,
                 lambdaBirth, lambdaDeath, NA_REAL);

}

void HDCluster::reportHierarchy(
		int& clusterCnt,
		vector<int>& nodeMembership, // The clusterid of the immediate parent for each point
		vector<double>& lambdas,
		vector<int>& clusterParent,
		vector<bool>& clusterSelected,
		vector<double>& clusterStability,
		vector<double>& lambdaBirth,
		vector<double>& lambdaDeath,
		const int parentCluster) const {
	int thisCluster = clusterCnt++;
	std::for_each(fallenPoints.begin(), fallenPoints.end(),
               [&nodeMembership, &lambdas, &thisCluster](const std::pair<arma::uword, double>& it) {
               	nodeMembership[it.first] = thisCluster;
               	lambdas[it.first] = it.second;
               });
	clusterParent.emplace_back(parentCluster);
	clusterSelected.push_back(selected);
	clusterStability.emplace_back(stability);
	lambdaBirth.emplace_back(lambda_birth);
	lambdaDeath.emplace_back(lambda_death);
	if (left != nullptr) left->reportHierarchy(clusterCnt, nodeMembership, lambdas,
     clusterParent, clusterSelected, clusterStability, lambdaBirth, lambdaDeath, thisCluster);
	if (right != nullptr) right->reportHierarchy(clusterCnt, nodeMembership, lambdas,
     clusterParent, clusterSelected, clusterStability, lambdaBirth, lambdaDeath,  thisCluster);
}







HDCluster::~HDCluster() {
	if (left != nullptr) delete left;
	if (right != nullptr) delete right;
}

HDCluster::HDCluster(const arma::uword& id) : sz(1), id(id) { }

HDCluster::HDCluster(HDCluster* a, HDCluster* b, const arma::uword& id, const double& d) :
	sz(a->sz + b->sz), lambda_birth(0), lambda_death(1/d), id(id), rank(max(a->rank, b->rank) + 1) {
#ifdef DEBUG
	if (lambda_death == INFINITY) throw Rcpp::exception("death is infiinity.");
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



