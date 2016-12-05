// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "progress.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;
//#define DEBUG
//#define DEBUG2

template<class T>
void findAndErase(set<T>& theSet, T& theItem) {
	auto it = theSet.find(theItem);
	if (it != theSet.end()) theSet.erase(it);
};

class HDCluster {
private:
	HDCluster* parent = nullptr;
	HDCluster* left = nullptr;
	HDCluster* right = nullptr;
	arma::uword sz; // Size at top of cluster
	arma::uword id;
	set< std::pair<arma::uword, double >> fallenPoints; // Points that leave cluster betweeen top and split
	double lambda_birth = 0; // 1 / Distance at which splits from parent cluster
	double lambda_death = INFINITY; // 1 / Distance at which cluster splits
	double sum_lambda_p = 0; // sum of lambda_p for all points in cluster, fallen and split
	double stability = 0;
	bool selected = false;

	void deselect();

	HDCluster* getRoot();

	void condenseSingleton(HDCluster* p);
	void condenseTooSmall(const bool& l);
	void extract( double* ret, arma::uword& selectedClusterCnt, arma::uword currentSelectedCluster) const;
	void reportHierarchy(
							 arma::uword& clusterCnt,
               vector<arma::uword>& nodeMembership, // The clusterid of the immediate parent for each point
               vector<double>& lambdas,
               vector<arma::uword>& clusterParent,
               vector<bool>& clusterSelected,
               vector<double>& clusterStability,
               const arma::uword parentCluster) const;

public:
	~HDCluster();
	explicit HDCluster(const arma::uword& id);

	HDCluster(HDCluster* point1, HDCluster* point2, set<HDCluster*>& roots, const double& d);

	void condense(const unsigned int& minPts);
	void determineStability(const unsigned int& minPts);

	void extract(
			double* ret, // An N * 2 array where for each point n * 2 is the cluster id for the point and n * 2 + 1 is lambda_p.
			arma::uword& selectedClusterCnt
	) const;

	void reportHierarchy(
			arma::uword& clusterCnt,
			vector<arma::uword>& nodeMembership, // The clusterid of the immediate parent for each point
			vector<double>& lambdas,
			vector<arma::uword>& clusterParent,
			vector<bool>& clusterSelected,
			vector<double>& clusterStability);
};

class HDBSCAN {
private:
  arma::uword N;
  Progress p;
  set<HDCluster*> roots;
  double* coreDistances;

  void buildHierarchy(const vector<pair<double, arma::uword>>& container, const arma::uword* minimum_spanning_tree);
  void condense(const unsigned int& minPts) const;
  void determineStability(const unsigned int& minPts) const;
  void extractClusters(double* ret) const;
public:
	HDBSCAN(const unsigned long long& N, const bool& verbose);
	~HDBSCAN();

	void makeCoreDistances(const sp_mat& edges,
                                 const IntegerMatrix& neighbors,
                                 const int& K);
	IntegerVector build( const unsigned int& K,
                               const sp_mat& edges,
                               const IntegerMatrix& neighbors);
	void condenseAndExtract(const unsigned int& minPts, double* clusters) const;
	Rcpp::List getHierarchy() const;
};