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

class HDCluster {
private:
	HDCluster* parent = nullptr;
public:
	HDCluster* left = nullptr;
private:
	HDCluster* right = nullptr;
public:
	const arma::uword sz; // Size at top of cluster
public:
	const arma::uword id;
private:
	list< std::pair<arma::uword, double >> fallenPoints; // Points that leave cluster betweeen top and split
public:
	arma::uword rank = 0;
private:
	double lambda_birth = 0; // 1 / Distance at which splits from parent cluster
	double lambda_death = INFINITY; // 1 / Distance at which cluster splits
	double sum_lambda_p = 0; // sum of lambda_p for all points in cluster, fallen and split
	double stability = 0;
	bool selected = false;

	void deselect();

	void mergeUp();
	void innerCondense(const unsigned int minPts);
	void condenseSingleton();
	void condenseTooSmall();
	void extract( double* ret, arma::uword& selectedClusterCnt, arma::uword currentSelectedCluster, Progress& p) const;
	void reportHierarchy(
							 arma::uword& clusterCnt,
               vector<arma::uword>& nodeMembership, // The clusterid of the immediate parent for each point
               vector<double>& lambdas,
               vector<arma::uword>& clusterParent,
               vector<bool>& clusterSelected,
               vector<double>& clusterStability,
               const arma::uword parentCluster) const;

public:
	void newparent(vector<HDCluster*>& points, HDCluster* newparent);
	~HDCluster();
	explicit HDCluster(const arma::uword& id);

	HDCluster* getRoot();

	HDCluster(HDCluster* a, HDCluster* b, const arma::uword& id, const double& d);

	double determineStability(const unsigned int& minPts, Progress& p);
	void determineSubStability(const unsigned int& minPts, Progress& p);

	void extract(
			double* ret, // An N * 2 array where for each point n * 2 is the cluster id for the point and n * 2 + 1 is lambda_p.
			arma::uword& selectedClusterCnt,
			Progress& p
	) const;

	void reportHierarchy(
			arma::uword& clusterCnt,
			vector<arma::uword>& nodeMembership, // The clusterid of the immediate parent for each point
			vector<double>& lambdas,
			vector<arma::uword>& clusterParent,
			vector<bool>& clusterSelected,
			vector<double>& clusterStability);

	void condense(const unsigned int minPts, unsigned int level);
};

template<typename Tval>
struct MyTemplatePointerHash1 {
	size_t operator()(const Tval* val) const {
		static const size_t shift = (size_t)log2(1 + sizeof(Tval));
		return (size_t)(val) >> shift;
	}
};

typedef unordered_map<arma::uword, HDCluster*> Rootset;

class HDBSCAN {
private:
  arma::uword N;
  Progress p;
  Rootset roots;
  double* coreDistances;

  void buildHierarchy(const vector<pair<double, arma::uword>>& mergeSequence,
                      const unsigned int& minPts,
                      const arma::uword* minimum_spanning_tree);
  void determineStability(const unsigned int& minPts);
  void extractClusters(double* ret);
  void condense(const unsigned int& minPts);
  void condense(const unsigned int& minPts, vector<HDCluster*>& points);
public:
	HDBSCAN(const arma::uword& N, const bool& verbose);
	~HDBSCAN();

	void makeCoreDistances(	const sp_mat& edges,
                          const IntegerMatrix& neighbors,
                          const int& K);
	IntegerVector build( 		const unsigned int& K,
                        	const sp_mat& edges,
                        	const unsigned int& minPts,
                        	const IntegerMatrix& neighbors);
	void condenseAndExtract(const unsigned int& minPts, double* clusters);
	Rcpp::List getHierarchy() const;
};