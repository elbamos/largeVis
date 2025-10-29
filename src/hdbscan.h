#include <RcppArmadillo.h>
#include "progress.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;
//#define DEBUG
//#define DEBUG2

class HDCluster {
private:
	HDCluster* parent = nullptr;

	HDCluster* right = nullptr;

	list< std::pair<arma::uword, double >> fallenPoints; // Points that leave cluster betweeen top and split
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
	void extract( int* clusters, double* lambas, int& selectedClusterCnt, int currentSelectedCluster, Progress& p) const;
	void reportHierarchy(
							 int& clusterCnt,
               vector<int>& nodeMembership, // The clusterid of the immediate parent for each point
               vector<double>& lambdas,
               vector<int>& clusterParent,
               vector<bool>& clusterSelected,
               vector<double>& clusterStability,
               vector<double>& lambdaBirth,
               vector<double>& lambdaDeath,
               const int parentCluster) const;

public:
	HDCluster* left = nullptr;
	const arma::uword sz; // Size at top of cluster
	const arma::uword id;
	arma::uword rank = 0;

	void newparent(vector<HDCluster*>& points, HDCluster* newparent);
	~HDCluster();
	explicit HDCluster(const arma::uword& id);

	HDCluster* getRoot();

	HDCluster(HDCluster* a, HDCluster* b, const arma::uword& id, const double& d);

	void condense(const unsigned int minPts, unsigned int level);
	double determineStability(const unsigned int& minPts, Progress& p);
	void determineSubStability(const unsigned int& minPts, Progress& p);

	void extract(
			int* clusters,
			double* lambdas, // An N * 2 array where for each point n * 2 is the cluster id for the point and n * 2 + 1 is lambda_p.
			int& selectedClusterCnt,
			Progress& p
	) const;

	void reportHierarchy(
			int& clusterCnt,
			vector<int>& nodeMembership, // The clusterid of the immediate parent for each point
			vector<double>& lambdas,
			vector<int>& clusterParent,
			vector<bool>& clusterSelected,
			vector<double>& clusterStability,
			vector<double>& lambdaBirth,
			vector<double>& lambdaDeath);
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
  void extractClusters(int* clusters, double* lambdas);
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
	void condenseAndExtract(const unsigned int& minPts, int* clusters, double* lambdas);
	Rcpp::List getHierarchy() const;
};