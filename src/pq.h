// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "minindexedpq.h"
#include <queue>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "progress.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;
//#define DEBUG
//#define DEBUG2

template<class VIDX>
class UF {
private:
  VIDX reservesize; // used during initialization
	// Used after buildHierarchy to manage condensation, stability extraction, and cluster identification
  VIDX* parents;
  double* lambda_deaths;
  double* stabilities;
  VIDX* sizes;
  set< VIDX >* fallenPointses;
  set< VIDX >* goodChildrens;
  bool* selected;

protected:
  VIDX N;
  Progress p;
  double* lambda_births;

  VIDX starterIndex = 0;
  PairingHeap<VIDX, double> Q;

  VIDX counter = 0;

  set< VIDX > survivingClusters;

  // Used after condensation to speed-up remaining phases
  set< VIDX > roots;

  UF(VIDX N, bool verbose, bool pq) : N{N}, p(Progress(10 * N, verbose)),
  	Q(PairingHeap<VIDX, double>(N)) {
#ifdef DEBUG
  	setupTest();
#endif
  }
  virtual ~UF() {
  	delete[] fallenPointses;
  	delete[] goodChildrens;
  	delete[] selected;
  	delete[] parents;
  	delete[] sizes;
  	delete[] lambda_births;
  	delete[] lambda_deaths;
  	delete[] stabilities;
  }

  /*
   * The hierarchy construction phase
   */
  VIDX getRoot(const VIDX& p) const {
    return (parents[p] == p) ? p : getRoot(parents[p]);
  }

  VIDX add() {
    VIDX newid = counter++;
    parents[newid] = newid;
    sizes[newid] = 1;
    lambda_births[newid] = 0;
    lambda_deaths[newid] = 0;
    return newid;
  }

	/*
	 * The union part of union-find
	 */
  virtual void agglomerate(const VIDX& a, const VIDX& b, const double& d) {
    if (b == -1 || d == INFINITY) {
      lambda_births[a] = (lambda_births[a] > 0) ? lambda_births[a] : -1;
      return;
    } else if (lambda_births[b] == -1) lambda_births[b] = 1 / d;
    VIDX n_a, n_b;
    n_a = getRoot(a);
    n_b = getRoot(b);
    if (n_a == n_b) return;
    VIDX parent = add();
    parents[n_a] = parent;
    parents[n_b] = parent;
    sizes[parent] = sizes[n_a] + sizes[n_b];
    if (d == 0 || (1 / d) == INFINITY) stop("Infinite lambda");
    lambda_births[n_a] = lambda_births[n_b] = lambda_deaths[parent] = 1 / d;
  }

  void buildHierarchy(const vector<pair<double, VIDX>>& container, const VIDX* minimum_spanning_tree) {
  	reservesize = 2 * N + 1;
  	parents = new VIDX[reservesize];
  	sizes = new VIDX[reservesize];
  	lambda_births = new double[reservesize];
  	lambda_deaths = new double[reservesize];

  	for (VIDX n = 0; n != N; ++n) add();

  	for (auto it = container.begin(); it != container.end(); it++) {
  		VIDX n = it -> second;
  		agglomerate(n, minimum_spanning_tree[n], it -> first);
  	}
  }

  /*
   * The condensation phase.
   *
   * First, breadthwise search for leaves with < minPts nodes. While doing that,
   * identify the roots.
   *
   * Then, top-down search for nodes with only one child, iterating through roots.
   *
   */
  void condenseUp(const VIDX& p, const VIDX& mergeTarget) {
    if (p == mergeTarget) return;
    sizes[p] = -1;
    for (auto it = fallenPointses[p].begin();
         it != fallenPointses[p].end();
         it++) {
      parents[*it] = mergeTarget;
      fallenPointses[mergeTarget].insert(*it);
    }
    for (auto it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) {
      goodChildrens[mergeTarget].insert(*it);
      parents[*it] = mergeTarget;
    }
    fallenPointses[p].clear();
    goodChildrens[p].clear();
    lambda_deaths[parents[p]] = max(lambda_deaths[mergeTarget],
                                    lambda_deaths[p]);
    parents[p] = -1;
    return;
  }

  void condenseOne(const VIDX& p, const int& minPts) {
  	if (sizes[p] == -1) return;
    if (parents[p] ==  p) {
      roots.insert(p);
      return;
    }
    if (sizes[p] < minPts) {
      condenseUp(p, parents[p]);
    } else if (parents[p] != p) goodChildrens[parents[p]].insert(p);
  }

  void condenseTwo(const VIDX& p) {
    if (sizes[p] == -1) stop("Recursion error");
    while (goodChildrens[p].size() == 1) {
      VIDX onlychild = *(goodChildrens[p].begin());
      goodChildrens[p].erase(goodChildrens[p].find(onlychild));
      condenseUp(onlychild, p);
    }
    for (auto it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) condenseTwo(*it);
  }

  void condense(const int& minPts) {
    roots = set< VIDX >();
    goodChildrens = new set< VIDX >[2 * N + 1];
    fallenPointses = new set< VIDX >[2 * N + 1];
    for (VIDX n = 0; n != N; n++) {
      fallenPointses[n] = set< VIDX >();
      fallenPointses[n].insert(n);
    }
    for (VIDX n = N; n != counter; n++) {
      fallenPointses[n] = set< VIDX >();
      goodChildrens[n] = set< VIDX >();
    }
    for (VIDX n = 0; n != counter; ++n) if (p.increment()) {
      condenseOne(n, minPts);
    }
    for (auto it2 = roots.begin();
         it2 != roots.end();
         it2++) {
    	if (p.increment(sizes[*it2])) condenseTwo(*it2);
    	}
  }

  /*
   * Cluster identification phase.
   * Determine stability, starting with roots and recursing.
   * Extract clusters, starting with roots and recursing.
   *
   */
  void determineStability(const VIDX& p, const int& minPts) {
  	if (sizes[p] < minPts && p != parents[p]) stop("Condense failed");
    if (sizes[p] < minPts && p != parents[p]) {
      stop("Condense failed");
    }
    if (goodChildrens[p].size() == 1) stop("Only child");
    double stability = 0;
    double lambda_birth = lambda_births[p];

    for (auto it = fallenPointses[p].begin();
         it != fallenPointses[p].end();
         it++) {
      stability += lambda_births[*it] - lambda_birth;
    }
    VIDX descendantCount = 0;
    for (auto it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) {
      descendantCount += sizes[*it];
      determineStability(*it, minPts);
    }
    stability += descendantCount * (lambda_deaths[p] - lambda_births[p]);
    stabilities[p] = stability;
  }

  void determineStability(int minPts) {
    stabilities = new double[counter];
    selected = new bool[counter];
    for (auto it = roots.begin();
         it != roots.end();
         it++) {
    	if (p.increment(sizes[*it])) determineStability(*it, minPts);
    }
  }

  void deselect(const VIDX& p) {
    selected[p] = false;
    for (auto it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) deselect(*it);
  }

  void extractClusters(const VIDX& p, const int& minPts) {
    survivingClusters.insert(p);
    double childStabilities = 0;
    for (auto it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) {
      extractClusters(*it, minPts);
      childStabilities += stabilities[*it];
    }
    // Cluster has too few points but isn't a root - this should never happen
    // because if it has too few it should be consolidated upwards.
    if (sizes[p] < minPts && parents[p] != p) {
    	stop("child with too few points");
    }
    // If children are superior, don't select, and we're done.
    if (childStabilities > stabilities[p]) { // Children are superior
    		selected[p] = false;
    		stabilities[p] = childStabilities;
    		return;
	  }
    // Cluster is selectable, but apply logic to prevent sole root from ever being a cluster.
    if (parents[p] == p && roots.size() == 1) {
    	selected[p] = false;
    	return;
    }
    selected[p] = true;
  	for (auto it = goodChildrens[p].begin();
     it != goodChildrens[p].end();
     it++) deselect(*it);
  }

  void extractClusters(const VIDX& minPts) {
    survivingClusters = set< VIDX >();
    for (auto it = roots.begin();
         it != roots.end();
         it++) {
    	if (p.increment(sizes[*it])) 	extractClusters(*it, minPts);
    }
  }

  void getClusters(double* ret, const VIDX& n, const int& cluster, const double& cluster_death) const {
  	static int clusterCount = 1;
  	int thisCluster = cluster;
    double thisDeath = cluster_death;
    if (selected[n]) {
#ifdef DEBUG
	Rcout << "\n selecting cluster " << clusterCount << " node " << n;
#endif
      thisCluster = clusterCount++;
      thisDeath = lambda_deaths[n];
    }
    for (auto it = fallenPointses[n].begin();
         it != fallenPointses[n].end();
         it++) {
      ret[*it * 2] = thisCluster;
      ret[(*it * 2) + 1] = min(lambda_births[n], thisDeath) - lambda_births[*it];
    }
    for (auto it = goodChildrens[n].begin();
         it != goodChildrens[n].end();
         it++) getClusters(ret, *it, thisCluster, thisDeath);
  }

  double* getClusters(double* clusters) const {
    for (auto it = roots.begin();
         it != roots.end();
         it++) {
      getClusters(clusters, *it, -1, lambda_deaths[*it]);
    }
    return clusters;
  }

  VIDX reportClusterCnt = 0;

  /*
   * Fetch the hierarchy.
   *
   */
  void reportAHierarchy(const VIDX& last,
                        const VIDX& oldidx,
                        vector<int>& newparent,
                        vector<int>& nodeMembership,
                        vector<double>& newstabilities,
                        vector<int>& newselected,
                        const int& level,
                        const int& childNum) {
    VIDX newidx = reportClusterCnt++;
    if (last == oldidx) newparent[newidx] = newidx;
    else newparent[newidx] = last;

    for (auto it = fallenPointses[oldidx].begin();
         it != fallenPointses[oldidx].end();
         it++) nodeMembership[*it] = newidx;
    newstabilities[newidx] = stabilities[oldidx];
    newselected[newidx] = selected[oldidx];
    int child = 0;
    for (auto it = goodChildrens[oldidx].begin();
         it != goodChildrens[oldidx].end();
         it++) {
      reportAHierarchy(newidx, *it, newparent, nodeMembership,
                       newstabilities, newselected, level + 1, child++);
    }
  }
};

template<class VIDX, class D>
class PrimsAlgorithm {
private:
	const VIDX N;
	PairingHeap<VIDX, D> Q;
	VIDX starterIndex = 0;
	VIDX* minimum_spanning_tree;
	const D* coreDistances;

	void updateVWD(const VIDX& v, const VIDX& w, const double& d) {
		if (!Q.contains(w)) return;
		double dist = max(coreDistances[v], coreDistances[w]);
		dist = max(d, dist);
		if (Q.decreaseIf(w, dist)) minimum_spanning_tree[w] = v;
		//  || w == starterIndex
	}

public:
	PrimsAlgorithm(const VIDX& N,
                 const D* coreDistances) : N{N}, Q(PairingHeap<VIDX,D>(N)), coreDistances{coreDistances} {
		minimum_spanning_tree = new VIDX[N];
	}

	~PrimsAlgorithm() {
		delete[] minimum_spanning_tree;
	}

	VIDX* run(const sp_mat& edges, const IntegerMatrix& neighbors,
           Progress& p, const VIDX& start) {
		starterIndex = start;
		for (VIDX n = 0; n != N; ++n) minimum_spanning_tree[n] = -1;
		Q.batchInsert(N, start);
		Q.decreaseIf(starterIndex, -1);
		p.increment(N);
		VIDX v;
		while (! Q.isEmpty()) {
			if (Q.size() < 0) stop("bad");
			v = Q.pop();
			if (! p.increment()) break;
			if (Q.keyOf(v) == INFINITY || Q.keyOf(v) == -1) starterIndex = v;
			IntegerVector vNeighbors = neighbors.column(v);
			for (auto it = vNeighbors.begin();
        it != vNeighbors.end() && *it != -1;
        it++) {
				updateVWD(v, *it, edges(v, *it));
			}
			for (auto it = edges.begin_col(v);
        it != edges.end_col(v);
        it++) {
				updateVWD(v, it.row(), *it);
			}
		}
		return minimum_spanning_tree;
	}

	vector< pair<D, VIDX> > getMergeSequence() const {
		vector< pair<D, VIDX> > container = vector< pair<D, VIDX> >(N);
		for (VIDX n = 0; n != N; ++n) {
			container[n].second = n;
			container[n].first = Q.keyOf(n);
		}
		sort(container.begin(), container.end());
		return container;
	}
};