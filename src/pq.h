// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
#include "largeVis.h"
#include "minindexedpq.h"

//#define DEBUG
//#define DEBUG2

template<class VIDX>
class UF {
private:
  VIDX reservesize; // used during initialization
	// Used after buildHierarchy to manage condensation, stability extraction, and cluster identificaTION
  arma::Col< VIDX > parents;
  arma::vec lambda_deaths;
  arma::vec stabilities;
  arma::Col< VIDX > sizes;
  std::unique_ptr< std::set< VIDX >[] > fallenPointses;
  std::unique_ptr< std::set< VIDX >[] > goodChildrens;
  std::unique_ptr< bool[] > selected;

protected:
  typedef std::pair<VIDX, double> iddist;
  class CompareDist {
	public:
  	bool operator()(iddist n1, iddist n2) const {
  		return n1.second > n2.second;
  	}
  };

  typedef std::priority_queue<iddist,
                              std::vector< iddist >,
                              CompareDist> DistanceSorter;
  typedef typename std::set<VIDX>::iterator Vidxerator;

  VIDX N;
  Progress p;
  arma::vec lambda_births;

  arma::vec coreDistances;

  VIDX starterIndex = 0;
  PQ<VIDX, double>* Q;
  std::unique_ptr< VIDX[] >   	minimum_spanning_tree;

  VIDX counter = 0;

  std::set< VIDX > survivingClusters;

  int clusterCount = 1;

  // Used after condensation to speed-up remaining phases
  std::set< VIDX > roots;

  UF(VIDX N, bool verbose, bool pq) : N{N}, p(Progress(10 * N, verbose)) {
	  if (pq) Q = new PairingHeap<VIDX, double>(N);
	  else Q = new MinIndexedPQ<VIDX, double>(N);
  	minimum_spanning_tree = std::unique_ptr<VIDX[]>(new VIDX[N]);
#ifdef DEBUG
  	setupTest();
#endif
  }

//  ~UF() {
//    delete Q;
//  }

#ifdef DEBUG
  std::set< VIDX > testers = std::set< VIDX >();
  void setupTest() {
    testers.insert(1);
  	testers.insert(2);
  	testers.insert(3);
  	testers.insert(107);
  	testers.insert(108);
  	testers.insert(208);
  	testers.insert(209);
  	testers.insert(210);
  }
  bool trace(const VIDX p) const {
    bool ret = testers.find(p) != testers.end();
    if (ret) {
      Rcout << "\ntrace " << p << ": ";
      VIDX lastp = p;
      while (lastp != parents[lastp]) {
        Rcout << lastp << " ";
        lastp = parents[lastp];
      }
    }
    return ret;
  }
#endif

  inline double getMRD(const long long i, const long long j, const double dist) const {
  	double d = max(coreDistances[i], coreDistances[j]);
  	d = max(d, dist);
  	return d;
  }

  void updateVWD(VIDX v, VIDX w, double d) {
  	if (Q -> contains(w) || w == starterIndex ) if (Q -> decreaseIf(w, getMRD(v, w, d))) minimum_spanning_tree[w] = v;
  }

  void primsAlgorithm(const arma::sp_mat& edges, VIDX start) {
  	starterIndex = start;
  	for (VIDX n = 0; n != N; n++) minimum_spanning_tree[n] = -1;
  	Q -> batchInsert(N);
  	Q -> decreaseIf(starterIndex, -1);
  	p.increment(N);
  	VIDX v;
  	while (! Q -> isEmpty()) {
  		v = Q -> pop();
  		if (! p.increment()) break;
	  	if (Q -> keyOf(v) == INFINITY || Q -> keyOf(v) == -1) starterIndex = v;
		  for (auto it = edges.begin_row(v);
         it != edges.end_row(v);
         it++) {
		  	updateVWD(v, it.col(), *it);
		  }
		  for (auto it = edges.begin_col(v);
         it != edges.end_col(v);
         it++) {
		  	updateVWD(v, it.row(), *it);
		  }
  	}
  }

  void primsAlgorithm(const arma::sp_mat& edges, const IntegerMatrix& neighbors, VIDX start) {
  	starterIndex = start;
  	for (VIDX n = 0; n != N; n++) minimum_spanning_tree[n] = -1;
  	Q -> batchInsert(N);
  	Q -> decreaseIf(starterIndex, -1);
  	p.increment(N);
  	VIDX v;
  	while (! Q -> isEmpty()) {
  		v = Q -> pop();
  		if (! p.increment()) break;
  		if (Q -> keyOf(v) == INFINITY || Q -> keyOf(v) == -1) starterIndex = v;
  		IntegerVector vNeighbors = neighbors.column(v);
  		for (auto it = vNeighbors.begin();
         	 it != vNeighbors.end() && *it != -1;
         	 it++) {Rcout << ","; updateVWD(v, *it, max(edges(v, *it), edges(*it, v)));}
  	}
  }

  /*
   * The hierarchy construction phase
   */
  VIDX getRoot(VIDX p) const {
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
  void agglomerate(VIDX a, VIDX b, double d) {
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
    lambda_births[n_a] = lambda_births[n_b] = lambda_deaths[parent] = 1 / d;
  }

  void buildHierarchy() {
  	reservesize = 2 * N + 1;
  	parents = arma::Col< VIDX >(reservesize);
  	sizes = arma::Col< VIDX >(reservesize);
  	lambda_births = arma::vec(reservesize);
  	lambda_deaths = arma::vec(reservesize);

  	for (VIDX n = 0; n != N; n++) add();

  	std::vector< iddist > container = std::vector< iddist >();
  	container.reserve(N);
  	typename std::vector< iddist >::iterator adder = container.end();
  	for (VIDX n = 0; n != N; n++) {
  		container.emplace(adder++, n, Q -> keyOf(n));
  		if (n % 50 == 0 ) if (!p.increment(50)) return;
  	}
  	DistanceSorter srtr = DistanceSorter(CompareDist(), container);
  	while (! srtr.empty() && p.increment()) {
  		const iddist ijd = srtr.top();
  		srtr.pop();
  		agglomerate(ijd.first, minimum_spanning_tree[ijd.first], ijd.second);
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
  void condenseUp(VIDX p, VIDX mergeTarget) {
    if (p == mergeTarget) return;
    sizes[p] = -1;
    for (Vidxerator it = fallenPointses[p].begin();
         it != fallenPointses[p].end();
         it++) {
      parents[*it] = mergeTarget;
      fallenPointses[mergeTarget].insert(*it);
    }
    for (Vidxerator it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) {
      goodChildrens[mergeTarget].insert(*it);
      parents[*it] = mergeTarget;
    }
    fallenPointses[p].clear();
    goodChildrens[p].clear();
    lambda_deaths[parents[p]] = max(lambda_deaths[mergeTarget],
                                    lambda_deaths[p]);
    return;
  }

  void condenseOne(VIDX p, int minPts) {
    if (parents[p] ==  p) {
      roots.insert(p);
      return;
    }
    if (sizes[p] < minPts) {
      condenseUp(p, parents[p]);
    } else if (parents[p] != p) goodChildrens[parents[p]].insert(p);
  }

  void condenseTwo(VIDX p) {
    if (sizes[p] == -1) stop("Recursion error");
    while (goodChildrens[p].size() == 1) {
      VIDX onlychild = *(goodChildrens[p].begin());
      goodChildrens[p].erase(goodChildrens[p].find(onlychild));
      condenseUp(onlychild, p);
    }
    for (Vidxerator it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) condenseTwo(*it);
  }

  void condense(int minPts) {
    roots = std::set< VIDX >();
    goodChildrens = std::unique_ptr< std::set< VIDX >[] >(new std::set< VIDX >[2 * N + 1]);
    fallenPointses = std::unique_ptr< std::set< VIDX >[] >(new std::set< VIDX >[2 * N + 1]);
    for (VIDX n = 0; n != N; n++) {
      fallenPointses[n] = std::set< VIDX >();
      fallenPointses[n].insert(n);
    }
    for (VIDX n = N; n != counter; n++) {
      fallenPointses[n] = std::set< VIDX >();
      goodChildrens[n] = std::set< VIDX >();
    }
    for (VIDX n = 0; n != counter; n++) if (p.increment()) {
      condenseOne(n, minPts);
    }
    Vidxerator it2;
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp single private(it2)
{
#endif
    for (it2 = roots.begin();
         it2 != roots.end();
         it2++)
#ifdef _OPENMP
#pragma omp task
#endif
		 {if (p.increment(sizes[*it2])) {condenseTwo(*it2);}}
#ifdef _OPENMP
}
}
#endif
  }

  /*
   * Cluster identification phase.
   * Determine stability, starting with roots and recursing.
   * Extract clusters, starting with roots and recursing.
   *
   */

  void determineStability(VIDX p, int minPts) {
    if (sizes[p] < minPts && p != parents[p]) {
      stop("Condense failed");
    }
    if (goodChildrens[p].size() == 1) stop("Only child");
    double stability = 0;
    double lambda_birth = lambda_births[p];

    for (Vidxerator it = fallenPointses[p].begin();
         it != fallenPointses[p].end();
         it++) {
      stability += lambda_births[*it] - lambda_birth;
    }
    VIDX descendantCount = 0;
    for (Vidxerator it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) {
      descendantCount += sizes[*it];
      determineStability(*it, minPts);
    }
    stability += descendantCount * (lambda_deaths[p] - lambda_births[p]);
    stabilities[p] = stability;
  }

  void determineStability(int minPts) {
    stabilities = arma::Col< double >(counter);
    selected = std::unique_ptr< bool[] >(new bool[counter]);
    Vidxerator it;
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp single private(it)
{
#endif
    for (it = roots.begin();
         it != roots.end();
         it++)
#ifdef _OPENMP
#pragma omp task
#endif
    	{if (p.increment(sizes[*it])) {determineStability(*it, minPts);}}
#ifdef _OPENMP
}}
#endif
  }

  void deselect(VIDX p) {
    selected[p] = false;
    for (Vidxerator it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) deselect(*it);
  }

  void extractClusters(VIDX p) {
#ifdef DEBUG3
    Rcout << "\nextract " << p;
#endif
    survivingClusters.insert(p);
    double childStabilities = 0;
    for (Vidxerator it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) {
      extractClusters(*it);
      childStabilities += stabilities[*it];
    }
    if (childStabilities > stabilities[p]) {
#ifdef DEBUG
    	Rcout << "\npoint " << p << " children are more stable.";
#endif
    	selected[p] = false;
      stabilities[p] = childStabilities;
    } else {
#ifdef DEBUG
Rcout << "\n point " << p << " more stable than children ";
#endif
      // If this is the only root, and there are children, don't let it be selected,
      // because then we'd have only one cluster.
      if (parents[p] == p && roots.size() == 1) {
#ifdef DEBUG
      	Rcout << "not selecting top parent";
#endif
        selected[p] = false;
        return;
      } else {
#ifdef DEBUG
      	Rcout << "selecting, deselcting children.";
#endif
      	selected[p] = true;
      	for (Vidxerator it = goodChildrens[p].begin();
            it != goodChildrens[p].end();
            it++) deselect(*it);
      	}
    }
  }

  void extractClusters() {
    survivingClusters = std::set< VIDX >();
  	Vidxerator it;
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp single private(it)
{
#endif
    for (it = roots.begin();
         it != roots.end();
         it++)
#ifdef _OPENMP
#pragma omp task
#endif
    	   {if (p.increment(sizes[*it])) 	{extractClusters(*it);}}
#ifdef _OPENMP
}}
#endif
  }

  void getClusters(arma::mat& ret, VIDX n, int cluster, double cluster_death) {
    int thisCluster = cluster;
    double thisDeath = cluster_death;
    if (selected[n])
#ifdef _OPENMP
#pragma omp critical
#endif
		{
#ifdef DEBUG
	Rcout << "\n selecting cluster " << clusterCount << " node " << n;
#endif
      thisCluster = clusterCount++;
      thisDeath = lambda_deaths[n];
    }
    for (Vidxerator it = fallenPointses[n].begin();
         it != fallenPointses[n].end();
         it++) {
      ret(0, *it) = thisCluster;
      ret(1, *it) = min(lambda_births[n], thisDeath) - lambda_births[*it];
    }
    for (Vidxerator it = goodChildrens[n].begin();
         it != goodChildrens[n].end();
         it++) getClusters(ret, *it, thisCluster, thisDeath);
  }

  arma::mat getClusters() {
    arma::mat ret = arma::mat(2, N, fill::zeros);
  	Vidxerator it;
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp single private(it)
{
#endif
    for (it = roots.begin();
         it != roots.end();
         it++)
#ifdef _OPENMP
#pragma omp task
#endif
    	{
      getClusters(ret, *it, -1, lambda_deaths[*it]);
    }
#ifdef _OPENMP
}}
#endif
    return ret;
  }

  /*
   * Fetch the hierarchy.
   *
   */

  VIDX reportClusterCnt = 0;

  void reportAHierarchy(VIDX last, VIDX oldidx,
                       IntegerVector& newparent,
                       IntegerVector& nodeMembership,
                       NumericVector& newstabilities,
                       IntegerVector& newselected,
                       int level,
                       int childNum) {
#ifdef DEBUG
    Rcout << "\n";
    for (int i = 0; i < level; i++) Rcout << " |";
    if (level != 0 && childNum == 0) Rcout << "\\_";
    else if (level != 0) Rcout << "+-";
    Rcout << oldidx;
    if (level == 0) Rcout << "\t";
    Rcout << "\tstability: " << stabilities[oldidx] <<
      "  selected: " << selected[oldidx] << "  fallen: " <<
        fallenPointses[oldidx].size() << " sz: " << sizes[oldidx];
    Rcout << " lambda " << lambda_births[oldidx] << "/" << lambda_deaths[oldidx];
#endif
    VIDX newidx = reportClusterCnt++;
    if (last == oldidx) newparent[newidx] = newidx;
    else newparent[newidx] = last;

    for (Vidxerator it = fallenPointses[oldidx].begin();
         it != fallenPointses[oldidx].end();
         it++) nodeMembership[*it] = newidx;
    newstabilities[newidx] = stabilities[oldidx];
    newselected[newidx] = selected[oldidx];
    int child = 0;
    for (Vidxerator it = goodChildrens[oldidx].begin();
         it != goodChildrens[oldidx].end();
         it++) {
      reportAHierarchy(newidx, *it, newparent, nodeMembership,
                       newstabilities, newselected, level + 1, child++);
    }
  }
};