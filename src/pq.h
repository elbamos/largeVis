// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
#include "largeVis.h"
//#define DEBUG
//#define DEBUG2
// MinIndexedPQ Basedon https://github.com/kartikkukreja/blog-codes/blob/master/src/Indexed%20Min%20Priority%20Queue.cpp

template<class VIDX, class D>
class MinIndexedPQ {
private:
  VIDX NMAX, N;
  std::unique_ptr< VIDX[] > heap, index;
  std::unique_ptr< D[] > keys;

  void swap(VIDX i, VIDX j) {
    VIDX t = heap[i];
    heap[i] = heap[j];
    heap[j] = t;
    index[heap[i]] = i;
    index[heap[j]] = j;
  }

  void bubbleUp(VIDX k)    {
    while(k > 1 && keys[heap[k/2]] > keys[heap[k]])   {
      swap(k, k/2);
      k = k/2;
    }
  }

  void bubbleDown(VIDX k)  {
    VIDX j;
    while(2*k <= N) {
      j = 2*k;
      if(j < N && keys[heap[j]] > keys[heap[j+1]])
        j++;
      if(keys[heap[k]] <= keys[heap[j]])
        break;
      swap(k, j);
      k = j;
    }
  }

public:
  // Create an empty MinIndexedPQ which can contain atmost NMAX elements
  MinIndexedPQ(VIDX NMAX) : NMAX{NMAX},
                                     N(0) {
    keys = std::unique_ptr< D[] >(new D[NMAX + 1]);
    heap = std::unique_ptr< VIDX[] >(new VIDX[NMAX + 1]);
    index = std::unique_ptr< VIDX[] >(new VIDX[NMAX + 1]);
    for(VIDX n = 0; n != NMAX; n++) index[n] = -1;
  }

  // check if the PQ is empty
  bool isEmpty()  {
    return N == 0;
  }

  // check if i is an index on the PQ
  bool contains(VIDX i)    {
    return index[i] != -1;
  }

  // return the number of elements in the PQ
  VIDX size()  {
    return N;
  }

  // associate key with index i; 0 < i < NMAX
  void insert(VIDX i, D key) {
    N++;
    index[i] = N;
    heap[N] = i;
    keys[i] = key;
    bubbleUp(N);
  }

  // returns the index associated with the minimal key
  VIDX minIndex()  {
    return heap[1];
  }

  // returns the minimal key
  D minKey()    {
    return keys[heap[1]];
  }

  // delete the minimal key and return its associated index
  // Warning: Don't try to read from this index after calling this function
  VIDX deleteMin() {
    VIDX min = heap[1];
    swap(1, N--);
    bubbleDown(1);
    index[min] = -1;
    heap[N+1] = -1;
    return min;
  }

  // returns the key associated with index i
  D keyOf(VIDX i)    {
    return keys[i];
  }

  // change the key associated with index i to the specified value
  void changeKey(VIDX i, D key)  {
    keys[i] = key;
    bubbleUp(index[i]);
    bubbleDown(index[i]);
  }

  // delete the key associated with index i
  void remove(VIDX i)   {
    VIDX ind = index[i];
    swap(ind, N--);
    bubbleUp(ind);
    bubbleDown(ind);
    index[i] = -1;
  }
};

template<class VIDX>
class UF {
private:
  VIDX reservesize;
  arma::Col< VIDX > parents;
  arma::Col< double > lambda_deaths;
  arma::Col< double > stabilities;
  arma::Col< VIDX > sizes;
  std::unique_ptr< std::set< VIDX >[] > fallenPointses;
  std::unique_ptr< std::set< VIDX >[] > goodChildrens;
  arma::Col< int > selected;
protected:
  VIDX N;
  typedef std::pair<VIDX, double> iddist;
  Progress p;
  arma::Col< double > lambda_births;

  UF(VIDX N, VIDX nedges, bool verbose) : N{N},
                                          p(Progress((9 * N) + (nedges), verbose)) {

  }

#ifdef DEBUG
  std::set< VIDX > testers = std::set< VIDX >();
  void setupTest() {
    testers.insert(209);
    testers.insert(234);
    testers.insert(1541);
    testers.insert(1546);
    testers.insert(1548);
    testers.insert(1553);
    testers.insert(1560);
    testers.insert(1569);
  }
  bool trace(const VIDX p) {
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

  void setup() {
#ifdef DEBUG
    setupTest();
#endif
    reservesize = 2 * N + 1;
    parents = arma::Col< VIDX >(reservesize);
    sizes = arma::Col< VIDX >(reservesize);
    lambda_births = arma::Col< double >(reservesize);
    lambda_deaths = arma::Col< double >(reservesize);
    for (VIDX n = 0; n != N; n++) add();
  }

  VIDX getRoot(VIDX p) {
    return (parents[p] == p) ? p : getRoot(parents[p]);
  }

  VIDX counter = 0;

  VIDX add() {
    VIDX newid = counter++;
    parents[newid] = newid;
    sizes[newid] = 1;
    lambda_births[newid] = 0;
    lambda_deaths[newid] = 0;
    return newid;
  }

  void agglomerate(VIDX a, VIDX b, double d) {
    if (b == -1) {
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

  std::set< VIDX > roots;

  void condenseUp(VIDX p, VIDX mergeTarget) {
    if (p == mergeTarget) return;
    sizes[p] = -1;
    for (typename std::set<VIDX>::iterator it = fallenPointses[p].begin();
         it != fallenPointses[p].end();
         it++) {
      parents[*it] = mergeTarget;
      fallenPointses[mergeTarget].insert(*it);
    }
    for (typename std::set<VIDX>::iterator it = goodChildrens[p].begin();
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
    for (typename std::set< VIDX >::iterator it = goodChildrens[p].begin();
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
    for (typename std::set< VIDX >::iterator it = roots.begin();
         it != roots.end();
         it++) condenseTwo(*it);
  }

  void determineStability(VIDX p, int minPts) {
    if (sizes[p] < minPts && p != parents[p]) {
      stop("Condense failed");
    }
    if (goodChildrens[p].size() == 1) stop("Only child");
    double stability = 0;
    double lambda_birth = lambda_births[p];

    for (typename std::set< VIDX >::iterator it = fallenPointses[p].begin();
         it != fallenPointses[p].end();
         it++) {
      stability += lambda_births[*it] - lambda_birth;
    }
    VIDX descendantCount = 0;
    for (typename std::set< VIDX >::iterator it = goodChildrens[p].begin();
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
    selected = arma::Col< int >(counter);

    for (typename std::set< VIDX >::iterator it = roots.begin();
         it != roots.end();
         it++) determineStability(*it, minPts);
  }

  void deselect(VIDX p) {
    selected[p] = false;
    for (typename std::set< VIDX >::iterator it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) deselect(*it);
  }

  std::set< VIDX > survivingClusters;

  void extractClusters(VIDX p) {
    survivingClusters.insert(p);
    double childStabilities = 0;
    for (typename std::set< VIDX >::iterator it = goodChildrens[p].begin();
         it != goodChildrens[p].end();
         it++) {
      extractClusters(*it);
      childStabilities += stabilities[*it];
    }
    if (childStabilities > stabilities[p]) {
      stabilities[p] = childStabilities;
    } else {
      // If this is the only root, and there are children, don't let it be selected, 
      // because then we'd have only one cluster.
      if (parents[p] == p && roots.size() == 1 && goodChildrens[p].size() > 0) {
        selected[p] = false;
        return;
      }
      selected[p] = true;
      for (typename std::set< VIDX >::iterator it = goodChildrens[p].begin();
           it != goodChildrens[p].end();
           it++) deselect(*it);
    }
  }

  void extractClusters() {
    survivingClusters = std::set< VIDX >();
    for (typename std::set< VIDX >::iterator it = roots.begin();
         it != roots.end();
         it++) extractClusters(*it);
  }

  int clusterCount = 0;

  void getClusters(arma::mat& ret, VIDX n, int cluster, double cluster_death) {
    int thisCluster = cluster;
    double thisDeath = cluster_death;
    if (selected[n]) {
      thisCluster = clusterCount++;
      thisDeath = lambda_deaths[n];
    }
    for (typename std::set< VIDX >::iterator it = fallenPointses[n].begin();
         it != fallenPointses[n].end();
         it++) {
      ret(0, *it) = thisCluster;
      ret(1, *it) = min(lambda_births[n], thisDeath) - lambda_births[*it];
    }
    for (typename std::set< VIDX >::iterator it = goodChildrens[n].begin();
         it != goodChildrens[n].end();
         it++) getClusters(ret, *it, thisCluster, thisDeath);
  }

  arma::mat getClusters() {
    arma::mat ret = arma::mat(2, N, fill::zeros);
    for (typename std::set< VIDX >::iterator it = roots.begin();
         it != roots.end();
         it++) {
      getClusters(ret, *it, -1, lambda_deaths[*it]);
    }
    return ret;
  }

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
    newparent[newidx] = last;
    for (typename std::set< VIDX >::iterator it = fallenPointses[oldidx].begin();
         it != fallenPointses[oldidx].end();
         it++) nodeMembership[*it] = newidx;
    newstabilities[newidx] = stabilities[oldidx];
    newselected[newidx] = selected[oldidx];
    int child = 0;
    for (typename std::set< VIDX >::iterator it = goodChildrens[oldidx].begin();
         it != goodChildrens[oldidx].end();
         it++) {
      reportAHierarchy(newidx, *it, newparent, nodeMembership,
                       newstabilities, newselected, level + 1, child++);
    }
  }

};