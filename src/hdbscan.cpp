// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "pq.h"

class HDBSCAN {
private:
  typedef std::pair<long long, double> iddist;
  class CompareDist {
  public:
    bool operator()(iddist n1, iddist n2) {
      return n1.second > n2.second;
    }
  };
  typedef std::priority_queue<iddist, 
                              std::vector<iddist>, 
                              CompareDist> DistanceSorter;
  
  typedef std::pair<std::pair<long long, long long>,double> pairdist;
  class ComparePairDist {
  public:
    bool operator()(pairdist n1, pairdist n2) {
      return n1.second > n2.second;
    }
  };
  typedef std::priority_queue<pairdist, 
                              std::vector<pairdist>, 
                              ComparePairDist> PairDistanceSorter;
  
  const int                 K;
  const long long           N;
  
  long long starterIndex;
  
  arma::vec makeCoreDistances(const arma::sp_mat& edges) {
   arma::vec coreDistances = arma::vec(N);
   double bestOverall = INFINITY;
    for (long long n = 0; n < N; n++) {
      DistanceSorter srtr = DistanceSorter();
      for (auto it = edges.begin_col(n);
           it != edges.end_col(n);
           it++) srtr.emplace(iddist(it.row(), *it));
     double best = INFINITY;
      for (int i = 0; i != K && ! srtr.empty(); i++) {
        best = srtr.top().second;
        if (best < bestOverall) {
          bestOverall = best;
          starterIndex = n;
        }
        srtr.pop();
      }
      coreDistances[n] = best;
    }
    return coreDistances;
  }
  
  arma::sp_mat mrdMatrix;
  
  void makeMRD(const arma::sp_mat& edges) {
    mrdMatrix = arma::sp_mat(edges);
    arma::vec coreDistances = makeCoreDistances(edges);
    for (auto it = mrdMatrix.begin(); 
         it != mrdMatrix.end();
         it++) {
      long long i = it.row();
      long long j = it.col();
      *it = (coreDistances[i] > coreDistances[j]) ? coreDistances[i] : coreDistances[j]; 
      *it = (*it > edges(i, j)) ? *it : edges(i, j);
    }
  }
  
  long long*   minimum_spanning_tree;
  double*      minimum_spanning_distances;
  
  void primsAlgorithm() {
    double* Cv = minimum_spanning_distances;
    long long* Ev = minimum_spanning_tree;
    MinIndexedPQ<long long, double> Q = MinIndexedPQ<long long, double>(N);
    for (long long n = 0; n != N; n++) {
      Cv[n] = (n == starterIndex) ? -1 : INFINITY;
      Ev[n] = -1;
      Q.insert(n, Cv[n]);
    }

    long long v;
    while (! Q.isEmpty()) {
      v = Q.deleteMin();
      for (auto it = mrdMatrix.begin_col(v);
           it != mrdMatrix.end_col(v);
           it++) {
        long long w = it.row();
        if (Ev[w] == -1) {
          Cv[w] = *it;
          Ev[w] = v;
        }
        if (Q.contains(w) && *it < Cv[w]) {
          Q.changeKey(w, *it);
          Cv[w] = *it;
          Ev[w] = v;
        }
      }
    }
    for (int n = 0; n < N; n++) Rcout << Cv[n] << " ";
    Rcout << "\n";
    for (int n = 0; n < N; n++) Rcout << Ev[n] << " ";
  }
  UF<long long>* hierarchy = nullptr;
  void buildHierarchy(int minPts) {
    DistanceSorter srtr = DistanceSorter();
    for (long long n = 0; n != N; n++) {
      // FIX ME - This should never be used.
      double distance = (minimum_spanning_distances[n] == -1) ?
                        minimum_spanning_distances[minimum_spanning_tree[n]] :
                        minimum_spanning_distances[n];
      srtr.push(iddist(n, distance));
    }
    hierarchy = new UF<long long>(N);
    while (! srtr.empty()) {
      const iddist ijd = srtr.top();
      srtr.pop();
      hierarchy -> agglomerate(ijd.first, minimum_spanning_tree[ijd.first], ijd.second);
    }
    hierarchy -> condense(minPts);
  }
  
    
  
public:
  HDBSCAN(const int N,
          const int K) : K{K}, N{N} {
    minimum_spanning_tree = new long long[N];
    minimum_spanning_distances = new double[N];
  }
  
  arma::mat process(const arma::sp_mat& edges, int minPts) {
    makeMRD(edges);
    primsAlgorithm();
    buildHierarchy(minPts);
    hierarchy -> determineStability();
    return hierarchy -> extractClusters(minPts);
  }
};

// [[Rcpp::export]]
arma::mat hdbscan(const arma::sp_mat& edges, int K, int minPts) {
  HDBSCAN object = HDBSCAN(edges.n_cols, K);
  return object.process(edges, minPts);
}

template<class VIDX>
long long UF_node<VIDX>::counter = 0;