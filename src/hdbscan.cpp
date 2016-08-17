// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "pq.h"

class HDBSCAN {
public:
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
  
  typedef std::pair<std::pair<long long, long long>, double> pairdist;
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
    for (long long n = 0; n < N; n++) {
      DistanceSorter srtr = DistanceSorter();
      for (auto it = edges.begin_col(n);
           it != edges.end_col(n);
           it++) srtr.emplace(iddist(it.row(), *it));
      for (int k = 0; k != K && srtr.size() > 1; k++) srtr.pop();
      coreDistances[n] = srtr.top().second;
    }
    return coreDistances;
  }
  
  arma::sp_mat mrdMatrix;
  
  void makeMRD(const arma::sp_mat& edges) {
    mrdMatrix = arma::sp_mat(edges);
    arma::vec coreDistances = makeCoreDistances(edges);
    double bestOverall = INFINITY;
    for (auto it = mrdMatrix.begin(); 
         it != mrdMatrix.end();
         it++) {
      long long i = it.row();
      long long j = it.col();
      double d = max(coreDistances[i], coreDistances[j]);
      d = (d > *it) ? d : *it;
      if (d < bestOverall) {
        starterIndex = j;
        bestOverall = d;
      }
      *it = d;
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
      Rcout << "\n" << Q.minKey();
      v = Q.deleteMin();
      Rcout << " " << v;
      for (auto it = mrdMatrix.begin_row(v);
           it != mrdMatrix.end_row(v);
           it++) {
        long long w = it.col();
        Rcout << " " << w;
      //  if (Ev[w] == -1) {
      //    Cv[w] = *it;
      //    Ev[w] = v;
      //  }
        if (Q.contains(w) && *it < Cv[w]) {
          Rcout << "a";
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
  void buildHierarchy() {
    DistanceSorter srtr = DistanceSorter();
    for (long long n = 0; n != N; n++) {
      // FIX ME - Handle the roots from the prior phase
      double distance = minimum_spanning_distances[n];
      srtr.push(iddist(n, distance));
    }
    hierarchy = new UF<long long>(N);
    while (! srtr.empty()) {
      const iddist ijd = srtr.top();
      srtr.pop();
      hierarchy -> agglomerate(ijd.first, minimum_spanning_tree[ijd.first], ijd.second);
    }
  }
  
  HDBSCAN(const int N,
          const int K) : K{K}, N{N} {
            minimum_spanning_tree = new long long[N];
            minimum_spanning_distances = new double[N];
          }
  
  arma::mat process(const arma::sp_mat& edges, int minPts) {
    makeMRD(edges);
    primsAlgorithm();
    buildHierarchy();
    hierarchy -> condense(minPts);
    hierarchy -> determineStability(minPts);
    hierarchy -> extractClusters();
    return hierarchy -> getClusters();
  }
};

// [[Rcpp::export]]
List hdbscanc(const arma::sp_mat& edges, int K, int minPts) {
  HDBSCAN object = HDBSCAN(edges.n_cols, K);
  arma::mat clusters = object.process(edges, minPts);
  arma::ivec tree = arma::ivec(edges.n_cols);
  for (int n = 0; n != edges.n_cols; n++) {
    tree[n] = object.minimum_spanning_tree[n];
  }
  return Rcpp::List::create(Rcpp::Named("clusters") = clusters,
                            Rcpp::Named("tree") = tree);
}