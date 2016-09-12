// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include <Rmath.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

typedef pair<long long, double> iddist;

class CompareDist {
public:
  bool operator()(iddist n1, iddist n2) {
    return n1.second > n2.second;
  }
};

typedef std::priority_queue<iddist,
                            vector<iddist>,
                            CompareDist> NNheap;
typedef std::vector<iddist> NNlist;

class Cluster {
protected:
	const long long N;
	int K;
	const double radius2;
	const int minPts;
	int D;
  arma::mat* data = NULL;
  arma::imat* neighbors = NULL;
  const arma::sp_mat* edges;
  const bool hasEdges;
  bool hasData;

  std::vector<bool> visited;

  Progress progress;

  inline double dist2(const double * x_i,
                      long long id) const {
    double dist = 0;
    const double * x_j = data->colptr(id);
    for (int i = 0; i != D; i++) dist += (x_i[i] - x_j[i]) * (x_i[i] - x_j[i]);
    return dist;
  }

  Cluster(arma::imat& neighbors,
          double eps,
          int minPts,
          bool verbose) : N(neighbors.n_cols),
          								K(neighbors.n_rows),
          								radius2{eps * eps},
          								edges(nullptr),
          								minPts{minPts},
                          neighbors{&neighbors},
                        	hasEdges(false),
                        	hasData(false),
                        	visited(std::vector<bool>(N, false)),
                        	progress(Progress(N, verbose)) {}

  Cluster(arma::sp_mat& edges,
          double eps,
          int minPts,
          bool verbose) : N(edges.n_cols),
          								radius2{eps * eps}, minPts{minPts},
          								edges{&edges},
                          hasEdges(true),
                          hasData(false),
                          visited(std::vector<bool>(N, false)),
                          progress(Progress(N, verbose))
                          { }

  void frNNrecurse(double distStart,
                   const double * x_i,
                   long long id,
                   NNheap& found,
                   std::set<long long>& checked) const {
    double dist;
    for (auto it = neighbors ->begin_col(id);
         it != neighbors->end_col(id) && *it != -1;
         it++) {
      if (checked.insert(*it).second) {
        dist = dist2(x_i, *it);
        if (dist < radius2) {
          found.push(iddist(*it, dist));
          frNNrecurse(distStart + dist, x_i, *it, found, checked);
        }
      }
    }
  }

  void frNNrecurseEdge(double distStart,
                        const double * x_i,
                        long long id,
                        NNheap& found,
                        std::set<long long>& checked) const {
    double dist;
    for (auto it = edges -> begin_col(id);
         it != edges -> end_col(id);
         it++) {
      if (checked.insert(it.row()).second && distStart + *it < radius2) {
        dist = dist2(x_i, it.row());
        if (dist < radius2) {
          found.push(iddist(it.row(), dist));
          frNNrecurseEdge(distStart + dist,
                      x_i,
                      it.row(),
                      found,
                      checked);
        }
      }
    }
  }

  NNlist fixedRadiusNearestNeighbors(long long id) const {
    NNheap found = NNheap();
    std::set<long long> checked = std::set<long long>();

    found.push(iddist(id, 0));
    checked.insert(id);

    const double * x_i = (hasData) ? data->colptr(id) : NULL;
    if (hasEdges) {
      auto end = edges -> end_col(id);
      for (auto it = edges -> begin_col(id); it != end; it++) {
        if (checked.insert(it.row()).second) {
          found.push(iddist(it.row(), *it));
          if (hasData)frNNrecurseEdge(0, x_i, it.row(), found, checked);
        }
      }
    } else for (auto it = neighbors->begin_col(id);
         it != neighbors->end_col(id) && *it != -1;
         it++) {
      frNNrecurse(0, x_i, *it, found, checked);
    }

    NNlist ret = NNlist();
    while (! found.empty()) {
      ret.push_back(found.top());
      found.pop();
    }
    return ret;
  }

public:
  void setData(arma::mat& data) {
    this -> data = &data;
    D = data.n_rows;
    hasData = true;
  }
  virtual IntegerVector run() = 0;
};

class OPTICS : public Cluster {
protected:
  std::vector<long long> orderedPoints, seeds;
  std::vector<long double> ds;
public:
  std::vector<long double> reachdist, coredist;
  OPTICS( arma::imat& neighbors,
          double eps,
          int minPts,
          bool verbose) : Cluster(neighbors, eps, minPts, verbose),
								          orderedPoints(std::vector<long long>()),
								          seeds(std::vector<long long>()),
          								ds(std::vector<long double>()),
          								reachdist(std::vector<long double>(N, INFINITY)),
								          coredist(std::vector<long double>(N, INFINITY)) {
    orderedPoints.reserve(N);
  }

  OPTICS( arma::SpMat<double>& edges,
          double eps,
          int minPts,
          bool verbose) : Cluster(edges, eps, minPts, verbose),
								          orderedPoints(std::vector<long long>()),
								          seeds(std::vector<long long>()),
								          ds(std::vector<long double>()),
								          reachdist(std::vector<long double>(N, INFINITY)),
								          coredist(std::vector<long double>(N, INFINITY)) {
    orderedPoints.reserve(N);
  }

  void update(NNlist& frNeighbors,
              long long p) {

    std::vector<long long>::iterator pos_seeds;
    long double newreachdist;
    long long o;
    long double o_d;

    while(!frNeighbors.empty()) {
      o = frNeighbors.back().first;
      o_d = frNeighbors.back().second;
      frNeighbors.pop_back();

      if(visited[o]) continue;

      newreachdist = max(coredist[p], o_d);

      if(reachdist[o] == INFINITY) {
        reachdist[o] = newreachdist;
        seeds.push_back(o);
      } else if(newreachdist < reachdist[o]) reachdist[o] = newreachdist;
    }
  }

  IntegerVector run() {
    if (hasData + hasEdges < 1) stop("Need either data or edges.");
    NNlist frNeighbors;
    for (long long n = 0; n < N; n++) if (progress.increment()) {
      if (visited[n]) continue;

      frNeighbors = fixedRadiusNearestNeighbors(n);
      visited[n] = true;

      // find core distance
      if(frNeighbors.size() >= (size_t) minPts)
        coredist[n] = frNeighbors[minPts-1].second;

      orderedPoints.push_back(n);

      if (coredist[n] == INFINITY) continue; // core-dist is undefined

      // update
      update(frNeighbors, n);

      long long q;
      while (!seeds.empty()) {
        // get smallest dist (to emulate priority queue). All should have already
        // a reachability distance <Inf from update().
        std::vector<long long>::iterator q_it = seeds.begin();
        for (std::vector<long long>::iterator it = seeds.begin();
             it!=seeds.end(); ++it) if (reachdist[*it] < reachdist[*q_it] ||
                                        (reachdist[*it] == reachdist[*q_it] &&
                                              *q_it < *it)) q_it = it;
        q = *q_it;
        seeds.erase(q_it);

        frNeighbors = fixedRadiusNearestNeighbors(q);
        visited[q] = true;

        if(frNeighbors.size() >= (size_t) minPts) {
          coredist[q] = frNeighbors[minPts-1].second;
        }

        orderedPoints.push_back(q);

        if(frNeighbors.size() < (size_t) minPts) continue; //  == q has no core dist.

        // update seeds
        update(frNeighbors, q);
      }
    }
    return IntegerVector(orderedPoints.begin(), orderedPoints.end())+1;
  }
};

List optics_assemble(OPTICS& opt) {
  List ret;
  IntegerVector vec = opt.run();
  ret["order"] = vec;
  ret["reachdist"] = sqrt(NumericVector(opt.reachdist.begin(), opt.reachdist.end()));
  ret["coredist"] = sqrt(NumericVector(opt.coredist.begin(), opt.coredist.end()));

  return ret;
}

// [[Rcpp::export]]
List optics_e(arma::sp_mat& edges,
                double eps,
                int minPts,
                bool verbose) {
  OPTICS opt = OPTICS(edges, eps, minPts, verbose);
  return optics_assemble(opt);
}

// [[Rcpp::export]]
List optics_ed(arma::sp_mat& edges,
                arma::mat& data,
                double eps,
                int minPts,
                bool verbose) {
  OPTICS opt = OPTICS(edges, eps, minPts, verbose);
  opt.setData(data);
  return optics_assemble(opt);
}

// [[Rcpp::export]]
List optics_nd(arma::imat& neighbors,
               arma::mat& data,
               double eps,
               int minPts,
               bool verbose) {
  OPTICS opt = OPTICS(neighbors, eps, minPts, verbose);
  opt.setData(data);
  return optics_assemble(opt);
}

