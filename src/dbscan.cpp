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

/*
 * The code below is based on code which bears the following copyright notice; portions may
 * remain subject to it:
#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

class Cluster {
protected:
  int D;
  long long N;
  int K;
  double radius2;
  int minPts;

  arma::mat* data = NULL;
  arma::imat* neighbors = NULL;
  arma::sp_mat* edges = NULL;
  bool hasEdges;
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
          bool verbose) :
                        neighbors{&neighbors},
                        minPts{minPts},
                        radius2{eps * eps},
                        N(neighbors.n_cols),
                        K(neighbors.n_rows),
                        hasEdges(false),
                        hasData(false),
                        progress(Progress(N, verbose)),
                        visited(std::vector<bool>(N, false)) {}

  Cluster(arma::sp_mat& edges,
          double eps,
          int minPts,
          bool verbose) :
                          edges{&edges},
                          hasEdges(true),
                          hasData(false),
                          minPts{minPts},
                          radius2{eps * eps},
                          N(edges.n_cols),
                          progress(Progress(N, verbose)),
                          visited(std::vector<bool>(N, false)) { }

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

class DBSCAN : public Cluster {
public:
  DBSCAN( arma::imat& neighbors,
          double eps,
          int minPts,
          bool verbose) : Cluster(neighbors, eps, minPts, verbose) {}

  DBSCAN( arma::SpMat<double>& edges,
          double eps,
          int minPts,
          bool verbose) : Cluster(edges, eps, minPts, verbose) {}

  IntegerVector run() {
    if (hasData + hasEdges < 1) stop("Need either data or edges.");
    std::vector< std::vector<long long> > clusters; // vector of vectors == list
    NNlist  frNN1, frNN2;

    for (long long n = 0; n != N; n++) if (progress.increment()) {
      if (visited[n]) continue;

      // start new cluster and expand
      frNN1 = fixedRadiusNearestNeighbors(n);
      if (frNN1.size() < minPts) continue;
      std::vector<long long> cluster;
      cluster.push_back(n);
      visited[n] = true;

      while (!frNN1.empty()) {
        long long j = frNN1.back().first;
        frNN1.pop_back();

        if (visited[j]) continue; // point already processed
        visited[j] = true;

        frNN2 = fixedRadiusNearestNeighbors(j);

        if (frNN2.size() >= minPts) { // expand neighborhood
          copy(frNN2.begin(), frNN2.end(),
               back_inserter(frNN1));
          cluster.push_back(j);
        }
      }
      clusters.push_back(cluster);
    }
    IntegerVector ret = IntegerVector(N,0);
    for (int c = 0; c < clusters.size(); c++) {
      std::vector<long long> cluster = clusters[c];
      for (int j=0; j < cluster.size(); j++) ret[cluster[j]] = c + 1 ;
    }
    return ret;
  }
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
                          reachdist(std::vector<long double>(N, INFINITY)),
                          coredist(std::vector<long double>(N, INFINITY)),
                          ds(std::vector<long double>()),
                          orderedPoints(std::vector<long long>()),
                          seeds(std::vector<long long>()) {
    orderedPoints.reserve(N);
  }

  OPTICS( arma::SpMat<double>& edges,
          double eps,
          int minPts,
          bool verbose) : Cluster(edges, eps, minPts, verbose),
          reachdist(std::vector<long double>(N, INFINITY)),
          coredist(std::vector<long double>(N, INFINITY)),
          ds(std::vector<long double>()),
          orderedPoints(std::vector<long long>()),
          seeds(std::vector<long long>()) {
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

      newreachdist = std::max(coredist[p], o_d);

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

        //N2 = regionQueryDist(q, dataPts, kdTree, eps2, approx);
        frNeighbors = fixedRadiusNearestNeighbors(q);
        visited[q] = true;

        // update core distance
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

// [[Rcpp::export]]
IntegerVector dbscan_e(arma::sp_mat& edges,
                         double eps,
                         int minPts,
                         bool verbose) {
  DBSCAN db = DBSCAN(edges, eps, minPts, verbose);
  return db.run();
}

// [[Rcpp::export]]
IntegerVector dbscan_ed(arma::sp_mat& edges,
                         arma::mat& data,
                         double eps,
                         int minPts,
                         bool verbose) {
  DBSCAN db = DBSCAN(edges, eps, minPts, verbose);
  db.setData(data);
  return db.run();
}

// [[Rcpp::export]]
IntegerVector dbscan_nd(arma::imat& neighbors,
                         arma::mat& data,
                         double eps,
                         int minPts,
                         bool verbose) {
  DBSCAN db = DBSCAN(neighbors, eps, minPts, verbose);
  db.setData(data);
  return db.run();
}

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

// [[Rcpp::export]]
void silhouetteDbscan(const arma::sp_mat& edges,
                      NumericMatrix sil) {
  long long N = edges.n_cols;
  const NumericVector clusters = sil.column(0);
  int K = max(clusters) + 1; // n clusters
  long long* counts = new long long[K];
  for (int k = 0; k != K; k++) counts[k] = 0;
  double* diC = new double[K * N];
  for (long long kn = 0; kn != K * N; kn++) diC[kn] = 0;

  /*
   * Loop through points, accumulating for each cluster the total distance of the point
   * to all points in that cluster.  Distances to points in the noise cluster are ignored.
   * (Distances from points in the noise cluster, are tracked.)
   *
   * Using a kNN matrix to start, we don't have all N^2 distances. When a distance is missing,
   * we use 2 * the distance from that point to its furthest known nearest neighbor.
   */
  for (long long n = 0; n != N; n++) {
    long long last = 0;
    auto end = edges.end_col(n);
    int ns_cluster = (int) clusters[n];
    counts[ns_cluster]++;
    double* diN = diC + (K * n);
    double maxD = edges.col(n).max();
    for (auto it = edges.begin_col(n);
         it != end;
         it++) {
      while (last != it.row()) {
        if (clusters[last] != 0 || ns_cluster == 0) {
          double d = edges(n, it.row());
          d = (d == 0) ? maxD : d;
          diN[(int) clusters[last]] += d;
        }
        last++;
      }
      if (it.row() == n) continue;
      if (clusters[last] != 0 || ns_cluster == 0) diN[(int) clusters[last]] += *it;
    }
    while (last != N) {
      if (clusters[last] != 0 || ns_cluster == 0) diN[(int) clusters[last]] += maxD;
      last++;
    }
  }

  for (long long n = 0; n != N; n++) {
    const int ns_cluster = (int) clusters[n];
    bool computeSn = TRUE;
    double* diN = diC + (K * n);

    for (long long k = 0; k != K; k++) {
      if (k == ns_cluster) {
        if (counts[k] == 1) computeSn = FALSE;
        else diN[k] = diN[k] / (counts[k] - 1);
      } else if (ns_cluster != 0) diN[k] = diN[k] / counts[k];
    }
    const double an = diN[ns_cluster];
    /*
     * Iterate through clusters to find nearest cluster-neighbor.
     * Start by assuming neighbor of cluster 1 is 2, and of all other clusters is 1.
     */
    int candidateNeighbor = (ns_cluster == 1) ? 2 : 1;
    double bn = diN[candidateNeighbor];
    for (int k = 1; k != K; k++) if (k != ns_cluster && diN[k] < bn) {
      bn = diN[k];
      candidateNeighbor = k;
    }
    sil(n, 1) = candidateNeighbor;
    if (counts[ns_cluster] == 1) bn = 0;
    else bn = (bn - an) / std::max(bn, an);
    sil(n, 2) = bn;
  }
}

