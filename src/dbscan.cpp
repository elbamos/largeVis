// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

inline double dist2(const arma::mat& data, 
                    const double * x_i, 
                    long long id, 
                    int D) {
  double dist = 0;
  const double * x_j = data.colptr(id);
  for (int i = 0; i != D; i++) dist += (x_i[i] - x_j[i]) * (x_i[i] - x_j[i]);
  return dist;
}

typedef pair<int, double> iddist;

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

void frNNrecurse(const arma::mat& data, 
                 arma::imat& neighbors,
                 const double * x_i, 
                 int D, 
                 long long id,
                 double radius2,
                 NNheap& found,
                 std::set<long long>& checked) {
  double dist;
  for (arma::ivec::iterator it = neighbors.begin_col(id);
       it != neighbors.end_col(id) && *it != -1;
       it++) {
    if (checked.insert(*it).second) {
      const double * x_j = data.colptr(*it);
      dist = dist2(data, x_i, *it, D);
      if (dist < radius2) {
        found.push(iddist(*it, dist));
        frNNrecurse(data, neighbors, x_i, D, *it, radius2, found, checked);
      }
    }
  }
}

NNlist frNN(const arma::mat& data, 
            arma::imat& neighbors, 
            long long id,
            double radius2) {
  long long N = neighbors.n_cols;
  int D = data.n_rows;
  const double * x_i = data.colptr(id);
  
  NNheap found = NNheap();
  std::set<long long> checked = std::set<long long>();
  
  checked.insert(id);
  
  int k = 0;
  for (arma::ivec::iterator it = neighbors.begin_col(id); 
       it != neighbors.end_col(id) && *it != -1;
       it++) {
    frNNrecurse(data, neighbors, x_i, D, *it, radius2, 
                found, checked);
  }
  
  NNlist ret = NNlist();
  while (! found.empty()) {
    ret.push_back(found.top());
    found.pop();
  }
  return ret;
}

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

/* 
 *        DBSCAN
 */
// [[Rcpp::export]]
IntegerVector dbscan_int(const arma::mat& data, 
                         arma::imat& neighbors, 
                         double eps, 
                         int minPts) {

  // kd-tree uses squared distances
  double eps2 = eps*eps;
  long long N = data.n_cols;
  int D = data.n_rows;
  int K = neighbors.n_rows;
  
  // DBSCAN
  std::vector<bool> visited(N, false);
  std::vector< std::vector<int> > clusters; // vector of vectors == list
  NNlist  frNN1, frNN2;
  
  for (long long n = 0; n != N; n++) {
    if (!(n % 100)) Rcpp::checkUserInterrupt();
    
    if (visited[n]) continue;
    
    // start new cluster and expand
    frNN1 = frNN(data, neighbors, n, eps2);
    if (frNN1.size() < minPts) continue;
    std::vector<int> cluster;
    cluster.push_back(n);
    visited[n] = true;
    
    while (!frNN1.empty()) {
      int j = frNN1.back().first;
      frNN1.pop_back();
      
      if (visited[j]) continue; // point already processed
      visited[j] = true;
      
      frNN2 = frNN(data, neighbors, j, eps2);

      if (frNN2.size() >= minPts) { // expand neighborhood
        copy(frNN2.begin(), frNN2.end(),
                  back_inserter(frNN1));
        cluster.push_back(j);
      }
    }
    clusters.push_back(cluster);
  }
  
  // prepare cluster vector
  // unassigned points are noise (cluster 0)
  IntegerVector ret = IntegerVector(N,0);
  for (int c = 0; c < clusters.size(); c++) {
    std::vector<int> cluster = clusters[c];
    for (int j=0; j < cluster.size(); j++) {
      ret[cluster[j]] = c + 1 ;
    }
  }
  return ret;
}

/* 
 *          OPTICS
 */

void update(NNlist& frNeighbors,
            int p,
            std::vector<int> &seeds,
            int minPts,
            std::vector <bool> &visited,
            std::vector<int> &orderedPoints,
            std::vector<double> &reachdist,
            std::vector<double> &coredist) {
  
  std::vector<int>::iterator pos_seeds;
  double newreachdist;
  int o;
  double o_d;
  
  while(!frNeighbors.empty()) {
    o = frNeighbors.back().first; 
    o_d = frNeighbors.back().second;
    frNeighbors.pop_back(); 
    
    if(visited[o]) continue;
    
    newreachdist = std::max(coredist[p], o_d);
    
    if(reachdist[o] == INFINITY) {
      reachdist[o] = newreachdist;
      seeds.push_back(o);
    } else {
      if(newreachdist < reachdist[o]) reachdist[o] = newreachdist;
    }
  }
}


// [[Rcpp::export]]
List optics_int(const arma::mat& data, 
                arma::imat& neighbors,
                double eps, 
                int minPts) {
  
  // kd-tree uses squared distances
  double eps2 = eps*eps;
  long long N = data.n_cols;
  // OPTICS
  std::vector<bool> visited(N, false);
  std::vector<int> orderedPoints; orderedPoints.reserve(N);

  std::vector<double> reachdist(N, INFINITY); // we used Inf as undefined
  std::vector<double> coredist(N, INFINITY);
  
  std::vector<int> seeds;
  std::vector<double> ds;
  
  NNlist frNeighbors;
  
  for (int p = 0; p < N; p++) {
    if (!(p % 100)) Rcpp::checkUserInterrupt();
    Rprintf("processing point %d\n", p+1);
    
    if (visited[p]) continue;
    
    // ExpandClusterOrder
    //N = regionQueryDist(p, dataPts, kdTree, eps2, approx);
    frNeighbors = frNN(data, neighbors, p, eps2);
    
    visited[p] = true;
    
    // find core distance
    if(frNeighbors.size() >= (size_t) minPts) {
      coredist[p] = frNeighbors[minPts-1].second;
    }
    
    orderedPoints.push_back(p);
    
    if (coredist[p] == INFINITY) continue; // core-dist is undefined

    // update
    update(frNeighbors, p, seeds, minPts, visited, orderedPoints,
           reachdist, coredist);
    
    int q;
    while (!seeds.empty()) {
      // get smallest dist (to emulate priority queue). All should have already
      // a reachability distance <Inf from update().
      std::vector<int>::iterator q_it = seeds.begin();
      for (std::vector<int>::iterator it = seeds.begin();
           it!=seeds.end(); ++it) {
        // Note: The second part of the if statement ensures that ties are
        // always broken consistenty (higher ID wins to produce the same
        // results as the elki implementation)!
        if (reachdist[*it] < reachdist[*q_it] ||
            (reachdist[*it] == reachdist[*q_it] && *q_it < *it)) q_it = it;
      }
      q = *q_it;
      seeds.erase(q_it);
      
      //N2 = regionQueryDist(q, dataPts, kdTree, eps2, approx);
      frNeighbors = frNN(data, neighbors, q, eps2);
      visited[q] = true;
      
      // update core distance
      if(frNeighbors.size() >= (size_t) minPts) {
        coredist[q] = frNeighbors[minPts-1].second;
      }
      
      orderedPoints.push_back(q);
      
      if(frNeighbors.size() < (size_t) minPts) continue; //  == q has no core dist.
      
      // update seeds
      update(frNeighbors, q, seeds, minPts, visited, orderedPoints,
             reachdist, coredist);
    }
  }
  
  // prepare results (R index starts with 1)
  List ret;
  ret["order"] = IntegerVector(orderedPoints.begin(), orderedPoints.end())+1;
  ret["reachdist"] = sqrt(NumericVector(reachdist.begin(), reachdist.end()));
  ret["coredist"] = sqrt(NumericVector(coredist.begin(), coredist.end()));
  //ret["predecessor"] = IntegerVector(pre.begin(), pre.end())+1;
  return ret;
}