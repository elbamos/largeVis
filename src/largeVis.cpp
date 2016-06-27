#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
#include <omp.h>
#include "progress.hpp"
#include <math.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <set>
using namespace Rcpp;
using namespace std;
using namespace arma;


struct heapObject {
  double d;
  int n;

  heapObject(double d, int n) : d(d), n(n) {}

  bool operator<(const struct heapObject& other) const {
    return d < other.d;
  }
};

typedef priority_queue<heapObject> maxHeap;
typedef vector< imat::col_iterator > positionVector;
typedef vector<int> neighborhood;

// The Euclidean distance between two vectors

double relDist(const arma::vec& i, const arma::vec& j) {
  const int lim = i.n_elem;
  double cnt = 0;
  for (int idx = 0; idx < lim; idx++) cnt += ((i[idx] - j[idx]) * (i[idx] - j[idx]));
  return cnt;
}

double dist(const arma::vec& i, const arma::vec& j) {
  return sqrt(relDist(i,j));
}

// Stolen directly from Erik B's annoylib.h
double cosDist(const arma::vec& i, const arma::vec& j) {
  int lim = i.n_elem;
  double pp = 0, qq = 0, pq = 0;
  for (int z = 0; z < lim; z++) {
    pp += (i[z]) * (i[z]);
    qq += (j[z]) * (j[z]);
    pq += (i[z]) * (j[z]);
  }
  double ppqq = pp * qq;
  if (ppqq > 0) return 2.0 - 2.0 * pq / sqrt(ppqq);
  else return 2.0; // cos is 0
}

double sparseDist(const sp_mat& i, const sp_mat& j) {
  return as_scalar(sqrt(sum(square(i - j))));
}

double sparseCosDist(const sp_mat& i, const sp_mat& j) {
  return 1 - (as_scalar((dot(i,j)) / as_scalar(norm(i,2) * norm(j,2))));
}

double sparseRelDist(const sp_mat& i, const sp_mat& j) {
  return as_scalar(sum(square(i - j)));
}

void searchTree(const int& threshold,
                const arma::ivec& indices,
                const arma::mat& data,
                neighborhood* heap[],
                                  const int& iterations,
                                  Progress& progress) {
  const int I = indices.size();
  const int D = data.n_rows;
  if (progress.check_abort()) return;
  if (I < 2) stop("Tree split failure.");
  if (I <= threshold || iterations == 0) {
    ivec neighbors = ivec(indices);
    neighborhood tmpStorage = neighborhood();
    ivec::iterator newEnd = neighbors.end();
#pragma omp critical
{
  for (ivec::iterator it = neighbors.begin();
       it != newEnd;
       it++) {
    tmpStorage.clear();
    tmpStorage.swap(*heap[*it]);
    heap[*it] -> reserve(tmpStorage.size() + I);
    ivec::iterator newIt = neighbors.begin();
    vector<int>::iterator oldIt = tmpStorage.begin();
    vector<int>::iterator oldEnd = tmpStorage.end();
    int last;
    int best = -1;
    while (oldIt != oldEnd || newIt != newEnd) {
      if (oldIt == oldEnd) best = *newIt++;
      else if (newIt == newEnd) best = *oldIt++;
      else best = (*newIt < *oldIt) ? *newIt++ : *oldIt++;
      if (best == last || best == *it) continue;
      heap[*it] -> push_back(best);
      last = best;
    }
  }
}
progress.increment(I);
return;
  }
  vec direction = vec(indices.size());
  {
    int x1idx, x2idx;
    vec v;
    vec m;
    do {
      const vec selections = randu(2) * (I - 1);
      x1idx = indices[selections[0]];
      x2idx = indices[selections[1]];
      if (x1idx == x2idx) x2idx = indices[((int)selections[1] + 1) % indices.size()];
      const vec x2 = data.col(x2idx);
      const vec x1 = data.col(x1idx);
      // Get hyperplane
      m =  (x1 + x2) / 2; // Base point of hyperplane
      const vec d = x1 - x2;
      v =  d / as_scalar(norm(d, 2)); // unit vector
    } while (x1idx == x2idx);

    for (int i = 0; i < indices.size(); i++) {
      const vec X = data.col(indices[i]);
      direction[i] = dot((X - m), v);
    }
  }
  // Normalize direction
  const double middle = median(direction);

  const uvec left = find(direction > middle);
  const uvec right = find(direction <= middle);
  if (left.size() >= 2 && right.size() >= 2) {
    searchTree(threshold, indices(left), data, heap, iterations - 1, progress);
    searchTree(threshold, indices(right), data, heap, iterations - 1, progress);
  } else { // Handles the rare case where the split fails because of equidistant points
    searchTree(threshold, indices.subvec(0, indices.size() / 2), data, heap, iterations - 1, progress);
    searchTree(threshold, indices.subvec(indices.size() / 2, indices.size() - 1), data, heap, iterations - 1, progress);
  }
};

// [[Rcpp::export]]
arma::imat searchTrees(const int& threshold,
                       const int& n_trees,
                       const int& K,
                       const int& max_recursion_degree,
                       const int& maxIter,
                       const arma::mat& data,
                       const std::string& distMethod,
                       bool verbose) {

  const int N = data.n_cols;

  double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = relDist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = cosDist;
  else distanceFunction = relDist;

  Progress p((N * n_trees) + (2 * N) + (N * maxIter), verbose);

  imat knns;
  {
    set<int>** treeHolder = new set<int>*[N];
    double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j);
    if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = relDist;
    else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = cosDist;
    else distanceFunction = relDist;

    neighborhood** treeNeighborhoods = new neighborhood*[N];
    for (int i = 0; i < N; i++) {
      int seed[] = {i};
      treeNeighborhoods[i] = new vector<int>(seed, seed + sizeof(seed) / sizeof(int));
    }
    const ivec indices = regspace<ivec>(0, N - 1);

#pragma omp parallel for shared(treeNeighborhoods)
    for (int t = 0; t < n_trees; t++) if (! p.check_abort()) {
      searchTree(threshold,
                 indices,
                 data,
                 treeNeighborhoods,
                 max_recursion_degree, // maximum permitted level of recursion
                 p
      );
    }
    if (p.check_abort()) return knns;
    // Reduce size from threshold * n_trees to top K, and sort
    maxHeap thisHeap = maxHeap();
#pragma omp parallel for shared(treeHolder, treeNeighborhoods) private(thisHeap)
    for (int i = 0; i < N; i++) if (p.increment()) {
      const vec x_i = data.col(i);
      vector<int> *neighborhood = treeNeighborhoods[i];
      for (vector<int>::iterator j = neighborhood -> begin();
           j != neighborhood -> end();
           j++) {

        const double d = distanceFunction(x_i, data.col(*j));
        if (d != 0) {
          thisHeap.emplace(d, *j);
          if (thisHeap.size() > K) thisHeap.pop();
        }
      }
      delete treeNeighborhoods[i];
      treeHolder[i] = new set<int>();
      while (! thisHeap.empty()) {
        treeHolder[i] -> emplace(thisHeap.top().n);
        thisHeap.pop();
      }
    }
    // Copy sorted neighborhoods into matrix. This is faster than
    // sorting in-place.
    knns = imat(K,N);
#pragma omp parallel for
    for (int i = 0; i < N; i++) if (p.increment()) {
      set<int>::iterator sortIterator = treeHolder[i] -> begin();
      set<int>::iterator end = treeHolder[i] -> end();
      int j = 0;
      while (sortIterator != end) knns(j++, i) = *sortIterator++;
      if (j == 0) stop("Tree failure.");
      while (j < K) knns(j++,i) = -1;
      delete treeHolder[i];
    }
  }

  if (p.check_abort()) return imat(0);
  imat old_knns  = imat(K,N);
  for (int T = 0; T < maxIter; T++) if (! p.check_abort()) {
    imat tmp = old_knns;
    old_knns = knns;
    knns = tmp;
    maxHeap thisHeap = maxHeap();
    set<int> sorter = set<int>();
#pragma omp parallel for shared(old_knns, knns) private(thisHeap, sorter)
    for (int i = 0; i < N; i++) if (p.increment()) {
      double d;

      const vec x_i = data.col(i);

      positionVector positions = positionVector();
      positionVector ends = positionVector();

      positions.reserve(K + 1);
      ends.reserve(K + 1);

      positions.push_back(old_knns.begin_col(i));
      ends.push_back(old_knns.end_col(i));

      for (imat::col_iterator it = old_knns.begin_col(i);
           it != ends[0] && *it != -1;
           it++) {
        positions.push_back(old_knns.begin_col(*it));
        ends.push_back(old_knns.end_col(*it));
      }

      int lastOne = N + 1;
      // This is a K + 1 vector merge sort running in O(K * N)
      positionVector::iterator theEnd = positions.end();
      while (true) {
        imat::col_iterator whch = 0;

        for (pair< positionVector::iterator,
             positionVector::iterator > it(positions.begin(),
                                           ends.begin());
             it.first != theEnd;
             it.first++, it.second++) while (*it.first != *it.second) { // For each neighborhood, keep going until
               // we find a non-dupe or get to the end

               if (**it.first == -1) advance(*it.first, distance(*it.first, *it.second));
               else if (**it.first == i || **it.first == lastOne) advance(*it.first, 1);
               else if (whch == 0 || **it.first < *whch) {
                 whch = *it.first;
                 break;
               } else break;
             }

             if (whch == 0) break;
             lastOne = *whch;
             advance(whch, 1);

             d = distanceFunction(x_i, data.col(lastOne));
             if (d > 0) {
               thisHeap.emplace(d, lastOne);
               if (thisHeap.size() > K) thisHeap.pop();
             }
      }

      sorter.clear();
      while (!thisHeap.empty()) {
        sorter.emplace(thisHeap.top().n);
        thisHeap.pop();
      }
      set<int>::iterator sortIterator = sorter.begin();
      int j = 0;
      while (sortIterator != sorter.end()) knns(j++, i) = *sortIterator++;
      if (j == 0) stop("Neighbor exploration failure.");
      while (j < K) knns(j++,i) = -1;
    }
  }
  return knns;
};

void searchTree(const int& threshold,
                const arma::vec& indices,
                const sp_mat& data,
                vector<vector<int>* >& heap,
                const int& iterations,
                Progress& progress) {
  const int I = indices.size();
  const int D = data.n_rows;
  if (progress.check_abort()) return;
  if (I < 2) stop("Tree split failure.");
  if (I == 2) {
#pragma omp critical
{
  heap[indices[0]] -> push_back(indices[1]);
  heap[indices[1]] -> push_back(indices[0]);
}
    return;
  }
  if (I < threshold || iterations == 0) {
#pragma omp critical
{
  for (int i = 0; i < I; i++) {
    heap[indices[i]] -> reserve(I - 1);
    for (int j = 0; j < I; j++) if (i != j) heap[indices[i]] -> push_back(indices[j]);
  }
}
    progress.increment(I);
    return;
  }
  vec direction = vec(indices.size());
  {
    int x1idx, x2idx;
    sp_mat v;
    sp_mat m;
    do {
      const vec selections = randu(2) * (I - 1);
      x1idx = indices[selections[0]];
      x2idx = indices[selections[1]];
      if (x1idx == x2idx) x2idx = indices[((int)selections[1] + 1) % indices.size()];
      const SpSubview<double> x2 = data.col(x2idx);
      const SpSubview<double> x1 = data.col(x1idx);
      // Get hyperplane
      m =  (x1 + x2) / 2; // Base point of hyperplane
      const sp_mat d = x1 - x2;
      const double dn = as_scalar(norm(d, 2));
      v =  d / dn; // unit vector
    } while (x1idx == x2idx);

    for (int i = 0; i < indices.size(); i++) {
      const int I = indices[i];
      const SpSubview<double> X = data.col(I);
      direction[i] = dot((X - m), v);
    }
  }
  // Normalize direction
  const double middle = median(direction);

  const uvec left = find(direction > middle);
  const uvec right = find(direction <= middle);
  if (left.size() >= 2 && right.size() >= 2) {
    searchTree(threshold, indices(left), data, heap, iterations - 1, progress);
    searchTree(threshold, indices(right), data, heap, iterations - 1, progress);
  } else {
    searchTree(threshold, indices.subvec(0, indices.size() / 2), data, heap, iterations - 1, progress);
    searchTree(threshold, indices.subvec(indices.size() / 2, indices.size() - 1), data, heap, iterations - 1, progress);
  }
};

arma::mat searchTreesSparse(const int& threshold,
                            const int& n_trees,
                            const int& K,
                            const int& max_recursion_degree,
                            const int& maxIter,
                            const sp_mat& data,
                            const std::string& distMethod,
                            bool verbose) {

  const int N = data.n_cols;

  double (*distanceFunction)(const sp_mat& x_i, const sp_mat& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = sparseRelDist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = sparseCosDist;
  else distanceFunction = sparseRelDist;

  Progress p((N * n_trees) + (N) + (N * maxIter), verbose);

  vector<vector<int>* > treeNeighborhoods = vector<vector<int>* >(N);
  for (int i = 0; i < N; i++) {
    int seed[] = {i};
    treeNeighborhoods[i] = new vector<int>(seed, seed + sizeof(seed) / sizeof(int));
  }

  { // Artificial scope to destroy indices
    vec indices = regspace<vec>(0, N - 1);

#pragma omp parallel for shared(indices,treeNeighborhoods)
    for (int t = 0; t < n_trees; t++) if (! p.check_abort()) {
      searchTree(threshold,
                 indices,
                 data,
                 treeNeighborhoods,
                 max_recursion_degree, // maximum permitted level of recursion
                 p
      );

      if (t > 0 && ! p.check_abort())
#pragma omp critical
{
  for (int i = 0; i < N; i++) {
    vector<int>* neighbors = treeNeighborhoods[i];
    sort(neighbors -> begin(), neighbors -> end());
    vector<int>::iterator theEnd = unique(neighbors -> begin(), neighbors -> end());
    neighbors -> erase(theEnd, neighbors -> end());
    if (neighbors -> size() < 3) stop("Tree failure.");
  }
}
    }
  }

  if (p.check_abort()) return mat(0);

  // Initialize the knn matrix, and reduce the number of candidate neighbors per node
  // to K.  Otherwise the first neighborhood exploration pass takes N * trees * (threshold + 1),
  // instead of (N * K), which is prohibitive of large thresholds.
  mat knns = mat(threshold,N);
  knns.fill(-1);
#pragma omp parallel for shared(knns)
  for (int i = 0; i < N; i++) if (p.increment()){
    const SpSubview<double> x_i = data.col(i);
    priority_queue<heapObject> maxHeap = priority_queue<heapObject>();
    vector<int>* stack = treeNeighborhoods[i];
    for (vector<int>::iterator it = stack -> begin(); it != stack -> end(); it++) {
      const double d = distanceFunction(x_i, data.col(*it));
      maxHeap.push(heapObject(d, *it));
      if (maxHeap.size() > threshold) maxHeap.pop();
    }
    int j = 0;
    do {
      knns(j,i) = maxHeap.top().n;
      maxHeap.pop();
      j++;
    } while (j < threshold && ! maxHeap.empty());
    if (j == 1 && knns(0,i) == -1) stop("Bad neighbor matrix.");
  }
  if (p.check_abort()) return mat(0);

  for (int T = 0; T < maxIter; T++) {
    mat old_knns = knns;
    knns = mat(K,N);
    knns.fill(-1);
#pragma omp parallel for shared(knns, treeNeighborhoods)
    for (int i = 0; i < N; i++) if (p.increment()) {
      double d;

      const vec neighborhood = old_knns.col(i);
      const SpSubview<double> x_i = data.col(i);

      priority_queue<heapObject> heap;
      vector<int> pastVisitors = *(treeNeighborhoods[i]);
      pastVisitors.reserve((K + 1) * K);
      // Loop through immediate neighbors of i
      for (int jidx = 0; jidx < old_knns.n_rows; jidx++) {
        const int j = neighborhood[jidx];
        if (j == -1) break;
        if (j == i) continue; // This should never happen
        d = distanceFunction(x_i, data.col(j));
        if (d == 0) continue; // duplicate
        heap.push(heapObject(d, j));
        if (heap.size() > K) heap.pop();

        // For each immediate neighbor j, loop through its neighbors
        const vec locality = old_knns.col(j);
        for (int kidx = 0; kidx < old_knns.n_rows; kidx++) {
          const int k = locality[kidx];
          if (k == -1) break;
          if (k == i) continue;
          // Check if this is a neighbor we've already seen.  O(log k)
          pair<vector<int>::iterator,
                    vector<int>::iterator > firstlast = equal_range(pastVisitors.begin(),
                                                                              pastVisitors.end(),
                                                                              k);
          if (*(firstlast.first) == k) continue; // Found

          if (firstlast.second == pastVisitors.end()) pastVisitors.push_back(k);
          else pastVisitors.insert(firstlast.second, k);

          d = distanceFunction(x_i, data.col(k));
          if (d == 0) continue;
          if (heap.size() < K) heap.push(heapObject(d,k));
          else if (d < heap.top().d) {
            heap.push(heapObject(d, k));
            if (heap.size() > K) heap.pop();
          }
        }
      }
      int j = 0;
      while (j < K && ! heap.empty()) {
        knns(j, i) = heap.top().n;
        heap.pop();
        j++;
      }
      if (j == 0) stop("Failure in neighborhood exploration - this should never happen.");
      vector<int>(pastVisitors).swap(pastVisitors); // pre-C++11 shrink
    }
  }
  return knns;
};

// [[Rcpp::export]]
arma::mat searchTreesCSparse(const int& threshold,
                             const int& n_trees,
                             const int& K,
                             const int& max_recursion_degree,
                             const int& maxIter,
                             const arma::uvec& i,
                             const arma::uvec& p,
                             const arma::vec& x,
                             const std::string& distMethod,
                             bool verbose) {
  const int N = p.size() -1;
  const sp_mat data = sp_mat(i,p,x,N,N);
  return searchTreesSparse(threshold,n_trees,K,max_recursion_degree,maxIter,data,distMethod,verbose);
}

// [[Rcpp::export]]
arma::mat searchTreesTSparse(const int& threshold,
                             const int& n_trees,
                             const int& K,
                             const int& max_recursion_degree,
                             const int& maxIter,
                             const arma::uvec& i,
                             const arma::uvec& j,
                             const arma::vec& x,
                             const std::string& distMethod,
                             bool verbose) {
  const umat locations = join_cols(i,j);
  const sp_mat data = sp_mat(locations,x);
  return searchTreesSparse(threshold,n_trees,K,max_recursion_degree,maxIter,data,distMethod,verbose);
}


/*
 * Fast calculation of pairwise distances with the result stored in a pre-allocated vector.
 */
// [[Rcpp::export]]
arma::vec fastDistance(const NumericVector is,
                       const NumericVector js,
                       const arma::mat& data,
                       const std::string& distMethod,
                       bool verbose) {

  Progress p(is.size(), verbose);
  vec xs = vec(is.size());
  double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = dist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = cosDist;

#pragma omp parallel for shared (xs)
  for (int i=0; i < is.length(); i++) if (p.increment()) xs[i] =
    distanceFunction(data.col(is[i]), data.col(js[i]));
  return xs;
};

arma::vec fastSparseDistance(const arma::vec& is,
                             const arma::vec& js,
                             const sp_mat& data,
                             const std::string& distMethod,
                             bool verbose) {

  Progress p(is.size(), verbose);
  vec xs = vec(is.size());
  double (*distanceFunction)(
      const sp_mat& x_i,
      const sp_mat& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = sparseDist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = sparseCosDist;

#pragma omp parallel for shared (xs)
  for (int i=0; i < is.size(); i++) if (p.increment()) xs[i] =
    distanceFunction(data.col(is[i]), data.col(js[i]));
  return xs;
};

// [[Rcpp::export]]
arma::vec fastCDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& p_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose) {
  const int N = p_locations.size() - 1;
  const sp_mat data = sp_mat(i_locations, p_locations, x, N, N);
  return fastSparseDistance(is,js,data,distMethod,verbose);
}

// [[Rcpp::export]]
arma::vec fastSDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& j_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose) {
  const umat locations = join_cols(i_locations, j_locations);
  const sp_mat data = sp_mat(locations, x);
  return fastSparseDistance(is,js,data,distMethod,verbose);
}

// Take four vectors (i indices, j indices, edge distances, and sigmas), and calculate
// p(j|i) and then w_{ij}.
// [[Rcpp::export]]
arma::sp_mat distMatrixTowij(
    const NumericVector is,
    const NumericVector js,
    const NumericVector xs,
    const NumericVector sigmas,
    const int N,
    bool verbose
) {

  Progress p(xs.size() * 2, verbose);
  vec rowSums = vec(N);
  vec pjis = vec(is.length());
  for (int idx=0; idx < N; idx++) rowSums[idx] = 0;
  // Compute pji, accumulate rowSums at the same time
#pragma omp parallel for shared(pjis, rowSums)
  for (int e=0; e < pjis.size(); e++) if (p.increment()){
    const int i = is[e];
    const double pji = exp(- pow(xs[e], 2)) / sigmas[i];
    pjis[e] = pji;
#pragma omp atomic
    rowSums[i] += pji;
  }
  if (p.check_abort()) return sp_mat(0);
  // Now convert p(j|i) to w_{ij} by symmetrizing.
  // Loop through the edges, and populate a location matrix and value vector for
  // the sp_mat batch insertion constructor.  Put all coordinates in the
  // lower triangle.  The constructor will automatically add duplicates.
  vec values = vec(pjis.size());
  umat locations = umat(2, pjis.size());
#pragma omp parallel for shared(locations, values)
  for (int e = 0; e < pjis.size(); e++) if (p.increment()) {
    int newi = is[e], newj = js[e];
    if (newi < newj) swap(newi, newj);
    values[e] =  ((pjis[e] / rowSums[is[e]]) / (2 * N));
    locations(1,e) = newi;
    locations(0,e) = newj;
  }
  sp_mat wij = sp_mat(
    true, // add_values
    locations,
    values,
    N, N // n_col and n_row
  );
  wij = wij + wij.t();
  return wij;
};


// [[Rcpp::export]]
double sigFunc(const double& sigma,
               const NumericVector& x_i,
               const double& perplexity) {
  const NumericVector xs = exp(- pow(x_i,2) / sigma);
  const NumericVector softxs = xs / sum(xs);
  const double p2 = - sum(log(softxs) / log(2)) / xs.length();
  return pow(perplexity - p2, 2);
};

/*
 * Some helper functions useful in debugging.
 */
void checkVector(const arma::vec& x,
                 const std::string& label) {
  if (x.has_nan() || x.has_inf())
    Rcout << "\n Failure at " << label;
};

double objective(const arma::mat& inputs, double gamma, double alpha) {
  double objective = log(1 / (1 + (alpha * dist(inputs.col(0), inputs.col(1)))));
  for (int i = 2; i < 7; i++)
    objective += gamma * log(1 - (1 / (1 + alpha * dist(inputs.col(0), inputs.col(i)))));
  return objective;
}

void checkGrad(const arma::vec& x,
               const arma::vec& y,
               const arma::vec& grad,
               bool together,
               const std::string& label) {
  double oldDist = dist(x,y);
  double newDist = dist(x + grad, y - grad);
  if (together && newDist > oldDist) Rcout << "\nGrad " << label << " yi " << x << " other " << y << " grad " << grad << " moved further apart.";
  else if (! together && newDist < oldDist) Rcout << "\nGrad " << label << " yi " << x << " other " << y << " grad " << grad << "moved closer together.";
};

/*
 * The stochastic gradient descent function. Asynchronicity is enabled by openmp.
 */

// [[Rcpp::export]]
arma::mat sgd(arma::mat coords,
              arma::ivec& is, // vary randomly
              const IntegerVector js, // ordered
              const IntegerVector ps, // N+1 length vector of indices to start of each row j in vector is
              const NumericVector ws, // w{ij}
              const double gamma,
              const double rho,
              const double minRho,
              const bool useWeights,
              const long nBatches,
              const int M,
              const double alpha,
              bool verbose) {

  Progress progress(nBatches, verbose);

  const int D = coords.n_rows;
  const int N = ps.size() - 1;
  const int E = ws.length();
  // Calculate negative sample weights, d_{i}^0.75.
  // Stored as a vector of cumulative sums, normalized, so it can
  // be readily searched using binary searches.
  vec negativeSampleWeights = pow(diff(ps), 0.75);
  const double scale = sum(negativeSampleWeights);
  negativeSampleWeights = negativeSampleWeights / scale;
  negativeSampleWeights = cumsum(negativeSampleWeights);

  // positive edges for sampling
  vec positiveEdgeWeights;
  if (! useWeights) {
    const double posScale = sum(ws);
    positiveEdgeWeights = vec(E);
    positiveEdgeWeights[0] = ws[0] / posScale;
    for (int idx = 1; idx < E; idx++)
      positiveEdgeWeights[idx] = positiveEdgeWeights[idx - 1] + (ws[idx] / posScale);
  }

  const int posSampleLength = ((nBatches > 1000000) ? 1000000 : (int) nBatches);
  vec positiveSamples = randu<vec>(posSampleLength);

  // Cache some variables to make memory allocation cheaper
  vec::iterator posEnd = positiveEdgeWeights.end();
  vec::iterator negBegin = negativeSampleWeights.begin();
  vec::iterator negEnd = negativeSampleWeights.end();
  ivec::iterator searchBegin = is.begin();
  // Iterate through the edges in the positiveEdges vector
#pragma omp parallel for shared(coords, positiveSamples) schedule(static)
  for (long eIdx=0; eIdx < nBatches; eIdx++) {
    if (progress.increment()) {
      const double posTarget = *(positiveSamples.begin() + (eIdx % posSampleLength));
      int k;
      int e_ij;
      if (useWeights) {
        e_ij = posTarget * (E - 1);
      } else {
        e_ij = distance(positiveEdgeWeights.begin(),
                             upper_bound(positiveEdgeWeights.begin(),
                                              posEnd,
                                              posTarget));
      }
      const int i = is[e_ij];
      const int j = js[e_ij];

      const double localRho = rho - ((rho - minRho) * eIdx / nBatches);

      //if ((randn<vec>(1))[0] < 0) swap(i, j);

      const vec y_i = coords.col(i);
      const vec y_j = coords.col(j);

      // wij
      const double w = (useWeights) ? ws[e_ij] : 1;

      const double dist_ij = dist(y_i, y_j);

      const vec d_dist_ij = (y_i - y_j) / sqrt(dist_ij);
      double p_ij;
      if (alpha == 0)   p_ij =   1 / (1 +      exp(dist_ij));
      else              p_ij =   1 / (1 + (alpha * dist_ij));

      vec d_p_ij;
      if (alpha == 0) d_p_ij =  d_dist_ij * -2 * dist_ij * exp(dist_ij) / pow(1 + exp(dist_ij), 2);
      else            d_p_ij =  d_dist_ij * -2 * dist_ij * alpha        / pow(1 +    (dist_ij * alpha),2);

      //double o = log(p_ij);
      const vec d_j = (1 / p_ij) * d_p_ij;
      // alternative: d_i - 2 * alpha * (y_i - y_j) / (alpha * sum(square(y_i - y_j)))
      vec d_i = d_j;

      // Setup negative search
      vec samples = sort(randu<vec>(M));
      vec::iterator targetIt = samples.begin();
      int sampleIdx = 0;
      // The indices of the nodes with edges to i
      int m = 0;
      vec::iterator negativeIterator = negBegin;
      while (m < M) {
        if (sampleIdx % M == 0) {
          samples.randu();
          samples = sort(samples);
          targetIt = samples.begin();
          negativeIterator = negBegin;
        }
        // binary search implementing weighted sampling
        const double target = targetIt[sampleIdx++ % M];
        int k;
        if (useWeights) k = target * (N - 1);
        else {
          negativeIterator = upper_bound(negativeIterator,
                                              negEnd,
                                              target);
          k = negativeIterator - negBegin;
        }

        if (k == i ||
            k == j ||
            binary_search(searchBegin + ps[i],
                               searchBegin + ps[i + 1] - 1,
                               k)) continue;
        const vec y_k = coords.col(k);

        const double dist_ik = dist(y_i, y_k);
        if (dist_ik == 0 || dist_ik > (14 * (1 + (int) (sampleIdx / M)))) continue; // Duplicates and distances too large to care

        const vec d_dist_ik = (y_i - y_k) / sqrt(dist_ik);

        double p_ik;
        if (alpha == 0) p_ik  =  1 - (1 / (1 +      exp(dist_ik)));
        else            p_ik  =  1 - (1 / (1 + (alpha * dist_ik)));

        vec d_p_ik;
        if (alpha == 0) d_p_ik =  d_dist_ik * 2 * dist_ik * exp(dist_ik) / pow(1 +      exp(dist_ik),2);
        else            d_p_ik =  d_dist_ik * 2 * dist_ik * alpha        / pow(1 + (alpha * dist_ik),2);
        //o += (gamma * log(p_ik));

        const vec d_k = (gamma / p_ik) * d_p_ik;
        // alternative:  d_k = 2 * alpha * (y_i - y_k) / (square(1 + (alpha * sum(square(y_i - y_k)))) * (1 - (1 / (alpha * sum(square(y_i - y_k))))))

        d_i += d_k;
        for (int idx = 0; idx < D; idx++) coords(idx,k) -= d_k[idx] * localRho * w;

        m++;
      }

      for (int idx = 0; idx < D; idx++) {
        coords(idx,j) -=  d_j[idx] * w * localRho;
        coords(idx,i) +=  d_i[idx] * w * localRho;
      }

      if (eIdx > 0 &&
          eIdx % posSampleLength == 0) positiveSamples.randu();
    }
  }
  return coords;
};
