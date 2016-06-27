// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.hpp"

using namespace Rcpp;
using namespace std;
using namespace arma;

/*
 * The recursive function for the annoy neighbor search
 * algorithm. Partitions space by a random hyperplane,
 * and calls itself recursively (twice) on each side.
 *
 * If called with fewer nodes than the threshold,
 *
 */
void searchTree(const int& threshold,
                const arma::ivec& indices,
                const arma::mat& data,
                neighborhood* heap[],
                                  const int& iterations,
                                  Progress& progress) {
  const int I = indices.size();
  // const int D = data.n_rows;
  if (progress.check_abort()) return;
  if (I < 2) stop("Tree split failure.");
  if (I <= threshold || iterations == 0) {
    ivec neighbors = ivec(indices);
    neighborhood tmpStorage = neighborhood();
    ivec::iterator newEnd = neighbors.end();
#ifdef _OPENMP
#pragma omp critical
#endif
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
#ifdef _OPENMP
#pragma omp parallel for shared(treeNeighborhoods)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for shared(treeHolder, treeNeighborhoods) private(thisHeap)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
#ifdef _OPENMP
#pragma omp parallel for shared(old_knns, knns) private(thisHeap, sorter)
#endif
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

// Sparse version
void searchTree(const int& threshold,
                const arma::vec& indices,
                const sp_mat& data,
                vector<vector<int>* >& heap,
                const int& iterations,
                Progress& progress) {
  const int I = indices.size();
  // const int D = data.n_rows;
  if (progress.check_abort()) return;
  if (I < 2) stop("Tree split failure.");
  if (I == 2) {
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      heap[indices[0]] -> push_back(indices[1]);
      heap[indices[1]] -> push_back(indices[0]);
    }
        return;
  }
  if (I < threshold || iterations == 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
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
#ifdef _OPENMP
#pragma omp parallel for shared(indices,treeNeighborhoods)
#endif
    for (int t = 0; t < n_trees; t++) if (! p.check_abort()) {
      searchTree(threshold,
                 indices,
                 data,
                 treeNeighborhoods,
                 max_recursion_degree, // maximum permitted level of recursion
                 p
      );

      if (t > 0 && ! p.check_abort())
#ifdef _OPENMP
#pragma omp critical
#endif
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
#ifdef _OPENMP
#pragma omp parallel for shared(knns)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for shared(knns, treeNeighborhoods)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for shared(pjis, rowSums)
#endif
  for (int e=0; e < pjis.size(); e++) if (p.increment()){
    const int i = is[e];
    const double pji = exp(- pow(xs[e], 2)) / sigmas[i];
    pjis[e] = pji;
#ifdef _OPENMP
#pragma omp atomic
#endif
    rowSums[i] += pji;
  }
  if (p.check_abort()) return sp_mat(0);
  // Now convert p(j|i) to w_{ij} by symmetrizing.
  // Loop through the edges, and populate a location matrix and value vector for
  // the sp_mat batch insertion constructor.  Put all coordinates in the
  // lower triangle.  The constructor will automatically add duplicates.
  vec values = vec(pjis.size());
  umat locations = umat(2, pjis.size());
#ifdef _OPENMP
#pragma omp parallel for shared(locations, values)
#endif
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
  // dxs_ds = xs * pow(x_i, 2)
  const NumericVector softxs = xs / sum(xs);
  // dsoftxs_ds =
  const double p2 = - sum(log(softxs) / log(2)) / xs.length();
  return pow(perplexity - p2, 2);
};

/*
 * The stochastic gradient descent function.
 */
// [[Rcpp::export]]
arma::mat sgd(arma::mat coords,
              arma::ivec& is, // vary randomly
              const IntegerVector js, // ordered
              const IntegerVector ps, // N+1 length vector of indices to start of each row j in vector is
              const arma::vec ws, // w{ij}
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
  if (D > 10) stop("Low dimensional space cannot have more than 10 dimensions.");
  const int N = ps.size() - 1;
  const int E = ws.size();
  double *coordsPtr = coords.memptr();

  double* negProb = new double[N];
  int* negAlias = new int[N];
  makeAliasTable(N, pow(diff(ps), 0.75), negProb, negAlias);
  double* posProb = new double[E];
  int* posAlias = new int[E];
  if (! useWeights) makeAliasTable(E, ws, posProb, posAlias);
  else for (int i = 0; i < E; i++) posAlias[i] = 1;

  const int posSampleLength = ((nBatches > 1000000) ? 1000000 : (int) nBatches);
  mat positiveSamples = randu<mat>(2, posSampleLength);
  double *posRandomPtr = positiveSamples.memptr();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) shared (coords, positiveSamples)
#endif
  for (long eIdx=0; eIdx < nBatches; eIdx++) if (progress.increment()) {

    const int e_ij = searchAliasTable(E, posRandomPtr,
                                      posProb, posAlias,
                                      eIdx % posSampleLength);

    const int i = is[e_ij];
    const int j = js[e_ij];

    // mix weight into learning rate
    const double localRho =  ((useWeights) ? ws[e_ij] : 1.0) * (rho - ((rho - minRho) * eIdx / nBatches));

    double *y_i = coordsPtr + (i * D);
    double *y_j = coordsPtr + (j * D);

    double firstholder[10];
    double secondholder[10];

    positiveGradient(y_i, y_j, firstholder, alpha, D);

    for (int d = 0; d < D; d++) y_j[d] -= firstholder[d] * localRho;

    mat negSamples = mat(2, M * 2);
    double *samplesPtr = negSamples.memptr();
    int sampleIdx = 0;
    ivec searchRange = is.subvec(ps[i], ps[i + 1] - 1);
    ivec::iterator searchBegin = searchRange.begin();
    ivec::iterator searchEnd = searchRange.end();
    int m = 0;
    int k;
    while (m < M) {
      if (sampleIdx % (M * 2) == 0) negSamples.randu();
      k = searchAliasTable(N, samplesPtr, negProb, negAlias, sampleIdx++ % (M * 2));
      // Check that the draw isn't one of i's edges
      if (k == i ||
          k == j ||
          binary_search( searchBegin,
                         searchEnd,
                         k)) continue;

      double *y_k = coordsPtr + (k * D);

      if (negativeGradient(y_i, y_k, secondholder, alpha, gamma, D)) continue;

      for (int d = 0; d < D; d++) firstholder[d] += secondholder[d];
      for (int d = 0; d < D; d++) y_k[d] -= secondholder[d] * localRho;

      m++;
      if (sampleIdx > M * 10) stop("Bad sampleidx");
    }
    for (int d = 0; d < D; d++) y_i[d] += firstholder[d] * localRho;

    if (eIdx > 0 &&
        eIdx % posSampleLength == 0) positiveSamples.randu();
  }
  return coords;
};
