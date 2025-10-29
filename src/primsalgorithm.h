#include <RcppArmadillo.h>
#include "progress.hpp"
#include "minindexedpq.h"

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
		dist = fmax(d, dist);
		if (Q.decreaseIf(w, dist)) minimum_spanning_tree[w] = v;
		//  || w == starterIndex
	}

public:
	PrimsAlgorithm(const VIDX& N, const D* coreDistances) :
							N{N}, Q(PairingHeap<VIDX,D>(N)), coreDistances{coreDistances} {
  	minimum_spanning_tree = new VIDX[N];
  }

	PrimsAlgorithm(const PrimsAlgorithm& p) : PrimsAlgorithm(p.N, p.coreDistances) {};

	~PrimsAlgorithm() {
		delete[] minimum_spanning_tree;
	}

	VIDX* run(const arma::sp_mat& edges,
            const IntegerMatrix& neighbors,
            Progress& p,
            const VIDX& start) {
		starterIndex = start;
		for (VIDX n = 0; n != N; ++n) minimum_spanning_tree[n] = NA_INTEGER;
		Q.batchInsert(N, start);
		Q.decreaseIf(starterIndex, -1);
		while (! Q.isEmpty()) {
			VIDX v = Q.pop();
			if (! p.increment()) break;
			IntegerVector vNeighbors = neighbors.column(v);
			for (auto it = vNeighbors.begin();
        it != vNeighbors.end() && *it != -1;
        ++it) {
				updateVWD(v, *it, edges(v, *it));
			}
			for (arma::sp_mat::const_col_iterator it = edges.begin_col(v);
        it != edges.end_col(v);
        ++it) {
				updateVWD(v, it.row(), *it);
        }
		}
		return minimum_spanning_tree;
	}

	std::vector< std::pair<D, VIDX> > getMergeSequence() const {
		std::vector< std::pair<D, VIDX> > container;
		container.reserve(N);
		for (VIDX n = 0; n != N; ++n) container.emplace_back(Q.keyOf(n), n);
		sort(container.begin(), container.end());
		return container;
	}
};