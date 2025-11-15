#include "largeVis.h"
#include <random>
#include <queue>

using namespace std;

template <class T, // number of aliases
          class C, // coordinate type, used for probabilities
          class D> // distance type, used to take weights on initialization
class AliasTable final {
private:
	C* probs;
	T* aliases;
	mt19937_64 mt;
	T N;

	// Platform-independent uniform random number in [0,1)
	inline C uniformRandom() {
		// Convert mt19937_64 output to uniform [0,1) in a consistent way across platforms
		// mt19937_64 produces uint64_t values in [0, 2^64-1]
		// Divide by 2^64 to get [0,1)
		return static_cast<C>(mt()) / static_cast<C>(18446744073709551616.0);
	}

public:
	explicit AliasTable(const T& N) : N{N} {
		probs = new C[N];
		aliases = new T[N];
	}

	~AliasTable() {
		delete[] probs;
		delete[] aliases;
	}

	void initialize(const D* weights, const T& newSize) {
		N = newSize;
		probs = new C[N];
		aliases = new T[N];
		initialize(weights);
	}

	void initialize(const D* weights) {
		const D sm = std::accumulate(weights, weights + N, 0);
		for (T i = 0; i != N; i++) probs[i] = weights[i] * N / sm;
		queue<T> small = queue<T>();
		queue<T> large = queue<T>();
		for (T i = 0; i < N; i++) ((probs[i] < 1) ? small : large).push(i);
		while (! large.empty() && ! small.empty()) {
			T big = large.front();
			large.pop();
			T little = small.front();
			small.pop();
			aliases[little] = big;
			probs[big] = probs[big] + probs[little] - 1;
			(probs[big] < 1 ? small : large).push(big);
		}
		C accu = 0;
		while (! large.empty()) {
			accu += 1 - probs[large.front()];
			probs[large.front()] = 1;
			large.pop();
		}
		while (! small.empty()) {
			accu += 1 - probs[small.front()];
			probs[small.front()] = 1;
			small.pop();
		}
		if (accu > 1e-5) Rcpp::warning("Numerical instability in alias table " + to_string(accu));
	};

	long initRandom(long seed) {
		mt = mt19937_64(seed);
		return mt();
	}

	void initRandom() {
		random_device seed;
		initRandom(seed());
	}

	T operator()(const C& random, const C& random2) const {
		const T candidate = random * N;
		return (random2 >= probs[candidate]) ? aliases[candidate] : candidate;
	}

	T operator()() {
		return (*this)(uniformRandom(), uniformRandom());
	}
};
