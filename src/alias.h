#include "largeVis.h"

template <class T, // number of aliases
          class C, // coordinate type, used for probabilities
          class D> // distance type, used to take weights on initialization
class AliasTable {
private:
	unique_ptr< C[] > probs;
	unique_ptr< T[] > aliases;
	uniform_real_distribution< C > rnd = uniform_real_distribution< C >();
	mt19937_64 mt;
	T N;

public:
	AliasTable() {
	}

	void initialize(const D* weights, T N) {
		this -> N = N;
		probs = unique_ptr< C[] >( new C[N] );
		aliases = unique_ptr< T[] >(new T[N]);
		D sm = 0;
		for (T i = 0; i != N; i++) sm += weights[i];
		for (T i = 0; i != N; i++) probs[i] = weights[i] * N / sm;
		queue<T> small = queue<T>();
		queue<T> large = queue<T>();
		for (T i = 0; i < N; i++) ((probs[i] < 1) ?
                               small :
                               large).push(i);
		while (! large.empty() & ! small.empty()) {
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
			accu += 1 - large.front();
			probs[large.front()] = 1;
			large.pop();
		}
		while (! small.empty()) {
			accu += 1 - small.front();
			probs[small.front()] = 1;
			small.pop();
		}
		if (accu > 1e-5) warning("Numerical instability in alias table " + to_string(accu));
	};

	T search(C random, C random2) const {
		T candidate = random * N;
		return (random2 >= probs[candidate]) ? aliases[candidate] : candidate;
	};

	long initRandom(long seed) {
		mt = mt19937_64(seed);
		return mt();
	}
	void initRandom() {
		random_device seed;
		initRandom(seed());
	}
	C getRand() {
		return rnd(mt);
	}
	T operator()() {
		return search(rnd(mt), rnd(mt));
	}
};