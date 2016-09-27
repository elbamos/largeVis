// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "gradients.h"

using namespace Rcpp;
using namespace std;
using namespace arma;


/*
* Gradients
*/
Gradient::Gradient(const distancetype& g,
          const dimidxtype& d) : gamma{g}, cap(5), D{d} {};

inline distancetype Gradient::distAndVector(const coordinatetype *x_i,
                                            const coordinatetype *x_j,
                                            coordinatetype *output) const {
	double cnt = 0;
	for (int d = 0; d < D; d++) {
		double t = x_i[d] - x_j[d];
		output[d] = t;
		cnt += t * t;
	}
	return cnt;
}

Gradient::~Gradient() {

}

inline void Gradient::multModify(coordinatetype *col, const coordinatetype& adj) const {
	for (dimidxtype i = 0; i != D; i++) col[i] = clamp(col[i] * adj);
}
inline void Gradient::multModifyPos(coordinatetype *col, const coordinatetype& adj) const {
	for (dimidxtype i = 0; i != D; i++) col[i] *= adj;
}
inline coordinatetype Gradient::clamp(const coordinatetype& val) const {
	return fmin(fmax(val, -cap), cap);
}

void Gradient::positiveGradient(const coordinatetype* i,
                             		const coordinatetype* j,
                             		coordinatetype* holder) const {
	const double dist_squared = distAndVector(i, j, holder);
	_positiveGradient(dist_squared, holder);
};

void Gradient::negativeGradient(const coordinatetype* i,
                             		const coordinatetype* k,
                             		coordinatetype* holder) const {
	const double dist_squared = distAndVector(i, k, holder) ;
	_negativeGradient(dist_squared, holder);
}


void AlphaGradient::_positiveGradient(const double& dist_squared,
                                      coordinatetype* holder) const {
	const distancetype grad = twoalpha / (1 + alpha * dist_squared);
	multModifyPos(holder, grad);
};
void AlphaGradient::_negativeGradient(const double& dist_squared,
                                      coordinatetype* holder) const {
	const distancetype adk = alpha * dist_squared;
	const distancetype grad = alphagamma / (dist_squared * (adk + 1));
	multModify(holder, grad);
};

AlphaGradient::AlphaGradient(const distancetype& a,
             const distancetype& g,
             const dimidxtype& D) : Gradient(g, D),
             alpha{a},
             twoalpha(alpha * -2),
             alphagamma(alpha * gamma * 2) { } ;

AlphaOneGradient::AlphaOneGradient(const distancetype& g,
                									 const dimidxtype& d) : AlphaGradient(1, g, d) {
}
void AlphaOneGradient::_positiveGradient(const distancetype& dist_squared,
                                         coordinatetype* holder) const {
	const distancetype grad = - 2 / (1 + dist_squared);
	multModifyPos(holder, grad);
}
void AlphaOneGradient::_negativeGradient(const distancetype& dist_squared,
                                         coordinatetype* holder) const {
	const distancetype grad = alphagamma / (1 + dist_squared) / (0.1 + dist_squared);
	multModify(holder, grad);
};

ExpGradient::ExpGradient(const distancetype& g,
                         const dimidxtype& d) : Gradient(g, d), gammagamma(gamma * gamma) {
		cap = gamma;
};
void ExpGradient::_positiveGradient(const distancetype& dist_squared,
                                    coordinatetype* holder) const {
	const distancetype expsq = exp(dist_squared);
	const distancetype grad = (dist_squared > 4) ? -1 :
		-(expsq / (expsq + 1));
	multModifyPos(holder, grad);
};

void ExpGradient::_negativeGradient(const distancetype& dist_squared,
                                    coordinatetype* holder) const {
	const distancetype grad = (dist_squared > gammagamma) ? 0 :
	gamma / (1 + exp(dist_squared));
	multModify(holder, grad);
};
