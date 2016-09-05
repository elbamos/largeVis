#ifndef _LARGEVIS
#define _LARGEVIS
//#ifdef _WIN32
//#define ARMA_32BIT_WORD
//#endif
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "progress.hpp"
#include <math.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <set>
#include <memory>

using namespace Rcpp;
using namespace std;
using namespace arma;
/*
 * Global types
 */
typedef double distancetype;
typedef double coordinatetype;

typedef arma::sword vertexidxtype;
typedef arma::sword edgeidxtype;
typedef arma::sword iterationtype;

typedef int dimidxtype;
typedef int kidxtype;

#ifdef _OPENMP
void checkCRAN(Rcpp::Nullable<Rcpp::NumericVector> threads);
#endif

/*
 * Distance Functions
 */
distancetype dist(const arma::vec& i, const arma::vec& j);
distancetype relDist(const arma::vec& i, const arma::vec& j);
distancetype cosDist(const arma::vec& i, const arma::vec& j);
distancetype sparseDist(const arma::sp_mat& i, const arma::sp_mat& j);
distancetype sparseRelDist(const arma::sp_mat& i, const arma::sp_mat& j);
distancetype sparseCosDist(const arma::sp_mat& i, const arma::sp_mat& j);

// Exported distance functions for high dimensional space
arma::vec fastDistance(const NumericVector is,
                       const NumericVector js,
                       const arma::mat& data,
                       const std::string& distMethod,
                       bool verbose);
arma::vec fastSparseDistance(const arma::vec& is,
                             const arma::vec& js,
                             const arma::sp_mat& data,
                             const std::string& distMethod,
                             bool verbose);
arma::vec fastCDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& p_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose);
arma::vec fastSDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& j_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose);

class Gradient {
protected:
	const distancetype gamma;
	distancetype cap;
	const dimidxtype D;
	Gradient(const distancetype g,
          const dimidxtype d);
	virtual void _positiveGradient(const distancetype dist_squared,
                                coordinatetype* holder) const = 0;
	virtual void _negativeGradient(const distancetype dist_squared,
                                coordinatetype* holder) const = 0;
	inline void multModify(coordinatetype *col, coordinatetype adj) const;
	inline coordinatetype clamp(coordinatetype val) const;

public:
	virtual void positiveGradient(const coordinatetype* i,
                               const coordinatetype* j,
                               coordinatetype* holder) const;;
	virtual void negativeGradient(const coordinatetype* i,
                               const coordinatetype* k,
                               coordinatetype* holder) const;
	inline distancetype distAndVector(const coordinatetype *x_i,
                                   const coordinatetype *x_j,
                                   coordinatetype *output) const;
};

class AlphaGradient: public Gradient {
	const coordinatetype alpha;
	const coordinatetype twoalpha;
protected:
	const coordinatetype alphagamma;
	virtual void _positiveGradient(const double dist_squared,
                                coordinatetype* holder) const;
	virtual void _negativeGradient(const double dist_squared,
                                coordinatetype* holder) const;
public:
	AlphaGradient(const distancetype a,
               const distancetype g,
               const dimidxtype D);
};

class AlphaOneGradient: public AlphaGradient {
public:
	AlphaOneGradient(const distancetype g,
                  const dimidxtype d);
protected:
	virtual void _positiveGradient(const distancetype dist_squared,
                                coordinatetype* holder) const;
	virtual void _negativeGradient(const distancetype dist_squared,
                                coordinatetype* holder) const;
};

class ExpGradient: public Gradient {
public:
	const coordinatetype gammagamma;
	ExpGradient(const distancetype g, const dimidxtype d);
protected:
	virtual void _positiveGradient(const distancetype dist_squared,
                                coordinatetype* holder) const;
	virtual void _negativeGradient(const distancetype dist_squared,
                                coordinatetype* holder) const;
};


#endif
