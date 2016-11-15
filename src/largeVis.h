#ifndef _LARGEVIS
#define _LARGEVIS
//#ifdef _WIN32
//#define ARMA_32BIT_WORD
//#endif
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;
/*
 * Global types
 */
typedef double distancetype;
typedef double coordinatetype;

typedef sword vertexidxtype;
typedef sword edgeidxtype;
typedef uword iterationtype;

typedef unsigned int dimidxtype;
typedef unsigned int kidxtype;

#ifdef _OPENMP
void checkCRAN(Rcpp::Nullable<Rcpp::NumericVector> threads);
#endif

#endif
