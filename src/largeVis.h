#ifndef _LARGEVIS
#define _LARGEVIS

#ifndef __llvm__
#ifndef __clang__
#if ((__GNUC__ < 5) && (__GNUC_MINOR__ < 9))
#error largeVis is incompatible with gcc < 4.9. Upgrade gcc or use llvm.
#endif
#endif
#endif

#include <RcppArmadillo.h>

/*
 * Global types
 */
typedef double distancetype;
typedef double coordinatetype;

typedef arma::sword vertexidxtype;
typedef arma::sword edgeidxtype;
typedef arma::uword iterationtype;

typedef unsigned int dimidxtype;
typedef unsigned int kidxtype;

#endif
