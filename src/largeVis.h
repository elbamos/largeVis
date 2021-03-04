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

#endif
