#include "largeVis.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
bool checkBits() {
	if (sizeof( vertexidxtype ) == 8) return true;
	else return false;
}
