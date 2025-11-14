#include <Rcpp.h>
using namespace Rcpp;

long countSize(IntegerMatrix x) {
	long ones = 0;
	for (long i = 0; i < x.ncol(); ++i) {
		for (long j = 0; j < x.nrow(); ++j) {
			if (x(j, i) == -1) ones++;
		}
	}
	return (x.ncol() * x.nrow()) - ones;
}

// [[Rcpp::export]]
List neighborsToVectors(IntegerMatrix x) {
	long cols = x.ncol();
	long rows = x.nrow();
	long sz = countSize(x);
	IntegerVector is(sz);
	IntegerVector js(sz);
	long pos = 0;
	for (long col = 0; col < cols; ++col) {
		for (long row = 0; row < rows; ++row) {
			if (x(row, col) == -1) break;
			is[pos] = col;
			js[pos++] = x(row, col);
		}
	}
	return List::create(Named("i") = is , Named("j") = js);
}
