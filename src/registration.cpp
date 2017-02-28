#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls *//*
extern SEXP largeVis_checkBits();
extern SEXP largeVis_checkOpenMP();
extern SEXP largeVis_dbscan_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_fastCDistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_fastDistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_fastSDistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_hdbscanc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_optics_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_referenceWij(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_searchTrees(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_searchTreesCSparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_searchTreesTSparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP largeVis_sgd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"largeVis_checkBits",          (DL_FUNC) &largeVis_checkBits,           0},
  {"largeVis_checkOpenMP",        (DL_FUNC) &largeVis_checkOpenMP,         0},
  {"largeVis_dbscan_cpp",         (DL_FUNC) &largeVis_dbscan_cpp,          5},
  {"largeVis_fastCDistance",      (DL_FUNC) &largeVis_fastCDistance,       8},
  {"largeVis_fastDistance",       (DL_FUNC) &largeVis_fastDistance,        6},
  {"largeVis_fastSDistance",      (DL_FUNC) &largeVis_fastSDistance,       8},
  {"largeVis_hdbscanc",           (DL_FUNC) &largeVis_hdbscanc,            6},
  {"largeVis_optics_cpp",         (DL_FUNC) &largeVis_optics_cpp,          6},
  {"largeVis_referenceWij",       (DL_FUNC) &largeVis_referenceWij,        5},
  {"largeVis_searchTrees",        (DL_FUNC) &largeVis_searchTrees,         9},
  {"largeVis_searchTreesCSparse", (DL_FUNC) &largeVis_searchTreesCSparse, 11},
  {"largeVis_searchTreesTSparse", (DL_FUNC) &largeVis_searchTreesTSparse, 11},
  {"largeVis_sgd",                (DL_FUNC) &largeVis_sgd,                15},
  {NULL, NULL, 0}
};*/

void R_init_largeVis(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, /*CallEntries*/ NULL, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}