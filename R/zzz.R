.onLoad <- function(libname, pkgname) {
	if (! largeVis:::checkBits()) warning("largeVis was compiled with 32-bit types. This will limit the size of the datasets it can process. Consider recompiling with -DARMA_64BIT_WORD")
	if (! largeVis:::checkOpenMP()) message("largeVis was compiled without OpenMP support.")
}