.onAttach<- function(libname, pkgname) {
	if (! largeVis:::checkBits()) packageStartupMessage("largeVis was compiled with 32-bit types. This will limit the size of the datasets it can process. Consider recompiling with -DARMA_64BIT_WORD")
	if (! largeVis:::checkOpenMP()) packageStartupMessage("largeVis was compiled without OpenMP support.")
}