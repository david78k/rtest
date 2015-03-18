library("Rcpp")

Rcpp.package.skeleton("markovchainTaeSeungKang", example_code = FALSE,
			cpp_files = c("1_functions4Fitting.cpp"))

sourceCpp("compileAttributes.cpp")

