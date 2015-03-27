// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
//using namespace RcppParallel;

#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// [[Rcpp::export]]
List ForLoop(List theList) {

  // allocate the matrix we will return
  //NumericMatrix mat(orig.nrow(), orig.ncol());
  int n = theList.size();
  List pmsBootStrapped(n);

  // transform it
  //std::transform(orig.begin(), orig.end(), mat.begin(), ::sqrt);
        for(int i = 0; i < n; i ++) {
                //Rcout << omp_get_thread_num() << " cores: " << omp_get_num_threads() << std::endl;
                //pmsBootStrapped[i] = createSequenceMatrix_cpp(theList[i], true, true);
		pmsBootStrapped[i] = theList[i];
        }

  // return the new matrix
  return pmsBootStrapped;
}


struct ForLoopWorker : public RcppParallel::Worker
{
   // source matrix
   //const RcppParallel::RMatrix<double> input;
   const List input;
   
   // destination matrix
   //RcppParallel::RMatrix<double> output;
   List output;
   
   // initialize with source and destination
   //ForLoopWorker(const NumericMatrix input, NumericMatrix output) 
   ForLoopWorker(const List input, List output) 
      : input(input), output(output) {}
   
/*
 #pragma omp parallel for
        for(int i = 0; i < n; i ++) {
                Rcout << omp_get_thread_num() << " cores: " << omp_get_num_threads() << std::endl;
                pmsBootStrapped[i] = createSequenceMatrix_cpp(theList[i], true, true);
        }
*/
   // take the square root of the range of elements requested
   void operator()(std::size_t begin, std::size_t end) {
   //void operator()(int i) {
	//Rcout << i << std::endl;
//	printf("i %d\n", i);
//	List x = clone(input);
	//int i = 0;
	//for( List::iterator it = (x.begin() + begin); it != (x.end() + end); ++it ) {
		output[begin] = input[begin];		
		//output[i] = x[i];		
	//	i ++;
//	}
/*
      std::transform(input.begin() + begin, 
                     input.begin() + end, 
                     output.begin() + begin, 
                     ::sqrt);
*/	
   }
};

// [[Rcpp::export]]
List parallelForLoop(List theList) {
  
  // allocate the output matrix
 // NumericMatrix output(x.nrow(), x.ncol());
  int n = theList.size();
  List pmsBootStrapped(n);
  
  // SquareRoot functor (pass input and output matrixes)
  ForLoopWorker forloop(theList, pmsBootStrapped);
  
  // call parallelFor to do the work
  parallelFor(0, n, forloop);
  //parallelFor(0, theList.length(), forloop);
  
  // return the output matrix
  return pmsBootStrapped;
}


/*************************************************/
// [[Rcpp::export]]
NumericMatrix matrixSqrt(NumericMatrix orig) {

  // allocate the matrix we will return
  NumericMatrix mat(orig.nrow(), orig.ncol());

  // transform it
  std::transform(orig.begin(), orig.end(), mat.begin(), ::sqrt);

  // return the new matrix
  return mat;
}

struct SquareRoot : public RcppParallel::Worker
{
   // source matrix
   const RcppParallel::RMatrix<double> input;
   
   // destination matrix
   RcppParallel::RMatrix<double> output;
   
   // initialize with source and destination
   SquareRoot(const NumericMatrix input, NumericMatrix output) 
      : input(input), output(output) {}
   
   // take the square root of the range of elements requested
   void operator()(std::size_t begin, std::size_t end) {
      std::transform(input.begin() + begin, 
                     input.begin() + end, 
                     output.begin() + begin, 
                     ::sqrt);
   }
};

// [[Rcpp::export]]
NumericMatrix parallelMatrixSqrt(NumericMatrix x) {
  
  // allocate the output matrix
  NumericMatrix output(x.nrow(), x.ncol());
  
  // SquareRoot functor (pass input and output matrixes)
  SquareRoot squareRoot(x, output);
  
  // call parallelFor to do the work
  parallelFor(0, x.length(), squareRoot);
  
  // return the output matrix
  return output;
}


/*** R
dat <- list( 
    1:50000, ## integer
    as.numeric(1:50000) ## numeric
)
#print(dat)

tmp <- ForLoop(dat)
#print(tmp)

tmp <- parallelForLoop(dat)
#print(tmp)

# allocate a matrix
m <- matrix(as.numeric(c(1:1000000)), nrow = 1000, ncol = 1000)

# ensure that serial and parallel versions give the same result
#stopifnot(identical(ForLoop(m), parallelForLoop(m)))
#stopifnot(identical(matrixSqrt(m), parallelMatrixSqrt(m)))

# compare performance of serial and parallel
library(rbenchmark)
res <- benchmark(ForLoop(dat),
                 parallelForLoop(dat),
                 order="relative")

#res <- benchmark(matrixSqrt(m),
#                 parallelMatrixSqrt(m),
#                 order="relative")
res[,1:4]
*/
