#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix createSequenceMatrix_cpp(int N) {
  NumericMatrix mat(N, 2);
  return (mat);
}

// .mcFitMle<-function(stringchar,byrow)
void _mcFitMle_cpp() {

}

// .mcFitLaplacianSmooth<-function(stringchar,byrow,laplacian=0.01)
void _mcFitLaplacianSmooth_cpp() {

}

// .bootstrapCharacterSequences<-function(stringchar, n, size=length(stringchar))
void _bootstrapCharacterSequences_cpp() {

}

// .fromBoot2Estimate<-function(listMatr)
void _fromBoot2Estimate_cpp() {

}

// .mcFitBootStrap<-function(data, nboot=10,byrow=TRUE, parallel=FALSE)
void _mcFitBootStrap_cpp() {

}

// .matr2Mc<-function(matrData,laplacian=0) 
void _matr2Mc_cpp() {

}


// markovchainFit<-function(data,method="mle", byrow=TRUE,nboot=10,laplacian=0, name, parallel=FALSE)
// [[Rcpp::export]]
NumericMatrix markovchainFit_cpp(int N, int thin) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;

  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      //x = rgamma(1, 3, 1 / (y * y + 4))[0];
      x = rgamma(1, 3, y * y + 4)[0];
      y = rnorm(1, 1 / (x + 1), 1 / sqrt(2 * (x + 1)))[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }

  return(mat);
}


//*** R 
/*
library(microbenchmark)
sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
microbenchmark(
  markovchainFit(data = sequence),
  markovchainFit_cpp(100, 10)
)
*/

