#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

NumericVector rowSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) 
      total += x(i, j);
    out[i] = total;
  }
  return out;
}

// createSequenceMatrix<-function(stringchar, toRowProbs=FALSE,sanitize=TRUE)
NumericMatrix createSequenceMatrix_cpp(CharacterVector stringchar, bool toRowProbs=false, bool sanitize=true) {
//NumericMatrix createSequenceMatrix_cpp(DataFrame stringchar, bool toRowProbs=false, bool sanitize=true) {
/* 
  elements <- sort(unique(stringchar))
  sizeMatr <- length(elements)
  freqMatrix <- zeros(sizeMatr)
  rownames(freqMatrix) <- elements
  colnames(freqMatrix) <- elements
  for(i in 1:(length(stringchar)-1))
  {
    posFrom <- which(rownames(freqMatrix)==stringchar[i])
    posTo <- which(rownames(freqMatrix)==stringchar[i+1])
    freqMatrix[posFrom,posTo]=freqMatrix[posFrom,posTo]+1
  }
  #sanitizing if any row in the matrix sums to zero by posing the corresponding diagonal equal to 1/dim
  if(sanitize==TRUE)
  {
	  if(any(rowSums(freqMatrix)==0))
	  {
		  indexesToBeSanitized<-which(rowSums(freqMatrix)==0)
		  for(i in indexesToBeSanitized) {
			  for(j in 1:sizeMatr) freqMatrix[i,j]<-1/sizeMatr
		  }
	  }
  }
  if(toRowProbs==TRUE)
  {
    freqMatrix<-freqMatrix/rowSums(freqMatrix)
  }
  return(freqMatrix)
*/
//  Rcout << stringchar[0] << endl;
  //Rf_PrintValue(stringchar);
//  [1] "a" "b" "a" "a" "a" "a" "b" "a" "b" "a" "b" "a" "a" "b" "b" "b" "a"
  
  CharacterVector elements = unique(stringchar).sort();
  //Rf_PrintValue(elements);
  int sizeMatr = elements.size();
  //Rcout << sizeMatr << endl;
  
  NumericMatrix freqMatrix(sizeMatr);
  rownames(freqMatrix) = elements;
  colnames(freqMatrix) = elements;
  //Rf_PrintValue(rownames(freqMatrix));
  //Rf_PrintValue(colnames(freqMatrix));
  //Rf_PrintValue(freqMatrix.names());
  CharacterVector rnames = rownames(freqMatrix);
//  Rf_PrintValue(freqMatrix);
  for(int i = 0; i < stringchar.size() - 1; i ++) {
    int posFrom = find(rnames.begin(), rnames.end(), stringchar[i]) - rnames.begin();
    int posTo = find(rnames.begin(), rnames.end(), stringchar[i + 1]) - rnames.begin();
    //Rcout << stringchar[i] << "->" << stringchar[i + 1] << ": " << posFrom << " " << posTo << endl;
    freqMatrix(posFrom,posTo)++;
    //freqMatrix[posFrom][posTo]=freqMatrix[posFrom][posTo]+1;
    //Rf_PrintValue(freqMatrix);
  }
 
  //Rf_PrintValue(freqMatrix);

  //sanitizing if any row in the matrix sums to zero by posing the corresponding diagonal equal to 1/dim
  if(sanitize==true)
  {
	for (int i = 0; i < sizeMatr; i++) {
    		double total = 0;
    		for (int j = 0; j < sizeMatr; j++) 
      			total += freqMatrix(i, j);
		if(total == 0)
    			for (int j = 0; j < sizeMatr; j++) 
      				freqMatrix(i, j) = 1/sizeMatr;
	}
  }
  if(toRowProbs==true)
  {
    //freqMatrix<-freqMatrix/rowSums(freqMatrix)
	for (int i = 0; i < sizeMatr; i++) {
    		double total = 0;
    		for (int j = 0; j < sizeMatr; j++) 
      			total += freqMatrix(i, j);
    		for (int j = 0; j < sizeMatr; j++) 
      			freqMatrix(i, j) /= total;
	}
  }

  return (freqMatrix);
}

//template <typename T>
//T transpose(const T & m) {      // tranpose for IntegerMatrix / NumericMatrix, see array.c in R
NumericMatrix transpose(NumericMatrix & m) {      // tranpose for IntegerMatrix / NumericMatrix, see array.c in R
  int k = m.rows(), n = m.cols();
  //Rcpp::Rcout << "Transposing " << n << " by " << k << std::endl;
  //T z(n, k);
  NumericMatrix z(n, k);
  //CharacterVector rows = rownames(m);
  //CharacterVector cols = colnames(m);
  //List dimnms = List::create(CharacterVector::create(rownames(m), colnames(m));
  z.attr("dimnames") = List::create(CharacterVector::create(rownames(m), colnames(m))); 
  //rownames(z) = rownames(m);
  //Rf_PrintValue(rows);
  //Rf_PrintValue(rownames(m));
  //Rf_PrintValue(rownames(z));
  int sz1 = n*k-1;
  //typename T::iterator mit, zit;
  NumericMatrix::iterator mit, zit;
  for (mit = m.begin(), zit = z.begin(); mit != m.end(); mit++, zit += n) {
  	if (zit >= z.end()) zit -= sz1;
        *zit = *mit;
  }
  return(z);
}

// .mcFitMle<-function(stringchar,byrow)
NumericMatrix _mcFitMle(CharacterVector stringchar, bool byrow) {
//NumericMatrix _mcFitMle(DataFrame stringchar, bool byrow) {
/*
  initialMatr<-createSequenceMatrix(stringchar=stringchar,toRowProbs=TRUE)
  outMc<-new("markovchain", transitionMatrix=initialMatr,name="MLE Fit")
  if(byrow==FALSE) outMc<-t(outMc)
  out<-list(estimate=outMc)
*/
  NumericMatrix out = createSequenceMatrix_cpp(stringchar, true);
  NumericMatrix initialMatr = createSequenceMatrix_cpp(stringchar, true);
  Rf_PrintValue(out);
  //rownames(out) = rownames(initialMatr);
  //Rf_PrintValue(rownames(initialMatr));
  //Rf_PrintValue(rownames(out));
  out = transpose(out);
  //rownames(out) = rownames(initialMatr);
  //colnames(out) = colnames(initialMatr);
  //Rf_PrintValue(rownames(initialMatr));
  //Rf_PrintValue(rownames(out));
  //Rf_PrintValue(out);
  //NumericMatrix outMc(initialMatr); //("markovchain", initialMatr,"MLE Fit");
//  if(byrow==false) outMc = transpose(outMc);
  //out = list(outMc);
  return out;
}

// .mcFitLaplacianSmooth<-function(stringchar,byrow,laplacian=0.01)
NumericMatrix _mcFitLaplacianSmooth(DataFrame data, bool byrow, double laplacian=0.01) {
  NumericMatrix mat(1, 1);
/*
  origNum<-createSequenceMatrix(stringchar=stringchar,toRowProbs=FALSE)
	sumOfRow<-rowSums(origNum)
	origDen<-matrix(rep(sumOfRow,length(sumOfRow)),byrow = FALSE,ncol=length(sumOfRow))
	newNum<-origNum+laplacian
	newSumOfRow<-rowSums(newNum)
	newDen<-matrix(rep(newSumOfRow,length(newSumOfRow)),byrow = FALSE,ncol=length(newSumOfRow))
	transMatr<-newNum/newDen
	outMc<-new("markovchain", transitionMatrix=transMatr,name="Laplacian Smooth Fit")
	if(byrow==FALSE) outMc<-t(outMc)
	out<-list(estimate=outMc)
*/
  return mat;
}

// .bootstrapCharacterSequences<-function(stringchar, n, size=length(stringchar))
void _bootstrapCharacterSequences() {
/*
  contingencyMatrix<-createSequenceMatrix(stringchar=stringchar)
  samples<-list()
  itemset<-rownames(contingencyMatrix)
  for(i in 1:n) #cicle to fill the samples
  {
    charseq<-character()
    char<-sample(x=itemset,size=1)
    charseq<-c(charseq,char)
    for(j in 2:size) #cicle to define the element in a list
    {
      probsVector<-contingencyMatrix[which(rownames(contingencyMatrix)==char),]
      char<-sample(x=itemset,size=1, replace=TRUE,prob=probsVector)
      charseq<-c(charseq,char)
    }
    samples[[length(samples)+1]]<-charseq #increase the list
  }
  return(samples)
*/
}

// .fromBoot2Estimate<-function(listMatr)
void _fromBoot2Estimate() {
/*
  sampleSize<-length(listMatr)
  matrDim<-nrow(listMatr[[1]])
  #create the estimate output
  matrMean<-zeros(matrDim)
  matrSd<-zeros(matrDim)
  #create the sample output
  for(i in 1:matrDim) #move by row
  {
    for(j in 1:matrDim) #move by cols
    {
      probsEstimated<-numeric()
      #fill the probs
      for(k in 1:sampleSize) probsEstimated<-c(probsEstimated,listMatr[[k]][i,j])
      muCell<-mean(probsEstimated)
      sdCell<-sd(probsEstimated)
      matrMean[i,j]<-muCell
      matrSd[i,j]<-sdCell
    }
  }
	out<-list(estMu=matrMean, estSigma=matrSd)
    return(out)
*/
}

// .mcFitBootStrap<-function(data, nboot=10,byrow=TRUE, parallel=FALSE)
NumericMatrix _mcFitBootStrap(DataFrame data, int nboot=10, bool byrow=true, bool parallel=false) {
  NumericMatrix mat(1, 1);
/*
  #create the list of bootstrap sequence sample
	theList<-.bootstrapCharacterSequences(stringchar=data, n=nboot)
	if(!parallel)
		#convert the list in a probability matrix
		pmsBootStrapped<-lapply(X=theList, FUN=createSequenceMatrix, toRowProbs=TRUE,sanitize=TRUE)
		 else {
		#require(parallel)
		type <- if(exists("mcfork", mode = "function")) "FORK" else "PSOCK"
		cores <- getOption("mc.cores", detectCores())
		cl <- makeCluster(cores, type = type)
			clusterExport(cl, varlist = c("createSequenceMatrix","zeros"))
			pmsBootStrapped<-parLapply(cl=cl, X=theList, fun="createSequenceMatrix", toRowProbs=TRUE,sanitize=TRUE)
		stopCluster(cl)
	}
 
  estimateList<-.fromBoot2Estimate(listMatr=pmsBootStrapped)
  #from raw to estimate
  temp<-estimateList$estMu
  transMatr<-sweep(temp, 1, rowSums(temp), FUN="/")
  estimate<-new("markovchain",transitionMatrix=transMatr, byrow=byrow, name="BootStrap Estimate")
  out<-list(estimate=estimate, standardError=estimateList$estSigma,bootStrapSamples=pmsBootStrapped)
  return(out)
*/
  return mat;
}

// .matr2Mc<-function(matrData,laplacian=0) 
NumericMatrix _matr2Mc(NumericMatrix matrData, double laplacian=0) {
/*
  #find unique values scanning the matrix
  nCols<-ncol(matrData)
  uniqueVals<-character()
  for(i in 1:nCols) uniqueVals<-union(uniqueVals,unique(as.character(matrData[,i])))
  uniqueVals<-sort(uniqueVals)
  #create a contingency matrix
  contingencyMatrix<-matrix(rep(0,length(uniqueVals)^2),ncol=length(uniqueVals))
  rownames(contingencyMatrix)<-colnames(contingencyMatrix)<-uniqueVals
  #fill the contingency matrix
  for(i in 1:nrow(matrData))
  {
    for( j in 2:nCols)
    {
      stateBegin<-as.character(matrData[i,j-1]);whichRow<-which(uniqueVals==stateBegin)
      stateEnd<-as.character(matrData[i,j]);whichCols<-which(uniqueVals==stateEnd)
      contingencyMatrix[whichRow,whichCols]<-contingencyMatrix[whichRow,whichCols]+1
    }
  }
  #add laplacian correction if needed
  contingencyMatrix<-contingencyMatrix+laplacian
  #get a transition matrix and a DTMC
  transitionMatrix<-contingencyMatrix/rowSums(contingencyMatrix)
  outMc<-new("markovchain",transitionMatrix=transitionMatrix)
*/
  NumericMatrix outMc;
 
  return(outMc);
}


// markovchainFit<-function(data,method="mle", byrow=TRUE,nboot=10,laplacian=0, name, parallel=FALSE)
// [[Rcpp::export]]
NumericMatrix markovchainFit_cpp(CharacterVector data, String method="mle", bool byrow=true, int nboot=10, double laplacian=0, String name="", bool parallel=false) {
//NumericMatrix markovchainFit_cpp(DataFrame data, String method="mle", bool byrow=true, int nboot=10, double laplacian=0, String name="", bool parallel=false) {
  NumericMatrix out;
  //if(class(data) %in% c("data.frame","matrix")) {
  if(data.attr("class") == "data.frame" || data.attr("class") == "matrix") {
    //#if data is a data.frame forced to matrix
    //if(data.attr("class") == "data.frame") data =as.matrix(data);
    CharacterMatrix data2;
    if(data.attr("class") == "data.frame") data2 =wrap(data);
    //if(data.attr("class") == "data.frame") data =as.matrix(data);
    //byrow assumes distinct observations (trajectiories) are per row
    //otherwise transpose
    if(!byrow) data = trans(data2);
    NumericMatrix outMc =_matr2Mc(data,laplacian);
    //out<-list(estimate=outMc)
  } else {
    if(method == "mle") out = _mcFitMle(data, byrow);
    if(method == "bootstrap") out = _mcFitBootStrap(data, nboot, byrow, parallel);
    if(method == "laplace") out = _mcFitLaplacianSmooth(data, byrow, laplacian);
  }
  return out;
}

NumericMatrix simplemc(int N, int thin) {
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


/*** R 
library(microbenchmark)
sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a", "b")
  markovchainFit_cpp(sequence)
*/
/*
microbenchmark(
  markovchainFit(data = sequence),
  markovchainFit_cpp(100, 10)
)
*/

