#include <Rcpp.h>
#include <omp.h>
#include <unistd.h>
//#include <Environment.h>
//#include <RcppArmadillo.h>
//// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

template <typename T>
T _transpose(T & m) {      // tranpose for IntegerMatrix / NumericMatrix, see array.c in R
//T transpose(const T & m) {      // tranpose for IntegerMatrix / NumericMatrix, see array.c in R
//NumericMatrix transpose(NumericMatrix & m) {      // tranpose for IntegerMatrix / NumericMatrix, see array.c in R
  int k = m.rows(), n = m.cols();
  //Rcpp::Rcout << "Transposing " << n << " by " << k << std::endl;
  T z(n, k);
  //NumericMatrix z(n, k);
  z.attr("dimnames") = List::create(colnames(m), rownames(m)); 
  //z.attr("dimnames") = List::create(rownames(m), colnames(m)); 
  int sz1 = n*k-1;
  //NumericMatrix::iterator mit, zit;
  typename T::iterator mit, zit;
  for (mit = m.begin(), zit = z.begin(); mit != m.end(); mit++, zit += n) {
  	if (zit >= z.end()) zit -= sz1;
        *zit = *mit;
  }
  return(z);
}

NumericMatrix _toRowProbs(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix out(nrow);

  for (int i = 0; i < nrow; i++) {
    double rowSum = 0;
    for (int j = 0; j < ncol; j++) 
      rowSum += x(i, j);
    for (int j = 0; j < ncol; j++) 
      out(i, j) = x(i, j)/rowSum;
  }
  out.attr("dimnames") = List::create(rownames(x), colnames(x)); 
  return out;
}

NumericVector _rowSumsC(NumericMatrix x) {
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
// [[Rcpp::export]]
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
  freqMatrix.attr("dimnames") = List::create(elements, elements); 
//  Rf_PrintValue(rownames(freqMatrix));
//  Rf_PrintValue(colnames(freqMatrix));
  //Rf_PrintValue(freqMatrix.names());
  CharacterVector rnames = rownames(freqMatrix);
//  Rf_PrintValue(freqMatrix);
  int posFrom, posTo;
  for(int i = 0; i < stringchar.size() - 1; i ++) {
	for (int j = 0; j < rnames.size(); j ++) {
		if(stringchar[i] == rnames[j]) posFrom = j;
		if(stringchar[i + 1] == rnames[j]) posTo = j;
	}
    //int posFrom = find(rnames.begin(), rnames.end(), stringchar[i]) - rnames.begin();
    //int posTo = find(rnames.begin(), rnames.end(), stringchar[i + 1]) - rnames.begin();
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
    		double rowSum = 0;
    		for (int j = 0; j < sizeMatr; j++) 
      			rowSum += freqMatrix(i, j);
		if(rowSum == 0)
    			for (int j = 0; j < sizeMatr; j++) 
      				freqMatrix(i, j) = 1/sizeMatr;
	}
  }
  if(toRowProbs==true)
  {
    //freqMatrix<-freqMatrix/rowSums(freqMatrix)
	freqMatrix = _toRowProbs(freqMatrix);
/*	for (int i = 0; i < sizeMatr; i++) {
    		double rowSum = 0;
    		for (int j = 0; j < sizeMatr; j++) 
      			rowSum += freqMatrix(i, j);
    		for (int j = 0; j < sizeMatr; j++) 
      			freqMatrix(i, j) /= rowSum;
	}
*/
  }

  return (freqMatrix);
}

// .mcFitMle<-function(stringchar,byrow)
List _mcFitMle(CharacterVector stringchar, bool byrow) {
//NumericMatrix _mcFitMle(DataFrame stringchar, bool byrow) {
/*
  initialMatr<-createSequenceMatrix(stringchar=stringchar,toRowProbs=TRUE)
  outMc<-new("markovchain", transitionMatrix=initialMatr,name="MLE Fit")
  if(byrow==FALSE) outMc<-t(outMc)
  out<-list(estimate=outMc)
*/
  //NumericMatrix out = createSequenceMatrix_cpp(stringchar, true);
  NumericMatrix initialMatr = createSequenceMatrix_cpp(stringchar, true);
//  Rf_PrintValue(initialMatr);
  
  //if(byrow==false) outMc = transpose(outMc);
  if(byrow==false) initialMatr = _transpose(initialMatr);

  //NumericMatrix outMc(initialMatr); //("markovchain", initialMatr,"MLE Fit");
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = initialMatr;
  outMc.slot("name") = "MLE Fit";  

  return List::create(_["estimate"] = outMc);
  //List out(1);
  //out[0] = outMc;
  //return out;
}

// .mcFitLaplacianSmooth<-function(stringchar,byrow,laplacian=0.01)
List _mcFitLaplacianSmooth(CharacterVector stringchar, bool byrow, double laplacian=0.01) {
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
  NumericMatrix origNum = createSequenceMatrix_cpp(stringchar, false);
  int nRows = origNum.nrow(), nCols = origNum.ncol();
  for(int i = 0; i < nRows; i ++) {
	double rowSum = 0;
	for(int j = 0; j < nCols; j ++) {
    		origNum(i,j) += laplacian;
    		rowSum += origNum(i,j);
  	}
  	//#get a transition matrix and a DTMC
	for(int j = 0; j < nCols; j ++) 
    		origNum(i,j) /= rowSum;
  }
  
  if(byrow==false) origNum = _transpose(origNum);
//  Rf_PrintValue(origNum);
 
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = origNum;
  outMc.slot("name") = "Laplacian Smooth Fit";  

  return List::create(_["estimate"] = outMc);
}

// .bootstrapCharacterSequences<-function(stringchar, n, size=length(stringchar))
//List _bootstrapCharacterSequences(CharacterVector stringchar, int n, int size=stringchar.size()) {
List _bootstrapCharacterSequences(CharacterVector stringchar, int n, int size=-1) {
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
  if(size == -1) size = stringchar.size();
  NumericMatrix contingencyMatrix = createSequenceMatrix_cpp(stringchar);
  List samples;
  CharacterVector itemset = rownames(contingencyMatrix);
  int itemsetsize = itemset.size();
  //Rf_PrintValue(itemset);

  Function sample("sample");
  srand(time(NULL));
  for(int i = 0; i < n; i ++) {
	CharacterVector charseq;	
	int rnd = rand()%itemsetsize;
 	String ch = itemset[rnd];
	//Rcout << rnd << " " << itemset[rnd] << endl;
	//Rf_PrintValue(ch);
	List res = sample(itemset, 1);
	CharacterVector cv = res[0];
	charseq.push_back(cv[0]);
	//charseq.push_back(ch);
	//Rf_PrintValue(charseq);
	for(int j = 1; j < size; j ++) {
		NumericVector probsVector;
		for(int k = 0; k < itemsetsize; k ++) 
			if((string)itemset[k] == (string) ch) {
				probsVector = contingencyMatrix(k, _);	
				//Rcout << k << " " << (string)ch << endl;
				//Rf_PrintValue(probsVector);
				break;
			}
  		//srand(time(NULL));
		rnd = rand()%itemsetsize;
 		ch = itemset[rnd];
		res = sample(itemset, 1, true, probsVector);
		//SEXP character = sample(itemset, 1, true, probsVector);
		//Rcout << res[0] << endl;
		//Rcout << "res[0]" << endl;
		//Rf_PrintValue(res[0]);
		CharacterVector v = res[0];
		//Rf_PrintValue(res[0]);
		//Rf_PrintValue(character);
		charseq.push_back(v[0]);
		//charseq.push_back(ch);
 	}
	//samples[[samples.size() + 1]] = charseq;
	samples.push_back(charseq);
	//samples = List::create(clone(samples), charseq);
  }

  return samples;
}

// .fromBoot2Estimate<-function(listMatr)
List _fromBoot2Estimate(List listMatr) {
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
  int sampleSize = listMatr.size();
  NumericMatrix firstMat = listMatr[0];
  int matrDim = firstMat.nrow();
  //Rcout << matrDim << endl;
  NumericMatrix matrMean(matrDim);
  NumericMatrix matrSd(matrDim);

  for(int i = 0; i < matrDim; i ++) { 
  	for(int j = 0; j < matrDim; j ++) { 
		NumericVector probsEstimated;
		for(int k = 0; k < sampleSize; k ++) {
			NumericMatrix mat = listMatr[k];
			probsEstimated.push_back(mat(i,j));
			matrMean(i,j) = mean(probsEstimated);
			matrSd(i,j) = sd(probsEstimated);
		}
  	}
  }
  return List::create(_["estMu"]=matrMean, _["estSigma"]=matrSd);
}

List lapply(List input, Function f) {
  int n = input.size();
  List out(n);

  for(int i = 0; i < n; i++) {
    out[i] = f(input[i]);
  }

  return out;
}

// .mcFitBootStrap<-function(data, nboot=10,byrow=TRUE, parallel=FALSE)
List _mcFitBootStrap(CharacterVector data, int nboot=10, bool byrow=true, bool parallel=false) {
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
  List theList = _bootstrapCharacterSequences(data, nboot);
  int n = theList.size();
  //Rcout << "theList.size() " << n << endl;
  List pmsBootStrapped(theList.size());
  //List pmsBootStrapped(theList.size());

  //if(!parallel) pmsBootStrapped = theList;
  if(!parallel) { //pmsBootStrapped = lapply(theList, createSequenceMatrix_cpp, true, true);
	for(int i = 0; i < n; i++) { 
		pmsBootStrapped[i] = createSequenceMatrix_cpp(theList[i], true, true);
//		Rf_PrintValue(pmsBootStrapped[i]);
	}
  } else {
	//String type;
//	if(Rf_findFun("mcfork", R_GlobalEnv)) 
	//if(exists("mcfork"))
	//SEXP nameSym = Rf_install("mcfork");
	//SEXP res = Rf_findVarInFrame( String::get__() , nameSym  ) ;
	//if(res != R_UnboundValue)
		//type = "FORK"; else "PSOCK";			
	int cores = sysconf(_SC_NPROCESSORS_ONLN);
//	Rcout << "cores: " << cores << endl;
	//omp_set_dynamic(0);
	omp_set_num_threads(cores);
	#pragma omp parallel for
	for(int i = 0; i < n; i ++) {
		pmsBootStrapped[i] = createSequenceMatrix_cpp(theList[i], true, true);
		//Rf_PrintValue(pmsBootStrapped[i]);
	}
  }
  //estimateList<-.fromBoot2Estimate(listMatr=pmsBootStrapped)
  List estimateList = _fromBoot2Estimate(pmsBootStrapped);
  //NumericMatrix temp = estimateList["estMu"];
  NumericMatrix transMatr = _toRowProbs(estimateList["estMu"]);
  //int size = temp.nrow();

  //NumericMatrix transMatr = _toRowProbs(temp);
  //NumericMatrix transMatr(size);//= sweep(temp, 1, rowSumsC(temp), FUN="/");
/*
  for(int i = 0; i < size; i ++) {
	double rowSum = 0;
  	for(int j = 0; j < size; j ++) 
		rowSum += temp(i, j);
  	for(int j = 0; j < size; j ++) 
		temp(i, j) /= rowSum;
  }
  transMatr = temp;
*/
  S4 estimate("markovchain");
  estimate.slot("transitionMatrix") = transMatr;
  estimate.slot("byrow") = byrow;
  estimate.slot("name") = "BootStrap Estimate";  

  //Rcout << "bootStrapSamples.size() " << pmsBootStrapped.size() << endl;

  return List::create(_["estimate"] = estimate
		, _["standardError"] = estimateList["estSigma"]
		, _["bootStrapSamples"] = pmsBootStrapped
		);
}

// .matr2Mc<-function(matrData,laplacian=0) 
S4 _matr2Mc(CharacterMatrix matrData, double laplacian=0) {
//NumericMatrix _matr2Mc(SEXP matrData, double laplacian=0) {
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
  //CharacterMatrix mat(matrData);
  int nRows = matrData.nrow(), nCols = matrData.ncol();
  //Rf_PrintValue(matrData);
  //Rcout << nRows << endl;
  set<string> uniqueVals;
  for(int i = 0; i < nRows; i++) 
  	for(int j = 0; j < nCols; j++) 
		uniqueVals.insert((string)matrData(i, j));	
  //Rcout << uniqueVals << endl;
  //Rf_PrintValue(uniqueVals);
/*
  for(set<string>::iterator it=uniqueVals.begin(); it!=uniqueVals.end(); ++it)
	Rcout << ' ' << *it;
  Rcout << endl;
*/
  int usize = uniqueVals.size();
  NumericMatrix contingencyMatrix (usize, usize);
  contingencyMatrix.attr("dimnames") = List::create(uniqueVals, uniqueVals); 
  
  set<string>::iterator it;
  int stateBegin, stateEnd;
  for(int i = 0; i < nRows; i ++) {
	for(int j = 1; j < nCols; j ++) {
		int k = 0;
  		for(it=uniqueVals.begin(); it!=uniqueVals.end(); ++it, k++) {
			if(*it == (string)matrData(i,j-1)) stateBegin = k;
			if(*it == (string)matrData(i,j)) stateEnd = k;
		}
    //Rcout << stringchar[i] << "->" << stringchar[i + 1] << ": " << posFrom << " " << posTo << endl;
    		contingencyMatrix(stateBegin,stateEnd)++;
	}
  }

  //#add laplacian correction if needed
  //contingencyMatrix=contingencyMatrix+laplacian
  //NumericMatrix transitionMatrix (usize, usize);
  for(int i = 0; i < usize; i ++) {
	double rowSum = 0;
	for(int j = 0; j < usize; j ++) {
    		contingencyMatrix(i,j) += laplacian;
    		rowSum += contingencyMatrix(i,j);
  	}
  	//#get a transition matrix and a DTMC
	for(int j = 0; j < usize; j ++) 
    		//transitionMatrix(i,j) = contingencyMatrix(i,j)/rowSum;
    		contingencyMatrix(i,j) /= rowSum;
  }
  //#get a transition matrix and a DTMC
  //transitionMatrix<-contingencyMatrix/rowSums(contingencyMatrix);
  //outMc<-new("markovchain",transitionMatrix=transitionMatrix);
  
  //Function newClass("new");
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = contingencyMatrix;
  //outMc = (NumericMatrix)newClass("markovchain", contingencyMatrix);
  //outMc =new("markovchain",contingencyMatrix);

  return(outMc);
}


// markovchainFit<-function(data,method="mle", byrow=TRUE,nboot=10,laplacian=0, name, parallel=FALSE)
// [[Rcpp::export]]
List markovchainFit_cpp(SEXP data, String method="mle", bool byrow=true, int nboot=10, double laplacian=0, String name="", bool parallel=false) {
//List markovchainFit_cpp(SEXP data, String method="mle", bool byrow=true, int nboot=10, double laplacian=0, String name="", bool parallel=true) {
  List out;
 // Rf_PrintValue(data);
  //if(class(data) %in% c("data.frame","matrix")) {
  if(Rf_inherits(data, "data.frame") || Rf_inherits(data, "matrix")) { 
	CharacterMatrix mat;
    	//#if data is a data.frame forced to matrix
    	//if(data.attr("class") == "data.frame") data =as.matrix(data);
  	if(Rf_inherits(data, "data.frame")) {
//  		Rcout << "data.frame to matrix" << endl;
		DataFrame df(data);
		//data2 = internal::convert_using_rfunction(df, "as.matrix");
		//Function asMatrix("as.matrix"); // 2x faster than convert_using_frunction
		//mat = asMatrix(df);
		mat = CharacterMatrix(df.nrows(), df.size());
		for(int i = 0; i < df.size(); i++)
			mat(_,i) = CharacterVector(df[i]);
//		Rf_PrintValue(mat);
 	} else {
//  		Rcout << "matrix" << endl;
		mat = data;
	}
    	//byrow assumes distinct observations (trajectiories) are per row
    	//otherwise transpose
  	if(!byrow) mat = _transpose(mat);
   	S4 outMc =_matr2Mc(mat,laplacian);
    	//out<-list(estimate=outMc)
 	out = List::create(_["estimate"] = outMc);
	//out = List(1);
	//out["estimate"] = outMc;
  } else {
    if(method == "mle") out = _mcFitMle(data, byrow);
    if(method == "bootstrap") out = _mcFitBootStrap(data, nboot, byrow, parallel);
    if(method == "laplace") out = _mcFitLaplacianSmooth(data, byrow, laplacian);
  }
//  if(!missing(name)) out$estimate@name<-name
  if(name != "") {
    S4 estimate = out["estimate"];
    estimate.slot("name") = name;
  }
  //((S4)(out["estimate"])).slot("name") = name;
  //if(!name.empty()) ((S4)out["estimate"]).slot("name") = name;
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
sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
#sequence <- data.frame(t(sequence))
#microbenchmark(
  #markovchainFit(data = sequence)#,
  #markovchainFit(data = sequence, method="laplace", laplacian=0.1)#,
  markovchainFit(data = sequence, method="bootstrap")#,
  #markovchainFit(data = sequence, byrow=FALSE)#,
  #markovchainFit_cpp(sequence)
  markovchainFit_cpp(sequence, "bootstrap")
  #markovchainFit_cpp(sequence, "laplace", laplacian=0.1)
  #markovchainFit_cpp(sequence, byrow=FALSE)
#)
*/
/*  markovchainFit_cpp(sequence)
  #markovchainFit_cpp(sequence, byrow=FALSE)
*/
/*
microbenchmark(
  markovchainFit(data = sequence),
  markovchainFit_cpp(100, 10)
)
*/

