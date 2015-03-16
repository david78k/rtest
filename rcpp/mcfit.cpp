#include <Rcpp.h>
#include <omp.h>
#include <unistd.h>

using namespace Rcpp;

template <typename T>
T _transpose(T & m) {      
  int k = m.rows(), n = m.cols();
  T z(n, k);
  z.attr("dimnames") = List::create(colnames(m), rownames(m)); 
  int sz1 = n*k-1;
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

// [[Rcpp::export]]
NumericMatrix createSequenceMatrix_cpp(CharacterVector stringchar, bool toRowProbs=false, bool sanitize=true) {
  CharacterVector elements = unique(stringchar).sort();
  int sizeMatr = elements.size();
  
  NumericMatrix freqMatrix(sizeMatr);
  freqMatrix.attr("dimnames") = List::create(elements, elements); 
  CharacterVector rnames = rownames(freqMatrix);

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
List _mcFitMle(CharacterVector stringchar, bool byrow, double confidencelevel=95.0) {
//List _mcFitMle(CharacterVector stringchar, bool byrow) {
//NumericMatrix _mcFitMle(DataFrame stringchar, bool byrow) {
/*
  initialMatr<-createSequenceMatrix(stringchar=stringchar,toRowProbs=TRUE)
  outMc<-new("markovchain", transitionMatrix=initialMatr,name="MLE Fit")
  if(byrow==FALSE) outMc<-t(outMc)
  out<-list(estimate=outMc)
*/

  // get initialMatr and freqMatr at the same time for speedup
//  NumericMatrix initialMatr = createSequenceMatrix_cpp(stringchar, true);
  CharacterVector elements = unique(stringchar).sort();
  int sizeMatr = elements.size();
  //Rf_PrintValue(elements);
  
  NumericMatrix initialMatr(sizeMatr);
  NumericMatrix freqMatr(sizeMatr);
  initialMatr.attr("dimnames") = List::create(elements, elements); 
  //CharacterVector rnames = rownames(initialMatr);
//  Rf_PrintValue(freqMatrix);
  int posFrom, posTo;
  for(int i = 0; i < stringchar.size() - 1; i ++) {
	for (int j = 0; j < elements.size(); j ++) {
		if(stringchar[i] == elements[j]) posFrom = j;
		if(stringchar[i + 1] == elements[j]) posTo = j;
	}
  	freqMatr(posFrom,posTo)++;
  }
  //initialMatr = freqMatr;
  //Rf_PrintValue(freqMatr);
 
  // sanitize and to row probs
  for (int i = 0; i < sizeMatr; i++) {
  	double rowSum = 0;
  	for (int j = 0; j < sizeMatr; j++) 
  		rowSum += freqMatr(i, j);
  	// toRowProbs
    	for (int j = 0; j < sizeMatr; j++) {
		if(rowSum == 0)
      			initialMatr(i, j) = 1/sizeMatr;
		else
      			initialMatr(i, j) = freqMatr(i, j)/rowSum;
	}
  }
  //Rf_PrintValue(initialMatr);

  if(byrow==false) initialMatr = _transpose(initialMatr);

  // confidence interval
  double criticalValue = 1 - (1 - confidencelevel)/2;
  //double score = 1.644853;
  double zscore = 1.96;
  //double tscore = 2.12; // _tscore(criticalValue, n - 1);
  //double tscore = 1.96;
  if(confidencelevel == 99.9) zscore = 3.3;
  else if(confidencelevel == 99.0) zscore = 2.577;
  else if(confidencelevel == 98.5) zscore = 2.43;
  else if(confidencelevel == 97.5) zscore = 2.243;
  else if(confidencelevel == 90.0) zscore = 1.645;
  else if(confidencelevel == 85.0) zscore = 1.439;
  else if(confidencelevel == 75.0) zscore = 1.151;

  int n = stringchar.size();
  //Rcout << "transition count: " << n << endl;
  //double lowerEndpoint = cellMean - criticalValue*sigma/sqrt(n);
  //double upperEndpoint = cellMean + criticalValue*sigma/sqrt(n);
  NumericMatrix lowerEndpointMatr = NumericMatrix(initialMatr.nrow(), initialMatr.ncol());
  NumericMatrix upperEndpointMatr = NumericMatrix(initialMatr.nrow(), initialMatr.ncol());

  for(int i = 0; i < initialMatr.nrow(); i ++) {
	for(int j = 0; j < initialMatr.ncol(); j ++) {
		lowerEndpointMatr(i,j) = std::max(0.0, std::min(1.0, initialMatr(i, j) - zscore * initialMatr(i,j) / sqrt(freqMatr(i,j))));
		upperEndpointMatr(i,j) = std::max(0.0, std::min(1.0, initialMatr(i, j) + zscore * initialMatr(i,j) / sqrt(freqMatr(i,j))));
		//lowerEndpointMatr(i,j) = std::max(0.0, std::min(1.0, initialMatr(i, j) - score * initialMatr(i,j) / sqrt(n)));
		//upperEndpointMatr(i,j) = std::max(0.0, std::min(1.0, initialMatr(i, j) + score * initialMatr(i,j) / sqrt(n)));
		//lowerEndpoint = tscore(criticalValue, n - 1) * initialMatr(i,j) / sqrt(n);
		//lowerEndpoint = initialMatr(i,j) - criticalValue * sigma / sqrt(n);		
  	}
  }
  lowerEndpointMatr.attr("dimnames") = List::create(elements, elements); 
  upperEndpointMatr.attr("dimnames") = List::create(elements, elements); 
  //Rf_PrintValue(lowerEndpointMatr);
  //Rf_PrintValue(upperEndpointMatr);

  //NumericMatrix outMc(initialMatr); //("markovchain", initialMatr,"MLE Fit");
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = initialMatr;
  outMc.slot("name") = "MLE Fit";  
  //outMc.slot("lowerEndpointMatrix") = lowerEndpointMatr;
  //outMc.slot("upperEndpointMatrix") = upperEndpointMatr;
  
  return List::create(_["estimate"] = outMc
		, _["confidenceInterval"] = List::create(lowerEndpointMatr, upperEndpointMatr)
	);
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
  //Rf_PrintValue(contingencyMatrix);

  List samples, res;
  CharacterVector itemset = rownames(contingencyMatrix);
  int itemsetsize = itemset.size();
  //Rf_PrintValue(itemset);

  Function sample("sample");
  srand(time(NULL));
  for(int i = 0; i < n; i ++) {
	CharacterVector charseq;	
	int rnd = rand()%itemsetsize;
 	String ch = itemset[rnd];
 	//Rcout << std::string(ch) << std::endl;
	//Rcout << rnd << " " << itemset[rnd] << endl;
	//Rf_PrintValue(ch);
//	res = sample(itemset, 1);
//	CharacterVector cv = res[0];
//	ch = cv[0];
//	charseq.push_back(cv[0]);
	charseq.push_back(ch);
	//Rf_PrintValue(charseq);
	for(int j = 1; j < size; j ++) {
		NumericVector probsVector;
		for(int k = 0; k < itemsetsize; k ++) {
			if((std::string)itemset[k] == (std::string) ch) {
				probsVector = contingencyMatrix(k, _);	
				//Rcout << k << " " << (std::string)ch << std::endl;
				//Rf_PrintValue(probsVector);
				break;
			}
		}
		//Rf_PrintValue(probsVector);
  		//srand(time(NULL));
		//rnd = rand()%itemsetsize;
 		//ch = itemset[rnd];
		res = sample(itemset, 1, true, probsVector);
		//Rcout << "sampled: "; 
	//	Rf_PrintValue(res[0]);
		CharacterVector v = res[0];
		//Rf_PrintValue(res[0]);
		//Rf_PrintValue(character);
		ch = v[0];
		charseq.push_back(ch);
 	}
	//samples[[samples.size() + 1]] = charseq;
	samples.push_back(charseq);
	//samples = List::create(clone(samples), charseq);
  }

  //Rf_PrintValue(samples);

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
  //Rcout << "sampleSize: " << sampleSize << endl;
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
		}
//		Rcout << "probEstimated" << endl;
//		Rf_PrintValue(probsEstimated);
		matrMean(i,j) = mean(probsEstimated);
		matrSd(i,j) = sd(probsEstimated);
  	}
  }
  matrMean.attr("dimnames") = List::create(rownames(firstMat), colnames(firstMat)); 
  matrSd.attr("dimnames") = matrMean.attr("dimnames");
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
  //List theList = bootstrapCharacterSequences(data, nboot);
  int n = theList.size();
  //Rcout << "theList.size() " << n << endl;
  List pmsBootStrapped(n);
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
  NumericMatrix transMatr = _toRowProbs(estimateList["estMu"]);
  //Rf_PrintValue(transMatr);
  //int size = temp.nrow();
  
  //NumericMatrix temp = estimateList["estMu"];
  //Function sweep("sweep");
  //NumericMatrix transMatr2 = sweep(temp, 1, _rowSumsC(temp), "/");
  //Rf_PrintValue(transMatr2);
  
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
  std::set<std::string> uniqueVals;
  for(int i = 0; i < nRows; i++) 
  	for(int j = 0; j < nCols; j++) 
		uniqueVals.insert((std::string)matrData(i, j));	
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
  
  std::set<std::string>::iterator it;
  int stateBegin, stateEnd;
  for(int i = 0; i < nRows; i ++) {
	for(int j = 1; j < nCols; j ++) {
		int k = 0;
  		for(it=uniqueVals.begin(); it!=uniqueVals.end(); ++it, k++) {
			if(*it == (std::string)matrData(i,j-1)) stateBegin = k;
			if(*it == (std::string)matrData(i,j)) stateEnd = k;
		}
    //Rcout << stringchar[i] << "->" << stringchar[i + 1] << ": " << posFrom << " " << posTo << endl;
    		contingencyMatrix(stateBegin,stateEnd)++;
	}
  }

  //add laplacian correction if needed
  for(int i = 0; i < usize; i ++) {
	double rowSum = 0;
	for(int j = 0; j < usize; j ++) {
    		contingencyMatrix(i,j) += laplacian;
    		rowSum += contingencyMatrix(i,j);
  	}
  	//get a transition matrix and a DTMC
	for(int j = 0; j < usize; j ++) 
    		//transitionMatrix(i,j) = contingencyMatrix(i,j)/rowSum;
    		contingencyMatrix(i,j) /= rowSum;
  }
  //#get a transition matrix and a DTMC
  //transitionMatrix<-contingencyMatrix/rowSums(contingencyMatrix);
  //outMc<-new("markovchain",transitionMatrix=transitionMatrix);
  
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = contingencyMatrix;

  return(outMc);
}

// [[Rcpp::export]]
List markovchainFit_cpp(SEXP data, String method="mle", bool byrow=true, int nboot=10, double laplacian=0, String name="", bool parallel=false, double confidencelevel=0.95) {
  List out;
  //if(class(data) %in% c("data.frame","matrix")) {
  if(Rf_inherits(data, "data.frame") || Rf_inherits(data, "matrix")) { 
	CharacterMatrix mat;
    	//if data is a data.frame forced to matrix
  	if(Rf_inherits(data, "data.frame")) {
		DataFrame df(data);
		mat = CharacterMatrix(df.nrows(), df.size());
		for(int i = 0; i < df.size(); i++)
			mat(_,i) = CharacterVector(df[i]);
 	} else {
		mat = data;
	}
    	//byrow assumes distinct observations (trajectiories) are per row
    	//otherwise transpose
  	if(!byrow) mat = _transpose(mat);
   	S4 outMc =_matr2Mc(mat,laplacian);
 	out = List::create(_["estimate"] = outMc);
  } else {
    if(method == "mle") out = _mcFitMle(data, byrow);
    if(method == "bootstrap") out = _mcFitBootStrap(data, nboot, byrow, parallel);
    if(method == "laplace") out = _mcFitLaplacianSmooth(data, byrow, laplacian);
  }

  if(name != "") {
    S4 estimate = out["estimate"];
    estimate.slot("name") = name;
    out["estimate"] = estimate;
  }
  return out;
}


/*** R 
library(microbenchmark)
sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
#sequence <- data.frame(t(sequence))
#microbenchmark(
  markovchainFit(data = sequence)
  #markovchainFit(data = sequence, method="laplace", laplacian=0.1),
  #markovchainFit(data = sequence, method="bootstrap")#,
  #mcfit(data = sequence, method="bootstrap"),
  #markovchainFit(data = sequence, byrow=FALSE)#,

  markovchainFit_cpp(sequence)
  #markovchainFit_cpp(sequence, "bootstrap")
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

