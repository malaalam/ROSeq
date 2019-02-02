##' Sample data
##'
##' Tung et al. 2017 dataset. Read count data corresponding to 288 single cells each was available from the three human induced pluripotent stem cell lines (NA19098, NA19101 and NA19239). Bulk RNA-seq data corresponding to these classes was present in the
##' form of three replicates each, for both the datasets. Here, the data corresponding to individuals (NA19101 and NA19239) is made available.
##'
##' @name Tung2017
##' @usage data("Tung2017")
##' @docType data
##' @references \url{https://www.nature.com/articles/srep39921}
##' @export
NULL

##' Sample data
##'
##' Trapnell 2014
##'
##' @name Trapnell2014
##' @usage data("Trapnell2014")
##' @docType data
##' @references \url{https://www.nature.com/articles/srep39921}
##' @export
NULL

##' Sample data
##'
##' Trapnell 2014 (curated)
##'
##' @name Trapnell2014_T0_T24
##' @usage data("Trapnell2014_T0_T24")
##' @docType data
##' @references \url{https://www.nature.com/articles/srep39921}
##' @export
NULL

##' @title Computes differential analysis for a filtered and normalized read count matrix
##'
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param scdata The normalised and filtered, read count matrix
##' @param scgroups The location of the two sub-populations
##' @param mc.cores The number of cores to be used
##' @return pValuesAdjusted A vector containing FDR adjusted p significance values

callROSeq<-function(scdatafn, scgroups, numCores, cOne, cTwo){

  geneIndex<-1:nrow(scdatafn)
  library("pbmcapply", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
  results <- pbmclapply(geneIndex, initiateAnalysis, scdata=scdatafn, scgroups=scgroups, classOne=cOne, classTwo=cTwo, mc.cores=numCores)
  return(results)
}



##' @title Computes differential analysis for a given gene
##'
##' @description Takes in the row index which corresponds to a gene and evaluates for differential expression across two cell types.
##' @param gene The row index of the normalised and filtered, read count matrix
##' @param scdata The normalised and filtered, read count matrix
##' @param scgroups The location of the two sub-populations
##' @return combinedResults A vector containing 12 values (gr1: a, g1: b, gr1: A, gr1: number of bins, gr1: R2, gr2: a, gr2: b, gr2: A, gr2: number of bins, gr2: R2, T, p)

initiateAnalysis<-function(gene, scdata, scgroups, classOne, classTwo){

  sp<-scdata[gene, ]
  spOne<-scdata[gene, which(scgroups==classOne)]
  spTwo<-scdata[gene, which(scgroups==classTwo)]
  geneStats<-getDataStatistics(sp, spOne, spTwo)


  results_groupOne <-findParams(spOne, geneStats)
  results_groupTwo <-findParams(spTwo, geneStats)
  T<-tryCatch({computeDEG(results_groupOne, results_groupTwo)}, warning=function(w) {NA}, error=function(esp) {NA})
  pValues<-pchisq(T, df=2, lower.tail=FALSE)
  combinedResults<-c(results_groupOne[1], results_groupOne[2], results_groupOne[3], results_groupOne[4], results_groupOne[5], results_groupTwo[1], results_groupTwo[2], results_groupTwo[3], results_groupTwo[4], results_groupTwo[5], T, pValues)
  return(combinedResults)

}
##' @title Evaluates statistics of the read counts corresponding to the gene
##'
##' @description Takes in the complete read count vector corresponding to the gene (sp) and also the data corresponding to the two sub-populations (sp1 and sp2)
##'
##' @param sp The complete (normalized and filtered) read count data corresponding to the gene in question
##' @param spOne The (normalized and filtered) read count data corresponding to the first sub-population
##' @param spTwo The (normalized and filtered) read count data corresponding to the second sub-population
##' @return geneStats A vector containing 7 values corresponding to the gene data (maximum, minimum, mean, standard deviation, upper multiple of standard deviation, lower multiple of standard deviation and log2(fold change))

getDataStatistics<-function(sp, spOne, spTwo){
  maxds<-max(sp)
  minds<-min(sp)
  meands<-mean(sp)
  stdds<-sd(sp)
  ceilds<-ceiling((maxds-meands)/stdds)
  floords<-floor((minds-meands)/stdds)
  log2FC<-abs(log2(mean(spOne)/mean(spTwo)))
  geneStats<-c(maxds, minds, meands, stdds, ceilds, floords, log2FC)
  return (geneStats)
}
##' @title Finds the optimal values of parameters a and b that model the probability distribution of ranks, by Maximising the Log-Likelihood
##'
##' @description Takes in as input the read count data corresponding to one sub-population and the typical gene statistics. Then it splits the entire range into equally sized bins of size \eqn{k * \sigma}, where k is a scalar with a default value of 0.05, and \eqn{\sigma} is the standard deviation of the pulled expression estimates across the cell-groups. Each of these bins corresponds to a rank. Therefore, for each group, cell frequency for each bin maps to a rank.  These frequencies are normalized group-wise by dividing by the total cell count within a concerned group.
##' @param ds The (normalized and filtered) read count data corresponding to a sub-population
##' @param geneStats A vector containing 7 values corresponding to the gene data (maximum, minimum, mean, standard deviation, upper multiple of standard deviation, lower multiple of standard deviation and log_{2}(fold change))
##' @return results A vector containing 5 values (a, b, A, number of bins, R2)


findParams<-function(ds, geneStats){
  #a<-25
  #b<-25
  step<-0.05
  #while(a > 20 || b > 20){


    meands<-geneStats[3]
    stdds<-geneStats[4]
    ceilds<-geneStats[5]
    floords<-geneStats[6]

    binNumber<-length(seq(floords, ceilds-step, step))
    rs<-c(rep(NA), binNumber)
    count<-1
    for (i in seq(floords, ceilds-step, step)){
      LL<- meands+i*stdds
      UL <-meands+(i+step)*stdds
      rs[count]<-length(intersect(which(ds<UL), which(ds>=LL)))
      count<-count+1
    }

    fds<-rs[!is.na(rs)]
    number_of_bins<-length(fds)
    rank<-1:number_of_bins
    read_count_sorted<-sort(fds, decreasing=TRUE)
    normalized_read_count_sorted<-read_count_sorted/sum(read_count_sorted)
    model<-optim(par = c(0.25, 3), minimizeNLL, r=rank, readCount=normalized_read_count_sorted, method = "Nelder-Mead")
    #model<-constrOptim(par = c(0.25, 3), minimizeNLL, r=rank, readCount=normalized_read_count_sorted, method = "Nelder-Mead", ui=rbind(c(0,-1)), ci=c(-20))
    if (is.na(model$value)){
      a<-NA
      b<-NA
    }else{
      a<-model$par[1]
      b<-model$par[2]
    }
    #step=step+0.01

  #}
  #a<-min(a, 20)
  #b<-min(b, 20)

  A<-1/sum((number_of_bins+1-rank)^b/(rank^a))
  f<-A*((number_of_bins+1-rank)^b)/(rank^a)
  SS_res<-sum((normalized_read_count_sorted-f)^2)
  SS_tot<-sum((normalized_read_count_sorted-mean(normalized_read_count_sorted))^2)
  R2<-1-SS_res/SS_tot
  results<-c(a, b, A, number_of_bins, R2)
  return(results)
}

##' @title Minimizes the Negative Log-Likelihood by iterating across values of parameters a and b
##'
##' @description Takes in as input a vector of values (coefficients), the number of bins and the normalized probability dsitribution of ranks
##' @param coefficients A vector containing two values for a and b
##' @param r The number of bins
##' @param readCount A vector of (normalized) frequency of read counts that fall within each bin
##' @return NLL Negative-Log Likelihood for the input coefficients
##' @seealso \code{\link{findParams}}

minimizeNLL<-function(coefficients, r, readCount){

  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  sumReadCount<-sum(readCount)
  A <-1/sum(((N+1-r)^b)/(r^a))
  NLL<-a*sum(readCount*log(r)) - b*sum(readCount*log(N+1-r)) - sumReadCount*log(A)
  return (NLL)
}

##' @title Computes differential expression for the gene in question, by comparing the optimal parameters for sub-populations one and two
##' @description  Uses the (asymptotically) optimum two-sample Wald test  based on the MLE of the parameters and its asymptotic variances given by the inverse of the Fisher information matrix
##' @param results_1 A vector corresponding to sub-population one and containing 5 values (a, b, A, number of bins, R2)
##' @param results_2 A vector corresponding to sub-population two and containing 5 values (a, b, A, number of bins, R2)
##' @return T  The Wald test statistic for testing the null hypothesis
##' @seealso \code{\link{getI}}, \code{\link{findParams}}



computeDEG<-function(results_1, results_2){

  I_1<-getI(results_1)
  I_2<-getI(results_2)

  I1<-matrix(c(I_1[1], I_1[2], I_1[3], I_1[4]), nrow = 2, ncol=2)
  I2<-matrix(c(I_2[1], I_2[2], I_2[3], I_2[4]), nrow = 2, ncol=2)

  V1<-solve(I1, tol = 1e-20)
  V2<-solve(I2, tol = 1e-20)
  m<-results_1[4]
  n<-results_2[4]
  a1<-results_1[1]
  b1<-results_1[2]
  a2<-results_2[1]
  b2<-results_2[2]
  #if(a1==20 || a2==20){
  #  T<-1000
  #}
  #else{
  w<-n/(m + n)
  T<-(m*n/(m+n))* t(matrix(c(a1-a2,b1-b2), nrow=2, ncol=1)) %*% solve(w *V1 + (1-w)* V2, tol = 1e-20) %*% matrix(c(a1-a2,b1-b2), nrow=2, ncol=1)
  #}
  return(T)
}

##' @title Computes the Fisher Information Matrix
##' @description The Fisher Information Matrix and its derivatives are essential to evulate the minima of log likelihood
##' @param results A vector containing 5 values (a, b, A, number of bins, R2)
##' @return I  The Fisher Information Matrix

getI<-function(results){

  rank<-1:results[4]
  coefficients<-c(results[1], results[2])


  u1<-getu1(coefficients, rank)
  v<-getv(coefficients, rank)
  u2<-getu2(coefficients, rank)


  du1da<-getdu1da(coefficients, rank)
  du1db<-getdu1db(coefficients, rank)
  du2da<-getdu2da(coefficients, rank)
  du2db<-getdu2db(coefficients, rank)
  dvda<-getdvda(coefficients, rank)
  dvdb<-getdvdb(coefficients, rank)

  d2logAda2<-getd2logAda2( u1, v, du1da, dvda)
  d2logAdb2<-getd2logAdb2( u2, v, du2db, dvdb)
  d2logAdbda<-getd2logAdbda( u1, v, du1db, dvdb)
  d2logAdadb<-getd2logAdadb( u2, v, du2da, dvda)

  I<-c(-d2logAda2, -d2logAdadb, -d2logAdbda, -d2logAdb2)
  return(I)
}


##' @title Computes u1
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return u1

getu1<-function(coefficients, r){

  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)

  num1<-(N+1-r)^b
  num2<-log(r)
  den1<-r^a
  u1<-sum(num1*num2/den1)
  return(u1)
}

##' @title Computes v
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return v

getv<-function( coefficients, r){

  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  den1<-r^a
  v<-sum(num1/den1)
  return(v)
}

##' @title Computes u2
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return u2

getu2<-function(coefficients, r){

  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-log(N+1-r)
  den1<-r^a
  u2<-(-sum(num1*num2/den1))
  return(u2)
}




##' @title Finds the first derivative of u1 with respect to a. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return du1da

getdu1da<-function(coefficients, r){

  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-(log(r))^2
  den1<-r^a
  du1da<- -sum(num1*num2/den1)
  return (du1da)
}

##' @title Finds the first derivative of u1 with respect to b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return du1db

getdu1db<-function(coefficients, r){
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-log(r)
  num3<-log(N+1-r)
  den1<-r^a
  du1db<-sum(num1*num2*num3/den1)
  return(du1db)
}

##' @title Finds the first derivative of u2 with respect to a. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return du2da

getdu2da<-function(coefficients, r){

  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)

  num1<-(N+1-r)^b
  num2<-log(r)
  num3<-log(N+1-r)
  den1<-r^a
  du2da<-sum(num1*num2*num3/den1)
  return(du2da)
}

##' @title Finds the first derivative of u2 with respect to b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return du2db

getdu2db<-function( coefficients, r){



  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)

  num1<-(N+1-r)^b
  num2<-(log(N+1-r))^2
  den1<-r^a
  du2db<-(-sum(num1*num2/den1))
  return(du2db)
}

##' @title Finds the first derivative of v with respect to a. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return dvda

getdvda<-function(coefficients, r){

  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)

  num1<-(N+1-r)^b
  num2<-log(r)
  den1<-r^a
  dvda<-(-sum(num1*num2/den1))
  return(dvda)
}

##' @title Finds the first derivative of v with respect to b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return dvdb

getdvdb<-function( coefficients, r){


  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)

  num1<-(N+1-r)^b
  num2<-log(N+1-r)
  den1<-r^a
  dvdb<-sum(num1*num2/den1)
  return(dvdb)
}



##' @title Finds the double derivative of A with with respect to a. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param u1 u1
##' @param v v
##' @param du1da First derivative of u1 with respect to parameter a
##' @param dvda First derivative of v with respect to parameter a
##' @return d2logAda2
getd2logAda2<-function(u1, v, du1da, dvda){

  num1<-v*du1da
  num2<-u1*dvda
  den1<-v^2
  d2logAda2<-(num1-num2)/den1
  return(d2logAda2)
}

##' @title Finds the double derivative of A with with respect to a and b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param u2 u2
##' @param v v
##' @param du2da First derivative of u2 with respect to a
##' @param dvda First derivative of v with respect to a
##' @return d2logAdadb

getd2logAdadb<-function( u2, v, du2da, dvda){

  num1<-v*du2da
  num2<-u2*dvda
  den1<-v^2
  d2logAdadb<-(num1-num2)/den1
  return(d2logAdadb)
}

##' @title Finds the double derivative of A with with respect to b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param u2 u2
##' @param v v
##' @param du2db First derivative of u2 with respect to b
##' @param dvdb First derivative of v with respect to
##' @return d2logAdb2

getd2logAdb2<-function( u2, v, du2db, dvdb){


  num1<-v*du2db
  num2<-u2*dvdb
  den1<-v^2
  d2logAdb2<-(num1-num2)/den1
  return(d2logAdb2)
}

##' @title Finds the double derivative of A with with respect to a and b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param u1 u1
##' @param v v
##' @param du1db First derivative of u1 with respect to b
##' @param dvdb First derivative of v with respect to b
##' @return d2logAdbda

getd2logAdbda<-function( u1, v, du1db, dvdb){

  num1<-v*du1db
  num2<-u1*dvdb
  den1<-v^2
  d2logAdbda<-(num1-num2)/den1
  return(d2logAdbda)
}

##' @title Generates a ROC curve for prediction by Wilcoxon Rank-Sum method
##'
##' @description Takes in pValues generated by Wilcoxon Method and produces an AUC plot
##' @param pVals
##' @return AUC :  area under the ROC curve

DE_Quality_AUC_plot<- function(pVals, newcolor) {
  pVals <- pVals[names(pVals) %in% GroundTruth$DE |
                   names(pVals) %in% GroundTruth$notDE]
  truth <- rep(1, times = length(pVals));
  truth[names(pVals) %in% GroundTruth$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 4, 0, 0))
  ROCR::plot(perf,xaxis.cex.axis=1.5, xaxis.col='black', xaxis.col.axis="black", yaxis.col='black',
             yaxis.cex.axis=1.5, yaxis.las=1.5, xaxis.lwd=1.5, yaxis.lwd=1.5, yaxis.col.axis="black", cex.lab=1.5, cex.main=1.5,
             lwd=10,col=newcolor, cex=1.5)
  aucObj <- ROCR::performance(pred, "auc")
  return(aucObj@y.values[[1]])
}

##' @title Generates a ROC curve for prediction by ROSeq method
##'
##' @description Takes in pValues generated by ROSeq Method and produces an AUC plot
##' @param pVals
##' @return AUC : area under the ROC curve
DE_Quality_AUC_plot_add <- function(pVals, newcolor) {
  pVals <- pVals[names(pVals) %in% GroundTruth$DE |
                   names(pVals) %in% GroundTruth$notDE]
  truth <- rep(1, times = length(pVals));
  truth[names(pVals) %in% GroundTruth$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  ROCR::plot(perf, add=TRUE, lwd=10, col=newcolor, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  aucObj <- ROCR::performance(pred, "auc")
  return(aucObj@y.values[[1]])
}


##' @title Runs MAST Method
##'
##' @description Takes in filtered data and the factors and returns pvalues
##' @param FilteredData filtered data
##' @param classOne class one
##' @param classTwo class two
##' @param scgroups
##' @return pVals

Run_MAST <- function(FilteredData, classOne, classTwo, scgroups){
  library("tibble", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
  library("MAST", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
  log_counts <- log(FilteredData+1)/log(2)
  fData = data.frame(names=rownames(log_counts))
  rownames(fData) = rownames(log_counts)
  cData = data.frame(cond=scgroups)
  rownames(cData) = colnames(log_counts)

  obj <- FromMatrix(as.matrix(log_counts), cData, fData)
  colData(obj)$cngeneson <- scale(colSums(assay(obj)>0))

  cond <- factor(colData(obj)$cond)
  zlmCond <- zlm(~cond + cngeneson, obj)
  summaryCond <- summary(zlmCond, doLRT=paste0('cond', classTwo))
  summaryDt <- summaryCond$datatable
  summaryDt <- as.data.frame(summaryDt)
  pVals <- unlist(summaryDt[summaryDt$component == "H",4]) # H = hurdle model
  names(pVals) <- unlist(summaryDt[summaryDt$component == "H",1])
  return(pVals)
}

##' @title Runs ROTS Method
##'
##' @description Takes in filteredData and the factors and returns pvalues
##' @param FilteredData
##' @param classOne
##' @param classTwo
##' @param scgroups
##' @param B
##' @param K
##' @param FDR
##' @return pVals
Run_ROTS <- function(FilteredData, classOne, classTwo, scgroups,  B, K , FDR){
  library("ROTS", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
  library("edgeR", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
  require(ROTS)

  # FilteredData is the input raw count table
  # B is the number of bootstraps
  # k is the number of top ranked genes
  # FDR is the fdr threshold for the final detections
  # TMM normalization over the raw count matrix:

  TMM_Filtered_Data = TMMnormalization(FilteredData)

  # Running ROTS over the filtered and normalized data:
  groups2<-ifelse(scgroups==classOne, 0, 1)
  ROTS_filtered_TMM = ROTS(data = TMM_Filtered_Data, groups2, B = 1000, K = 6000 )

  return(ROTS_filtered_TMM)

}

##' @title Normalize filtered data using TMM method
##'
##' @description Takes in filteredData and returns Normalized Data
##' @param countTable
##' @return Filtered Data
TMMnormalization <- function(countTable){
  ## TMM normalization based on edgeR package:
  require("edgeR")

  nf=calcNormFactors(countTable ,method= "TMM")
  nf= colSums(countTable)*nf
  scalingFactors = nf/mean(nf)
  countTableTMM <- t(t(countTable)/scalingFactors)

  return(countTableTMM)

}
##' @title Run SCDE
##'
##' @description Takes in filteredData and returns Normalized Data
##' @param countTable
##' @return Filtered Data
Run_SCDE <- function(scdatafiltered, classOne, classTwo, scgroups, samplesize){
  # Solving Stage - Testing results with SCDE. sampleSize = 96 cells extracted randomly from each sub population.
  # Normalization = Inbuilt within SCDE

  library("scde", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
  c1<-sample(which(scgroups==classOne), samplesize)
  c2<-sample(which(scgroups==classTwo), samplesize)
  sc1<-scdatafiltered[,  c1]
  sc2<-scdatafiltered[,  c2]
  sc<-cbind.data.frame(sc1, sc2)
  #sc.labels <- substr(colnames(sc),0,7)
  #sc_groups <- factor(sc.labels,levels=c(classOne,classTwo))
  sc_groups <- factor(gsub(paste0("(", classOne, '|', classTwo, ").*"), "\\1", colnames(sc)), levels = c(classOne, classTwo))
  sc_groups<-table(sc_groups)
  counts<-sc
  group<-sc_groups
  cnts <- apply(
    counts,
    2,
    function(x) {
      storage.mode(x) <- 'integer'
      return(x)
    }
  )
  names(group) <- 1:length(group)
  colnames(cnts) <- 1:length(group)
  o.ifm <- scde::scde.error.models(
    counts = cnts,
    groups = group,
    n.cores = 1,
    threshold.segmentation = TRUE,
    save.crossfit.plots = FALSE,
    save.model.plots = FALSE,
    verbose = 1,
    min.size.entries = 2
  )
  priors <- scde::scde.expression.prior(
    models = o.ifm,
    counts = cnts,
    length.out = 400,
    show.plot = FALSE
  )
  resSCDE <- scde::scde.expression.difference(
    o.ifm,
    cnts,
    priors,
    groups = group,
    n.randomizations = 100,
    n.cores = 1,
    verbose = 1
  )
  # Convert Z-scores into 2-tailed p-values
  SCDE_pVals <- pnorm(abs(resSCDE$cZ), lower.tail = FALSE) * 2
  names(SCDE_pVals)<-rownames(single_counts)
  return(SCDE_pVals)
}

##' @title Add Legend
##'
##' @description Takes in auc of mast, roseq, rots, scde, wilcoxon
##' @param aucmast
##' @param aucroseq
##' @param aucrots
##' @param aucscde
##' @param aucwilcoxon
addLegend<-function(aucmast, aucroseq, aucrots, aucscde, aucwilcoxon){
  legend("bottomright",c(paste0('MAST (AUC = ', aucmast, ')'), paste0('ROSeq (AUC = ', aucroseq, ')'), paste0('ROTS (AUC = ', aucrots, ')'), paste0('SCDE (AUC = ', aucscde, ')'), paste0('Wilcoxon (AUC = ', aucwilcoxon, ')')),col=c('blue', 'black',  'magenta', 'green', 'red'),
         lwd=10, cex=1)

}
