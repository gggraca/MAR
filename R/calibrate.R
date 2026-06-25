#' @title Calibrate  1H-NMR spectra
#' 
#' @description
#' Function to align spectra according to a particular signal.
#  To use after 'read1DNMR' function.
#' 
#' 
#' @author Goncalo Graca, Imperial College London 
#'
#' @param X A matrix of spectra containing ppm values in the first column
#' and intensity values for each sample in the remaining columns. If a 
#' matrix with samples x ppm is provided instead, it will be transposed
#' prior to plotting.
#' @param f1 Chemical shift of the leftmost ppm value to plot containing 
#' the spectral reference (e.g. ERETIC, TSP, DSS, glucose...).
#' @param f2 Chemical shift of the rightmost ppm value to plot containing 
#' the spectral reference (e.g. ERETIC, TSP, DSS, glucose...). 
#'
#' @return A matrix of 1H-NMR spectra with the calibrated ppm scale.
#' 
#' @examples
#' data("urine")
#' cs <- calibrate(urine, 0.1, -0.1)
#' 
#' @export
calibrate <- function(X,f1,f2){
    idx <- which(X[,1] > f2 & X[,1] < f1)
    ind <- apply(X[idx,2:ncol(X)], 2, which.max)
    med <- ind[1]
    d <- as.numeric(lapply(ind,function(X) X-med))
    n <- nrow(X)
    p <- ncol(X)
    M <- matrix(rep(NA,times=n*p) , nrow = n)
    colnames(M) <- colnames(X)
    M[,1] <- X[,1]
    for(i in 1:length(d)){
        if (d[i] == 0){
        M[,i+1]<-X[,i+1]
        next
    }
    if(d[i] < 0){
        M[,i+1] <- c(rep(0,times=abs(d[i])), X[1:(dim(X)[1]-abs(d[i])),i+1])
      next
    }
    if(d[i] > 0){
        M[,i+1]<-c(X[(d[i]+1):(dim(X)[1]),i+1],rep(0,times=abs(d[i])))
        next
    }
    }
    return(M)
 }
