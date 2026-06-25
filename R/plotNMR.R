#' @title Plot stacked 1H-NMR spectra
#' 
#' @description
#' R function for plotting  1H-NMR spectra.
#' 
#' 
#' @author Goncalo Graca, Imperial College London 
#'
#' @param X A matrix of spectra containing ppm values in the first column
#' and intensity values for each sample in the remaining columns. If a 
#' matrix with samples x ppm is provided instead, it will be transposed
#' prior to plotting.
#' @param f1 Chemical shift of the leftmost ppm value to plot.
#' @param f2 Chemical shift of the rightmost ppm value to plot. 
#' @param ppm Logical value. If \code{TRUE} the ppm values will be  
#' as the plot x scale (default option). The variable indexes will used if
#' set to \code{FALSE}.
#' @param title Title of the plot (optional).
#'
#' @return A stacked plot of the full spectra or spectra region.
#' @example 
#' @export
plotNMR<-function(X, f1 = max(X[,1]), f2 = min(X[,1]), ppm = TRUE, title = ""){
	 #if the matrix is not in the format [variables,samples] it will be transposed
    if(dim(X)[1] < dim(X)[2]){
        X<-t(X)
    }
    lim <-which(X[,1]>f2 & X[,1]<f1)
    l2 <-lim[1]
    l1 <-lim[length(lim)]
    matplot(X[l1:l2,1], X[l1:l2,2:dim(X)[2]], type = "l", 
            xlim=rev(range(X[l1:l2,1])),
            xlab="chemical shift (ppm)", 
            ylab="intensity (a.u.)", main=title, lty=1)
    if(!ppm){
        xscale <- lim
        matplot(xscale, X[l1:l2,2:dim(X)[2]], type="l",
            xlab="variable index", 
            ylab="intensity (a.u.)", 
            main=title, lty=1)
  }
}