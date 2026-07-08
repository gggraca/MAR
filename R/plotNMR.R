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
#' @param ppm If set to \code{TRUE} (default option) the ppm values will be  
#' added to the plot x scale. The variable indexes will used if
#' set to \code{FALSE}. The chemical shift values can also be added as a
#' numeric vector of the same size as the number of rows or columns of
#' \code{X}.
#' @param title Title of the plot (optional).
#'
#' @return A stacked plot of the full spectra or spectra region.
#' @examples 
#' data("urine")
#' plotNMR(urine)
#' 
#' @export
plotNMR <- function(X, f1=max(X[,1]), f2=min(X[,1]), ppm=TRUE, title=""){
	 #if the matrix is not in the format [variables,samples] it will be transposed
    if(dim(X)[1] < dim(X)[2]){
        X<-t(X)
    }
	if(length(ppm) > 1){
		ppm <- as.numeric(ppm)
        lim <- which(ppm > f2 & ppm < f1)
        l2 <-lim[1]
        l1 <-lim[length(lim)]
        matplot(ppm[l1:l2], X[l1:l2,], type = "l", 
            xlim=rev(range(ppm[l1:l2])),
            xlab="chemical shift (ppm)", 
            ylab="intensity (a.u.)", main=title, lty=1)
	} else {
        lim <- which(X[,1]>f2 & X[,1]<f1)
        l2 <- lim[1]
        l1 <- lim[length(lim)]
        matplot(X[l1:l2,1], X[l1:l2,2:dim(X)[2]], type = "l", 
                xlim=rev(range(X[l1:l2,1])),
                xlab="chemical shift (ppm)", 
                ylab="intensity (a.u.)", main=title, lty=1)
        if(isFALSE(ppm)){
            xscale <- lim
            matplot(xscale, X[l1:l2,2:dim(X)[2]], type="l",
                xlab="variable index", 
                ylab="intensity (a.u.)", 
                main=title, lty=1)
        }
	}
}
