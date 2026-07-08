#' @title 1D SHY using 1H-NMR spectra and other Y data 
#' 
#' @description
#' Function to run 1D Statistical Hetero Spectroscopy (SHY) using
#' using 1H-NMR spectra as X matrix and other matrix-like Y dataset,
#' for instance an XCMS output for an LC-MS run or a low-resolution
#' spectroscopy data (e.g. infrared or Raman spectra).
#' 
#' @author Goncalo Graca, Imperial College London 
#'
#' @param X A matrix of spectra containing ppm values in the first row
#' and intensity values for each sample in the remaining rows If a 
#' matrix with ppm x samples is provided instead, it will be transposed
#' prior to calculations and plotting.
#' @param Y  A matrix of samples as rows and columns as variables, e.g.
#' a processed LC-MS output or or a low-resolution spectroscopy data 
#' matrix (e.g. infrared, Raman or UV-VIS spectra).
#' 
#' @param YvarName Name of the variable to correlate with the NMR spectra as it
#' is represented in the Y matrix.   
#' 
#' @param expansion A vector containing the end and start of a spectral region
#' in ppm. The default value is \code{NULL}, which will plot the full chemical
#' shift range.
#'
#' @return A a plot of the chemical shift vs covariance between the selected Y
#' variable and the NMR spectra. The correlation between the Y selected variable
#' and the NMR matrix will be color-coded on the covariance line.
#' 
#' @examples
#' data("mtbls1")
#' X <- mtbls1[-1,]
#' Y <- mtbls1[-1,]
#' colnames(Y) <- mtbls1[1,]
#' shy1d(X, Y, YvarName="3.04955836737954")
#' 
#' @importFrom ggplot2 ggplot geom_line scale_x_reverse theme_light theme xlab 
#' @importFrom ggplot2 scale_color_gradientn
#' 
#' @export
shy1d <- function(X, Y, YvarName="", ppm=NULL, expansion=NULL, title="") {
	if(nrow(X) > ncol(X)){
        X <- t(X)
    }
	if(is.null(ppm)){
        ppm <- as.numeric(colnames(X))
        idx <- which(colnames(Y) == YvarName)
        c <- cor(X, Y[,idx])
        cv <- cov(X, Y[,idx])
	} else {
        ppm <- X[1,]
        idx <- which(colnames(Y) == YvarName)
        c <- cor(X[-1, ], Y[,idx])
        cv <- cov(X[-1, ], Y[,idx])
	}
    s <- data.frame(ppm=ppm, covariance=cv, correlation=c)
    jet.colors <- colorRampPalette(c("#00007F", 
                                     "blue", 
                                     "#007FFF", 
                                     "cyan", 
                                     "#7FFF7F", 
                                     "yellow", 
                                     "#FF7F00", 
                                     "red", 
                                     "#7F0000"))
    if (!is.null(expansion)) {
        ppmDiv <- (expansion[1]-expansion[2])/10 # unit for ppm scale expansion 
        p <- ggplot(s, aes(x=ppm, y=covariance, color=correlation)) + 
                geom_line() + 
                scale_x_reverse(breaks=seq(expansion[2], expansion[1], ppmDiv), 
                            limits=expansion) +
                theme_light() + 
                theme(panel.grid.major=element_blank(), 
                      panel.grid.minor=element_blank(), 
                      panel.background=element_blank(), 
                      axis.line=element_line(colour="black")) +
                xlab("chemical shift (ppm)") + 
                scale_color_gradientn(colours = jet.colors(256)) +
                ggtitle(title)
    } else {
        p <- ggplot(s, aes(x=ppm, y=covariance, color=correlation)) + 
                geom_line() + 
                scale_x_reverse() + 
                theme_light() + 
                theme(panel.grid.major=element_blank(), 
                      panel.grid.minor=element_blank(), 
                      panel.background=element_blank(), 
                      axis.line = element_line(colour="black")) +
                xlab("chemical shift (ppm)") + 
                scale_color_gradientn(colours=jet.colors(256)) +
    		    ggtitle(title)
    }
    suppressWarnings(print(p))
}

