#' @title 1D STOCSY for 1H-NMR spectra
#' 
#' @description
#' Function to run 1D Statistical total correlation spectroscopy (STOCSY).
#' 
#' @author Goncalo Graca, Imperial College London 
#'
#' @param spectra A matrix of spectra containing ppm values in the first row
#' and intensity values for each sample in the remaining rows If a 
#' matrix with ppm x samples is provided instead, it will be transposed
#' prior to calculations and plotting.
#' @param driver Chemical shift of the peak to which all other intensities
#' will be correlated. Covariance is also calculated between the driver peak
#' and all other points in the spectrum.
#' @param expansion A vector containing the end and start of a spectral region
#' in ppm. The default value is \code{NULL}, which will plot the full chemical
#' shift range.
#'
#' @return A a plot of 1d STOCSY, i.e. chemical shift vs covariance, 
#' color-coded by the correlation coefficient.
#' 
#' @examples
#' data("mtbls1")
#' stocsy1d(mtbls1, 3.04)
#' 
#' @export
stocsy1d <- function(spectra, driver, expansion=NULL) {
	if(nrow(spectra) > ncol(spectra)){
        spectra <- t(spectra)
    }
    library(ggplot2)
    ppm <- as.numeric(spectra[1,])
    idx <- which.min(abs(driver - ppm))
    c <- cor(spectra[-1, ], spectra[-1, idx])
    cv <- cov(spectra[-1, ], spectra[-1, idx])
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
        ggplot(s, aes(x=ppm, y=covariance, color=correlation)) + 
            geom_line() + 
            scale_x_reverse(breaks=seq(0, 10, 0.5), limits=expansion) +
            theme_light() + 
            theme(panel.grid.major=element_blank(), 
                  panel.grid.minor=element_blank(), 
                  panel.background=element_blank(), 
                  axis.line=element_line(colour="black")) +
            xlab("chemical shift (ppm)") + 
            scale_color_gradientn(colours = jet.colors(256))
    } else {
        ggplot(s, aes(x=ppm, y=covariance, color=correlation)) + 
            geom_line() + scale_x_reverse(breaks=seq(0, 10, 0.5)) + 
            theme_light() + 
            theme(panel.grid.major=element_blank(), 
                  panel.grid.minor=element_blank(), 
                  panel.background=element_blank(), 
                  axis.line = element_line(colour="black")) +
            xlab("chemical shift (ppm)") + 
            scale_color_gradientn(colours=jet.colors(256))
    }
}

