#' @title Save the results of the quantification of 1H-NMR spectra
#' 
#' @description
#' R function for plotting and exporting  1H-NMR spectra quantification results.
#' Plots are generated for each integrated peak and the quantification result
#' of each sample or groups of samples (if a group label is provided).
#' Plots are saved as \code{.jpeg} and the quantification table as \code{.csv}.
#' 
#' @author Goncalo Graca, Imperial College London 
#'
#' @param resultObject A list containing containing several objects: 
#' the spectra matrix \code{spectra}, data frame containing the integration 
#' regions, the number of protons and chemical shift of each metabolite signal,
#' \code{integrationRegions}, the integration and quantification results 
#' (\code{integrals} and \code{quantification}) and other the metadata such as 
#' \code{baseline}, reference name(\code{reference}) and concentration 
#' (\code{refConc}), the groups labels \code{groups}, sample names 
#' \code{spNames}, and project name \code{projectName}.
#' @param DirPath Path to the user-defined results folder.
#' 
#' @return An overlay of the integrated spectral region and boxplots are saved if 
#' sample group information is provided in the \code{groups} parameter or
#' otherwise barplots of each sample concentration are saved.  
#' Images and quantification results table are saved to the user-defined folder 
#' \code{DirPath}.
#' 
#' @examples
#' userDir <- tempdir()
#' saveQuantification(result, DirPath=userDir)
#' 
#' @export
saveQuantification <- function(resultObject, DirPath=""){
	# check if quantification data exists on resultObject
	quan <- resultObject$quantification
	if(!exists("quan")) {
		stop("No quantification result found in the results!")
	}
	# get spectra from resultObject
    M <- resultObject$spectra
    # if the matrix is not in the format [variables,samples] it will be transposed
    if(nrow(M) < ncol(M)){ 
        M <- t(M)
    }
    # get integrated regions from resultObject
    reg <- resultObject$integrationRegions
    
    for(i in seq_len(nrow(reg))){
        ints <- which(M[,1] > reg[i,"ppm.end"] & M[,1] < reg[i,"ppm.start"])
        a <- ints[1]
        b <- ints[length(ints)]
        if(isTRUE(resultObject$baseline)){
            T <- bas(M, reg[i,"ppm.start"], reg[i,"ppm.end"])
        }
        
        # get groups and sample names from resultObject
        groups <- resultObject$groups
        SpNames <- resultObject$SpNames
        
        if(is.na(groups)) {
    	    grp <- seq_len(nrow(M))
        }
        if(!is.null(SpNames)){
            grp_names <- SpNames
        } else grp_names <- seq_len(nrow(M))
    
        # Save plots of integrated regions and quantification results
        jpegPath <- file.path(DirPath, 
        	                paste(reg[i,"metabolite"], "_", reg[i,"ppm.end"], "_",
                                     reg[i,"ppm.start"], "_ppm", ".jpg", sep=""))
        jpeg(jpegPath, res=300, quality=100, height=8, width=18, units="cm")
        if(isTRUE(resultObject$baseline)){
            par(mfrow=c(1,3))
            matplot(T[,1], T[,2:dim(T)[2]], type="l", lty=1, xlab="chemical shift (ppm)", 
            	ylab="intensity (a.u.)", 
                col = grp, xlim = rev(range(T[,1])), main = "with baseline correction")
            if(is.na(groups)){
                barplot(quan[,i], 
                        main=paste(reg[i,"metabolite"], 
                        " ", reg[i,"ppm"], " ppm"), names=grp_names, 
                xlab="sample index", ylab="Concentration (mM)")
            } else {
                boxplot(quan[,i] ~ grp, 
                        main=paste(reg[i,"metabolite"], 
                        " ", reg[i,"ppm"], " ppm"), names=grp_names, 
                xlab="", ylab="Concentration (mM)", col=unique(grp))
            }
        } else {
            par(mfrow=c(1,2))
            matplot(M[a:b,1], M[a:b,2:dim(M)[2]], type = "l", lty = 1, 
                xlab = "chemical shift (ppm)", col = grp,
                ylab = "intensity (a.u.)", xlim = rev(range(M[a:b,1])), 
                main = "without baseline correction")
        if(is.na(groups)){
                barplot(quan[,i], main=paste(reg[i,"metabolite"], 
                " ", reg[i,"ppm"], " ppm"), names=grp_names, 
                xlab="sample index", ylab="Concentration (mM)")
            } else {
                boxplot(quan[,i] ~ grp, main=paste(reg[i,"metabolite"], 
                " ", reg[i,"ppm"], " ppm"), names=grp_names, 
                xlab="", ylab="Concentration (mM)", col=unique(grp))
            }
        }
    dev.off()
    }
    
    # Save a table with integrals per metabolite (rows) per sample (column)
    csvPath <- file.path(DirPath, 
                paste(resultObject$projectName, 
                "_", "quantification", ".csv", sep=""))
    write.csv(resultObject$quantification, csvPath, quote=FALSE)
}
