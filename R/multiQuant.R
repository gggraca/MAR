#' @title Integration and quantification of 1H-NMR spectra
#' 
#' @description
#' R function for baseline correction, integration and quantification
#' The \code{ptw} package is required to perform baseline correction. 
#' The main function performs integration on multiple spectra regions with or
#' without baseline correction and calculates the respective metabolite 
#' concentration given a quantification reference (e.g. TSP or ERETIC signal).
#' 
#' @author Goncalo Graca, Imperial College London 
#'
#' @param M Matrix of 1H-NMR spectra where samples are arranged in rows 
#' and chemical shift values (ppm) in columns. The matrix will be transposed
#' otherwise. Sample names should not be included in this matrix, only ppm
#' values.
#' @param reg Data frame containing as columns: \code{metabolite}, \code{ppm}, 
#' number of protons corresponding to the signal \code{nprotons}, start and 
#' end of the signal in ppm (\code{ppm.start} and \code{ppm.end}, respectively)
#' @param reference Exact name of the reference metabolite or compound, as
#' named in the \code{reg} data frame.
#' @param refConc The concentration of the reference compound/metabolite in mM
#' @param spNames A string containing the sample names. The default 
#' is \code{NULL}.
#' @param baseline A logical value to specify if baseline should be performed.
#' The default value is \code{TRUE}.
#' @param groups A string of values or characters that specify any sample 
#' groups (e.g. controls and disease). The default is \code{NA}.
#' @projectName The name of the project, which will be appended to the results
#' files upon saving.
#'
#' @return A list containing containing several objects: 
#' the spectra matrix \code{spectra}, data frame containing the integration 
#' regions, the number of protons and chemical shift of each metabolite signal,
#' \code{integrationRegions}, the integration and quantification results 
#' (\code{integrals} and \code{quantification}) and other the metadata such as 
#' \code{baseline}, reference name(\code{reference}) and concentration 
#' (\code{refConc}), the groups labels \code{groups}, sample names 
#' \code{spNames}, and project name \code{projectName}.
#' 
#' @example ##
#' @export
multiQuant <- function(M, reg, reference="ERETIC", refConc=10, spNames=NULL, 
                       baseline=TRUE, groups=NA, projectName="") {
    #if the matrix is not in the format [variables,samples] it will be transposed
    if(nrow(M) < ncol(M)){
        M <- t(M)
    }
	# create matrix to store integration results
    intg <- matrix(NA, ncol = ncol(M)-1, nrow = nrow(reg))
    # integrate all regions first
    if(baseline){
        for(i in 1:nrow(reg)){
            intg[i,] <- bintegral(M, reg[i,"ppm.start"], reg[i,"ppm.end"])
        }
    }
    if(!baseline){
        for(i in 1:dim(reg)[1]){
            intg[i,] <- integral(M, reg[i,"ppm.start"], reg[i,"ppm.end"])
        }
    }
    # perform quantification given a reference and 
    # the number of protons per metabolite
    refIdx <- grep(reference, reg$metabolite) 
    
    # run the quantification
    quan <- lapply(seq_len(nrow(intg)), function(x) {
    	(refConc * intg[x,]) / (reg[x,"nprotons"] * intg[refIdx,])
    })
    # convert from list to matrix
    quan <- do.call(rbind,quan)
    # remove the reference from the result
    quan <- quan[-refIdx,]
    
    # transpose both integrals and quantification matrices
    intg <- t(intg)
    quan <- t(quan)
    # add metabolites and sample name
    colnames(intg) <- reg[,"metabolite"]
    colnames(quan) <- reg[,"metabolite"]
    if(is.null(SpNames)){
    	rownames(intg) <- colnames(M)[2:ncol(M)]
        rownames(quan) <- colnames(M)[2:ncol(M)]
    } else {
    	rownames(intg) <- SpNames
        rownames(quan) <- SpNames
    }
    
    result <- list(spectra=M, integrationRegions=reg, 
                   integrals=intg, quantification=quan, baseline,
                   reference, refConc, groups, spNames, projectName)
    return(result)
}
