#' @title Integration of 1H-NMR spectra
#' 
#' @description
#' R function for baseline correction and integration of NMR peaks.
#' The \code{ptw} package is required to perform baseline correction. 
#' The main function performs integration on multiple spectra regions with or
#' without baseline correction.
#' 
#' @author Goncalo Graca, Imperial College London 
#'
#' @param M Matrix of 1H-NMR spectra where samples are arranged in rows 
#' and chemical shift values (ppm) in columns. The matrix will be transposed
#' otherwise. Sample names should not be included in this matrix, only ppm
#' values.
#' @param reg Data frame containing as columns: \code{metabolite}, the
#' chemical shift of the signal \code{ppm}, number of protons corresponding 
#' to the signal \code{nprotons} (optional), start and 
#' end of the signal in ppm (\code{ppm.start} and \code{ppm.end}, respectively)
#' @param spNames A string containing the sample names. The default 
#' is \code{NULL}.
#' @param baseline A logical value to specify if baseline should be performed.
#' The default value is \code{TRUE}.
#' @param groups A string of values or characters that specify any sample 
#' groups (e.g. controls and disease). The default is \code{NA}.
#' @param projectName The name of the project, which will be appended to the results
#' files upon saving.
#'
#' A list containing containing several objects: 
#' the spectra matrix \code{spectra}, data frame containing the integration 
#' regions, the number of protons and chemical shift of each metabolite signal,
#' \code{integrationRegions}, the integration results (\code{integrals} and 
#' other the metadata such as \code{baseline}, the groups labels \code{groups},
#'  sample names \code{spNames}, and project name \code{projectName}.
#' 
#' @examples
#' data("mtbls1")
#' reg <- data.frame(metabolite = c("Creatine", "Creatinine", "Glucose"), 
#'                   nprotons = c(3, 2, 1), ppm=c(3.04, 4.06, 5.23), 
#'                   ppm.start=c(3.046, 4.068, 5.26), 
#'                    ppm.end=c(3.040, 4.054, 5.23))
#' integrationResult <- multiIntegration(urine, reg, baseline=TRUE, SpNames=NULL, 
#'                                     groups=NA, projectName="MTBLS1")
#' @export
multiIntegration <- function(M, reg, baseline=TRUE, SpNames=NULL, 
                            groups=NA, projectName="") {
    #if the matrix is not in the format [variables,samples] it will be transposed
    if(nrow(M) < ncol(M)){
        M <- t(M)
    }
	# create matrix to store integration results
    intg <- matrix(NA, ncol = ncol(M)-1, nrow = nrow(reg))
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

    # transpose integrals table to samples as rows and metabolites as columns
    intg <- t(intg)
    # add metabolites and sample names (if available)
    colnames(intg) <- reg[,"metabolite"]
    if(is.null(SpNames)){
        rownames(intg) <- colnames(M)[2:ncol(M)]
    } else colnames(intg) <- SpNames
    
    intg <- data.frame(intg)
    
    result <- list(spectra=M, integrationRegions=reg, 
                   integrals=intg, baseline=baseline, SpNames=SpNames,
                   groups=groups, projectName=projectName)
    return(result)
}