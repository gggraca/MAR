# R functions for baseline correction, integration and quantification
# ptw package is required to perform baseline correction
# Multi-integration results are stored as a list
# Goncalo Graca, Imperial College London
# g.gomes-da-graca@imperial.ac.uk
# 23 June 2026
# to work with NMR data in table format with samples in rows and ppms in columns
# the first row contains the chemical shift values

# function to integrate multiple regions with or without baseline correction
# and calculate the concentration given a quantification reference

multiQuant <- function(M, reg, reference="ERETIC", refConc=1, spNames=NULL, 
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
    quan <- lapply(seq_len(nrow(intgr)), function(x) {
    	(refConc * intg[x,]) / (reg[x,"nprotons"] * intg[refIdx,])
    })
    # convert from list to matrix
    quan <- do.call(rbind,quan)
    # remove the reference from the result
    quan <- quan[-refIdx,]
    # add metabolites and sample name
    rownames(quan) <- reg[,"metabolite"]
    if(is.null(SpNames)){
        colnames(quan) <- colnames(M)[2:ncol(M)]
    } else colnames(quan) <- SpNames
    
    result <- list(spectra=M, integrationRegions=reg, 
                   integrals=intg, quantification=quan, baseline=baseline, 
                   groups=groups, projectName=projectName)
}


# save stacked plots of the integration windows and boxplots with the peak area values
# the code below takes the indexes of the chemical shifts of the integration/baseline regions:
saveQuantification <- function(resultObject){
    M <- resultObject$spectra
    # if the matrix is not in the format [variables,samples] it will be transposed
    if(nrow(M) < ncol(M)){ 
        M <- t(M)
    }
    reg <- resultObject$integrationRegions
    ints <- which(M[,1] > reg[i,"ppm.end"] & M[,1] < reg[i,"ppm.start"])
    a <- ints[1]
    b <- ints[length(ints)]
    if(isTRUE(resultObject$baseline)){
        T <- bas(M, reg[i,"ppm.start"], reg[i,"ppm.end"])
    }
    intg <- resultObject$intg
    if(is.na(groups)) {
    	grp <- seq_len(nrow(M))
    }
    if(!is.null(SpNames)){
        grp_names <- spNames
    } else grp_names <- seq_len(nrow(M))
    
    # plot results
    jpeg(paste(reg[i,"metabolite"], "_", reg[i,"ppm.end"], "_", 
            reg[i,"ppm.start"], "_ppm", ".jpg", sep=""), res=300, quality=100, 
            height=8, width=18, units="cm")
        if(isTRUE(resultObject$baseline)){
            par(mfrow=c(1,3))
            matplot(T[,1], T[,2:dim(T)[2]], type="l", lty=1, xlab="chemical shift (ppm)", ylab="intensity (a.u.)", 
                col = grp, xlim = rev(range(T[,1])), main = "with baseline correction")
            if(is.na(grp)){
                barplot(intg[i,], main=paste(reg[i,1], " ", reg[i,2], " ppm"), names=grp_names, 
                xlab="", ylab="Concentration (mM)")
            } else {
                boxplot(intg[i,] ~ grp, main=paste(reg[i,1], " ", reg[i,2], " ppm"), names=grp_names, 
                xlab="", ylab="Concentration (mM)", col=unique(grp))
            }
        } else {
            par(mfrow=c(1,2))
            matplot(M[a:b,1], M[a:b,2:dim(M)[2]], type = "l", lty = 1, xlab = "chemical shift (ppm)", col = grp,
                ylab = "intensity (a.u.)", xlim = rev(range(M[a:b,1])), main = "without baseline correction")
        if(is.na(grp)){
                barplot(intg[i,], main=paste(reg[i,1], " ", reg[i,2], " ppm"), names=grp_names, 
                xlab="", ylab="Concentration (mM)")
            } else {
                boxplot(intg[i,] ~ grp, main=paste(reg[i,1], " ", reg[i,2], " ppm"), names=grp_names, 
                xlab="", ylab="Concentration (mM)", col=unique(grp))
            }
        }
    dev.off()
    # Save a table with integrals per metabolite (rows) per sample (column)
    fileName <- paste(resultObject$projectName, "_", "quantification", ".csv", sep="")
    write.csv(resultObject$quantification, fileName, quote=FALSE)
}
