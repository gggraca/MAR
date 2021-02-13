# R functions for baseline correction and integration
# ptw package is required to perform baseline correction
# Goncalo Graca, Imperial College London
# g.gomes-da-graca@imperial.ac.uk
# 13 February 2018
# to work with NMR data in table format with samples in rows and ppms in columns
# the first row contains the chemical shift values

# function for integration using trapezoid rule
integral <- function(X,f1,f2) {
  if(dim(X)[1] < dim(X)[2]){ #if the matrix is not in the format [variables,samples] it will be transposed
    X<-t(X)
  }
        ints <- which(X[,1] > f2 & X[,1] < f1)
        M <- X[ints,2:dim(X)[2]]
        s <- rep(1,length(2:dim(X)[2]))
        for(i in 1:dim(M)[2]){
                y <- M[,i]
                n <- length(y)-1
                x <- X[ints,1]
                if(x[1] > x[n+1]) {
                  b <- x[1]
                  a <- x[n+1]
                } else {
                  a <- x[1]
                  b <- x[n+1]
                }
                h <- (b-a)/n
                s[i] <- h * (y[1]/2 + sum(y[2:n]) + y[n+1]/2)
        }
        return(s)
}


# function to perform baseline correction only
bas <- function(X, f1, f2) {
  if(dim(X)[1] < dim(X)[2]){ #if the matrix is not in the format [variables,samples] it will be transposed
    X<-t(X)
  }
        lim<-which(X[,1] > f2 & X[,1] < f1)
        int1<-lim[1]
        int2<-lim[length(lim)]
        ints<-seq(from = int1, to = int2, by = 1)
        M<-apply(X[ints,2:dim(X)[2]], 2, baseline.corr)
        M<-cbind(X[ints,1], M)
        return(M)
}

# function to perform baseline correction followed by integration using the trapezoid rule
bintegral <- function(X,f1,f2) {
  if(dim(X)[1] < dim(X)[2]){ #if the matrix is not in the format [variables,samples] it will be transposed
    X<-t(X)
  }
        ints <- which(X[,1] > f2 & X[,1] < f1)
        M <- apply(X[ints, 2:dim(X)[2]], 2, baseline.corr)
        s <- M[1,]
        for(i in 1:dim(M)[2]){
                y <- M[,i]
                n <- length(y)-1
                x <- X[ints,1]
                if(x[1] > x[n+1]) {
                  b <- x[1]
                  a <- x[n+1]
                } else {
                  a <- x[1]
                  b <- x[n+1]
                }
                h <- (b-a)/n
                s[i]<- h * (y[1]/2 + sum(y[2:n]) + y[n+1]/2)
        }
        return(s)
}

# function to plot stacked spectral regions
plotNMR <- function(X,f1,f2){
  if(dim(X)[1] < dim(X)[2]){ #if the matrix is not in the format [variables,samples] it will be transposed
    X<-t(X)
  }
  ints <- which(X[,1] > f2 & X[,1] < f1)
  a <- ints[1]
  b <- ints[length(ints)]
  matplot(X[a:b,1],X[a:b,2:dim(X)[2]],type="l",lty=1,xlab='chemical shift (ppm)',ylab='intensity (a.u.)',xlim=rev(range(X[a:b,1])))
}

# function to integrate and plot multiple regions with or without baseline correction
# jpeg figures will saved on the workspace folder corresponding to each integrated peak, before and after baseline correction
multiIntegration <- function(M, reg, plots = TRUE, baseline = TRUE, grp, grp_names, save.results = TRUE) {
  if(nrow(M) < ncol(M)){ #if the matrix is not in the format [variables,samples] it will be transposed
    M <- t(M)
  }
  intg <- matrix(NA, ncol = ncol(M)-1, nrow = nrow(reg)) # create matrix to store integration results
  if(baseline){
    for(i in 1:nrow(reg)){
      intg[i,] <- bintegral(M, reg[i,2], reg[i,3])
      ints <- which(M[,1] > reg[i,3] & M[,1] < reg[i,2])
      a <- ints[1]
      b <- ints[length(ints)]
      T <- bas(M, reg[i,2], reg[i,3])
      if(plots){
        jpeg(paste(reg[i,1], "_", reg[i,3], "_", reg[i,2], "_ppm", ".jpg", sep = ""), res = 300, quality = 100, 
             height = 8, width = 18, units = "cm")
        par(mfrow=c(1,3))
        matplot(M[a:b,1], M[a:b,2:dim(M)[2]], type = "l", lty = 1, xlab = "chemical shift (ppm)", col = grp,
                ylab = "intensity (a.u.)", xlim = rev(range(M[a:b,1])), main = "without baseline correction")
        matplot(T[,1], T[,2:dim(T)[2]], type = "l", lty = 1, xlab = "chemical shift (ppm)", ylab = "intensity (a.u.)", 
                col = grp, xlim = rev(range(T[,1])), main = "with baseline correction")
        boxplot(intg[i,] ~ grp, main = paste(reg[i,1], " ", reg[i,2], " ppm"), names = grp_names, 
                xlab = "", ylab = "peak area (a.u.)", col = unique(grp))
        dev.off()
        }
      }
  }
  if(!baseline){
    for(i in 1:dim(reg)[1]){
      intg[i,] <- integral(M, reg[i,2], reg[i,3])
      ints <- which(M[,1] > reg[i,3] & M[,1] < reg[i,2])
      a <- ints[1]
      b <- ints[length(ints)]
      if(plots){
        jpeg(paste(reg[i,1], "_", reg[i,2], "_", reg[i,3], "_ppm", ".jpg", sep = ""), res = 300, quality = 100, 
             height = 8, width = 16, units = "cm")
        par(mfrow = c(1,2))
        matplot(M[a:b,1], M[a:b,2:dim(M)[2]], type = "l", lty = 1, xlab = "chemical shift (ppm)", col = grp,
                ylab = "intensity (a.u.)", xlim = rev(range(M[a:b,1])), main = "")
        boxplot(intg[i,] ~ grp, main = paste(reg[i,1], " ", reg[i,2], " ppm"), names = grp_names,
                xlab = "", ylab = "peak area (a.u.)", col = unique(grp))
        dev.off()
      }
    }
  }
  # A table was generated with integrals per metabolite (rows) per sample (column)
  # Next line gives the names to the rows by getting the metabolite names from the regions table
  rownames(intg) <- reg[,1]
  colnames(intg) <- colnames(M)[2:ncol(M)]
  assign("integrals", intg, envir = .GlobalEnv)
  write.csv(intg, "results_integrals.csv", quote = FALSE)
}