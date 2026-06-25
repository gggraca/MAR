## R functions for baseline correction and integration
##
## Requires the 'ptw' package to perform baseline correction
## Multi-integration results are stored as a list
## Goncalo Graca, Imperial College London
##
## to work with NMR data in matrix format with samples in rows and 
## chemical shift (ppm) as columns
## the first row should contain the chemical shift values, otherwise
## the spectra matrix will be transposed
##
## These functions are not exported

## function for integration using trapezoid rule-------------------------------
integral <- function(X,f1,f2) {
# if the matrix is not in the format [variables,samples] it will be transposed
    if(dim(X)[1] < dim(X)[2]){ 
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


## function to perform baseline correction only--------------------------------
bas <- function(X, f1, f2) {
#if the matrix is not in the format [variables,samples] it will be transposed
    if(dim(X)[1] < dim(X)[2]){
        X<-t(X)
    }
        lim <- which(X[,1] > f2 & X[,1] < f1)
        int1 <- lim[1]
        int2 <- lim[length(lim)]
        ints <- seq(from = int1, to = int2, by = 1)
        M<- apply(X[ints,2:dim(X)[2]], 2, baseline.corr)
        M <- cbind(X[ints,1], M)
        return(M)
}

## function to perform baseline correction 
## followed by integration using the trapezoid rule----------------------------
bintegral <- function(X,f1,f2) {
# if the matrix is not in the format [variables,samples] it will be transposed
    if(dim(X)[1] < dim(X)[2]){ 
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