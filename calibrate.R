# Function to align spectra according to a particular signal 
# To use after 'readBrukerNMR' function
# f1, f2 left and right margins of the region in ppm 
# containing the spectral reference (ex. ERETIC, TSP, DSS, glucose...)
# Goncalo Graca 1 July 2025

calibrate <- function(X,f1,f2){
  idx <- which(X[,1] > f2 & X[,1] < f1)
  ind <- apply(X[idx,2:ncol(X)], 2, which.max)
  med <- ind[1]
  d <- as.numeric(lapply(ind,function(X) X-med))
  n <- nrow(X)
  p <- ncol(X)
  M <- matrix(rep(NA,times=n*p) , nrow = n)
  colnames(M) <- colnames(X)
  M[,1] <- X[,1]
  for (i in 1:length(d)){
    if (d[i] == 0){
      M[,i+1]<-X[,i+1]
      next
    }
    if (d[i] < 0){
      M[,i+1]<-c(rep(0,times=abs(d[i])),X[1:(dim(X)[1]-abs(d[i])),i+1])
      next
    }
    if (d[i] > 0){
      M[,i+1]<-c(X[(d[i]+1):(dim(X)[1]),i+1],rep(0,times=abs(d[i])))
      next
    }
  }
  return(M)
 }
