# Functions to import 1d spectra from Bruker format
# interpolate to get same chemical shift scale and export as data matrix
# Goncalo Graca, 1 July 2025

# Main function:
read1DNMR <- function(filepaths,  xmax, xmin, npoints = 64000, simplify = TRUE){
  s <- lapply(filepaths, readBruker1D)
  f <- interpolate_1DNMR(s, xmax, xmin, npoints = 64000)
  if(simplify){
    r <- rbind(f$axis, f$data_1r)
    r <- t(r)
  } else {
    r <- f
  }
  return(r)
}


# Auxiliary functions:

# read from processed 1r Bruker Topspin files and metadata 
readBruker1D <- function(filepath) {
  # this reads the parameters need to get the intensities and ppm scale
  # filepath <- paste(filepath, "pdata/1/", sep = "")
  procs <- readLines(paste(filepath,"procs", sep = ""))
  title <- readLines(paste(filepath,"title", sep = ""), warn = FALSE)
  # if(nchar(title) == 0) title <- filepath
  
  # build the ppm scale (direct dimension F2)
  OFFSETF2 <- procs[grep('OFFSET', procs)]
  OFFSETF2 <- substr(OFFSETF2, start = 12, stop = nchar(OFFSETF2))
  OFFSETF2 <- as.numeric(OFFSETF2)
  
  SWF2 <- procs[grep('SW_p', procs)]
  SWF2 <- substr(SWF2, start = 10, stop = nchar(SWF2))
  SWF2 <- as.numeric(SWF2)
  
  SF <- procs[grep('SF=', procs)]
  SF <- substr(SF, start = 8, stop = nchar(SF))
  SF <- as.numeric(SF)
  
  SIF2 <- procs[grep('$SI=', procs, fixed=TRUE)]
  SIF2 <- substr(SIF2, start=8, stop = nchar(SIF2))
  SIF2 <- as.numeric(SIF2)
  
  SWPF2 <- SWF2/SF
  dppm <- SWPF2/SIF2
  ippm <- OFFSETF2-SWPF2
  ppm <- seq(ippm , OFFSETF2, dppm)
  if(length(ppm) > SIF2) ppm <- ppm[-1]
  
  # get intensities
  spath <- paste(filepath, "1r", sep = "")
  to.read <-  file(spath,"rb")
  m <- readBin(to.read, integer(), n = SIF2, endian = "little")
  close(to.read)
  NCproc <- procs[grep('NC_proc', procs)]
  NCproc <- substr(NCproc, start = 13, stop = nchar(NCproc))
  NCproc <- as.numeric(NCproc)
  m <- m/2^-NCproc
  ppm <- rev(ppm)
  
  # intensity matrix is transposed so it has increments (F1) as rows
  result <- list(intensities = m, axis = ppm, title = title)
  return(result)
}

# interpolate spectra to get same chemical shift
interpolate_1DNMR <- function(spectra_object, xmax, xmin, npoints = 64000){
  
  intensities <- NULL
  spectra_names <- NULL
  
  for (i in 1:length(spectra_object)){
    ppm <- spectra_object[[i]]$axis
    I <- spectra_object[[i]]$intensities
    xmax_point <- which.min(abs(ppm - xmax))
    xmin_point <- which.min(abs(ppm - xmin))
    spectrum <- spline(ppm[xmin_point:xmax_point],
                       I[xmin_point:xmax_point],
                       xmin = xmin, xmax = xmax, n = npoints)
    intensities <- rbind(intensities, spectrum$y)
    spectra_names <- c(spectra_names, spectra_object[[i]]$title)
  }
  
  axis <- spectrum$x
  result <- list(data_1r = intensities, axis = axis, titles = spectra_names)
  return(result)
}
