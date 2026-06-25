# read from Bruker raw processed data
# filepaths: full datapaths to the sample/exp_no/pdata/proc_no
# Goncalo Graca, 25 September 2020, g.gomes-da-graca@imperial.ac.uk
readBrukerNMR <- function(filepaths){
  procs <- readLines(paste(filepaths[1],"/","procs", sep = ""))
  OFFSET <- procs[grep('OFFSET', procs)]
  OFFSET <- substr(OFFSET, start = 12, stop = nchar(OFFSET))
  OFFSET <- as.numeric(OFFSET)
  
  SW <- procs[grep('SW_p', procs)]
  SW <- substr(SW, start = 10, stop = nchar(SW))
  SW <- as.numeric(SW)
  
  SF <- procs[grep('SF=', procs)]
  SF <- substr(SF, start = 8, stop = nchar(SF))
  SF <- as.numeric(SF)
  
  SI <- procs[grep('$SI=', procs, fixed=TRUE)]
  SI <- substr(SI, start=8, stop = nchar(SI))
  SI <- as.numeric(SI)
  
  SWP <- SW/SF
  dppm <- SWP/SI
  ippm <- OFFSET-SWP 
  ppm <- seq(ippm , OFFSET, dppm)
  if(length(ppm) > SI) ppm <- ppm[-1]
  
  # build spectra matrix
  counter <- length(filepaths)
  m <- matrix(NA, nrow = SI, ncol = counter+1) 
  m[,1] <- ppm
  
  for (i in 1:counter){
    filepath <- paste(filepaths[i], "/1r", sep = "")
    procs <- readLines(paste(as.character(filepaths[i]),"/","procs", sep = ""))
    to.read <-  file(filepath,"rb")
    tmp <- readBin(to.read, integer(), n = SI, endian = "little")
    close(to.read)
    NCproc <- procs[grep('NC_proc', procs)]
    NCproc <- substr(NCproc, start = 13, stop = nchar(NCproc))
    NCproc <- as.numeric(NCproc)
    m[,i+1] <- tmp/2^-NCproc
  }
  m[,1] <- rev(ppm)
  return(m)
} 