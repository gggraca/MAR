# Import spectra

This vignette illustrates how to import 1D 1H-NMR spectra (Bruker format) for further processing.

As example data set this tutorial makes use of the publicly available serum/plasma dataset from the MetaboLightsrepository,
study MTBLS470 ["A Metabonomic Comparison of Human Serum and Plasma Subtypes Using Untargeted NMR Spectroscopy and UPLC-MS"](https://www.ebi.ac.uk/metabolights/editor/MTBLS470). 

From the dataset, 5 1D 1H-NMR spectra were downloaded and unzipped to a local folder ("/MTBLS470").

Firstly we start by laoding the package:
```
library(MAR)
```
Next we load the data:
```
data(mtbls1)
```
The first step to import the data is to set the correct file paths to the processed data ("/pdata/" folder in the Bruker NMR folder structure).
We do this by first listing all files under the "./MTBLS470/" folder and then selecting all '1r' file paths, which are then edited to obtain 
the "./pdata/" file paths for all samples:
```
spaths <- dir("./MTBLS470/", full.names=TRUE, recursive=TRUE)
idx <- grep("1r", spaths)
spaths <- spaths[idx]
spaths <- sapply(spaths, function(x) substr(x,1,nchar(x)-2))
spaths <- as.character(spaths)
```

The spectra are then imported using the "read1DNMR" function.
The option "simplify=TRUE" returns a matrix with ppm scale in the first column and each spectrum intensity as the remaining columns:
```
spectra <- read1DNMR(spaths,  xmax=18, xmin=-1, npoints=64000, simplify=TRUE)
```
If the option "simplify" is set to FALSE, the data will be returned as a list containing the ppm scale and sample names (from the spectra title, if available) 
as a vectors and the matrix of spectral intensities for all the samples.
The spectra limits (xmax and xmin) and number of points (npoints) are selected by the user according to the type of spectra and acquisition
parameters. In this example npoints=64000 is sufficient to maintain spectral resolution.

The spectra might need to be re-calibrated after being imported. This should be checked by plotting several regions of the spectrum,
e.g. the region of the alpha-glucose anomeric proton doublet:
```
plotNMR(spectra, 5.3, 5.2)
```

To recalibrate the spectra, one singlet signal or one isolated signal of a multiplet can be selected as reference for re-calibration. 
The "calibrate" function then calibrates all the spectra to the median chemical shift position of the peak appexes from all samples.
In this example, we use the alpha-glucose anomeric proton doublet around 5.25 ppm:

```
s <- calibrate(s, 5.28, 5.23)
```
The chemical shift scale might require aditional shifting for a more precise calibration.
Here it was necessary to shift the ppm scale by 0.03 ppm:

```
s[,1] <- s[,1] - 0.03
```

Finally, we can add the sample names and ppm to the columns of the spectra matrix:
```
colnames(s) <- c("ppm","10","30","40","50","70")
```
The spectra matrix is now ready for further processing.