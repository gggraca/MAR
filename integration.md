# Spectral integration

The following text demosntrate how to perform numeric integration of defined NMR spectral peaks.
It requires the functions contained in the file [integration_functions.R](https://github.com/gggraca/MAR/integration_functions.R) and also a function from the R package ['ptw'](https://cran.r-project.org/web/packages/ptw/index.html).

As example data set this tutorial makes use of the publically available dataset from Salek *et al.* 2007 
["A metabolomic study of urinary changes in type 2 diabetes in human compared to the control group"](https://www.ebi.ac.uk/metabolights/MTBLS1/),
available in the [MetaboLights](https://www.ebi.ac.uk/metabolights/) repository. 

The dataset is composed of 132 <sup>1</sup>H-NMR spectra from human urine of healthy and diabetic subjects.
The version of the data used here has already been aligned using the R package ['speaq'](https://cran.r-project.org/web/packages/speaq/index.html) and normalised using the PQN method.

Firstly, install the ptw package and load the [integration_functions.R](https://github.com/gggraca/MAR/integration_functions.R), which should be stored in the working directory.
```
install.packages("ptw")
source("integration_functions.R")
```

Next load the data table containing the preprocessed spectra:
```
spec <- read.csv("MTBLS1_aligned_normalised.csv", header = TRUE)
```
