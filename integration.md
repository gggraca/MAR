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
Stacked spectra can be plotted to check the data was correctly loaded using the toolbox function plotNMR:
```
plotNMR(spec, 10, 0)
```
![Stacked Spectra](https://github.com/gggraca/MAR/images/stacked_urine.png)

Now it is necessary to create a table that contains the metabolite regions to be integrated. The table should have the followin format:
Metabolite | ppm.start | ppm.end
---------- | ----------|---------
metabolite 1 | 8.0 | 7.0  
metabolite 2 | 4.2 | 3.0 
...|...|...
metabolite n | 1.5 | 0.5

This table can either be created in R or read from an external file.
Here is an example of table in R with three metabolites that will be integrated:
```
reg <- data.frame(metabolite = c("Creatine", "Creatinine", "Glucose"), 
    ppm.start = c(3.046,4.068, 5.26), 
    ppm.end = c(3.040,4.054, 5.23)) 
```
It is then possible to run the integration on the 3 peaks defined in "reg" across the 132 samples from the dataset.
The function "multiIntegration" from the toolbox will integrate the peaks, with or without baseline correction, and optionally it can save the plots of stacked spectral peaks and boxplots coloured by group. 
The following example will use these two options for plotting, hence it will be necessary to define the group variable which contains the colours. In this particular case, the dataset is composed of 48 diabetic subjects and 84 healthy controls. Then the group variable "grp" is defined as: 
```
grp <- c(rep("red", 48), rep("blue", 84))
```
A variable "grp_names" must also be defined to include the group names:
```
grp_names <- c("Diabetic", "Healthy Controls")
```
Then run the integration function:
```
multiIntegration(spec, reg, plots = TRUE, baseline = TRUE, grp, grp_names, save.results = TRUE)
```
The function will calculate the integrals, create an "integrals" object containing the integration results and will also save the results as a .csv table in the working directory.

The generated plots for the example are shown below:
