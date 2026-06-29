# Spectral integration

The following text demonstrates how to perform numeric integration (using the trapezoid method<sup>[1](#References)</sup>) of defined NMR spectral peaks.

As example data set this tutorial makes use of the publicly available dataset from Salek *et al.* 2007 <sup>[2](#References)</sup>
["A metabolomic study of urinary changes in type 2 diabetes in human compared to the control group"](https://www.ebi.ac.uk/metabolights/MTBLS1/),
available through the [MetaboLights](https://www.ebi.ac.uk/metabolights/) repository. 

The dataset is composed of 132 <sup>1</sup>H-NMR spectra from human urine of healthy and diabetic subjects.
The version of the data used in this vignette has already been aligned using the R package ['speaq'](https://cran.r-project.org/web/packages/speaq/index.html) and normalised using the PQN method<sup>[3](#References)</sup>.
The processed data is included as demo data for the package.

Firstly we start by laoding the package:
```
library(MAR)
```
Next we load the data:
```
data(mtbls1)
```
Stacked spectra can be plotted to check the data was correctly loaded using the toolbox function plotNMR:
```
plotNMR(mtbls1)
```

Now it is necessary to create a table that contains the metabolite regions to be integrated and/or quantitated. The table should have the following format:
Metabolite | nprotons| ppm | ppm.start | ppm.end
---------- | --------|-----|-----------|--------
metabolite 1 |   2   | 7.5 |    8.0    |   7.0  
metabolite 2 |   1   | 2.8 |    3.0    |   2.6
.............|.......|.....|...........|........
metabolite n |   3   | 0.5 |    1.0    |   0.0

This table can either be created in R or read from an external file.
Here is an example of table in R with three metabolites that will be integrated:
```
reg <- data.frame(metabolite = c("Creatine", "Creatinine", "Glucose"), 
                   nprotons = c(3, 2, 1), ppm=c(3.04, 4.06, 5.23), 
                   ppm.start=c(3.046, 4.068, 5.26), 
                   ppm.end=c(3.040, 4.054, 5.23))
```

It is then possible to run the integration on the 3 peaks defined in "reg" across the 132 samples from the dataset.
The function "multiIntegration" will integrate the peaks defined in the 'reg' data frame. This can be done with or without baseline correction.
Optionally the groups of samples (if any) can also be defined and used as metadata. Similarly, the sample names can also be added for later use. 
In the example, we will define the sample groups, in this case "diabetes" and healthy controls:

```
grp <- c(rep("Diabetes", 48), rep("Control", 84))
```

The integration function can now be run using the following command:
```
results <- multiIntegration(mtbls1, reg, baseline=TRUE, SpNames=NULL, 
                            groups=grp, projectName="MTBLS1")
```
The function will calculate the integrals after local baseline correction, and store the integration results and metadata on the variable 'results'.

The results from the integration, as well as stacked plots for the integrated regions and boxplots can now be exported to the project folder:
```
saveIntegration(results, DirPath="./project")
```

# Metabolite quantification

It is also possible to perform metabolite quantification if a quantification reference compound with known concentration is present
in the spectra (note that absolute quantification is only reliable if full nuclei relaxation is achieved). 
Using the data from the example above, it can consider glucose as the metabolite with known concentration (e.g. 1 mM) and
use it has reference in the quantification in the following way:

```
quanResult <- multiQuant(urine, reg, reference="Glucose", refConc=1, removeRef=TRUE,
                        SpNames=NULL, baseline=TRUE, groups=grp, 
                        projectName="MTBLS1")
```

To save the quantification results, the stacked plots for the integrated regions and boxplots with the concentration values to the project folder, 
the following command can now applied:
```
saveQuantification(quanResults, DirPath="./project")
```

<h3>References</h3>

1. Trapezoid Rule: https://en.wikipedia.org/wiki/Trapezoidal_rule

2. Salek, R. M., Maguire, M. L., Bentley, E., Rubtsov, D. V., Hough, T., Cheeseman, M., Nunez, D., Sweatman, B. C., Haselden, J. N., Cox, R. D., Connor, S. C., Griffin, J. L. A metabolomic comparison of urinary changes in type 2 diabetes in mouse, rat, and human. Physiological Genomics 2007 29:2, 99-108. https://doi.org/10.1152/physiolgenomics.00194.2006
         
3. Dieterle F, Ross A, Schlotterbeck G, Senn H. Probabilistic quotient normalization as robust method to account for dilution of complex biological mixtures. Application in 1H NMR metabonomics. Anal Chem. 2006 Jul 1;78(13):4281-90.  https://doi.org/10.1021/ac051632c
