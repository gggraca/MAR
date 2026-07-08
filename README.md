# Metabolomics And nmR (MAR) toolbox

A collection of functions and scripts to read, process and analyse proton 1D NMR spectra.

Some examples of application of the code from the toolbox can be found in the published papers. <sup>[1](#myfootnote1),</sup> <sup>[2](#myfootnote2),</sup> <sup>[3](#myfootnote3)</sup> (see references).

<h3> Description </h3>

Main functionalities include data import from [Bruker](https://www.bruker.com/en/products-and-solutions/mr.html) format, preprocessing (trimming, normalisation and scaling), peak integration and quantification. 

<h3> Installation </h3>

To install the 'Metabolomics And nmR toolbox' use the commands below:
```
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("https://github.com/gggraca/MAR/")
```
The package also depends on functionality from the 'ptw' R package to perform baseline corrections and some plotting functionality from the 'ggplot2' package. Both packages should be automatically installed, otherwise they can be installed with:
```
install.packages("ptw")
install.package("ggplot2")
```

<h3> Examples </h3>

Example of how to import 1D NMR data can be found in the [import_data](https://github.com/gggraca/MAR/blob/main/vignettes/import_data.md) vignette.
A tipical example of data integration and quantification can be found in the vignettes folder as [integration_examples](https://github.com/gggraca/MAR/blob/main/vignettes/integration_examples.md). 


<h3> References </h3>

<a name="myfootnote1">1</a>. Akhbari, P. et al. Differences between infected and noninfected synovial fluid: An observational study using metabolic phenotyping (or “metabonomics"). Bone Joint Res, 2021;10(2):85–95. https://doi.org/10.1302/2046-3758.101.BJR-2020-0285.R1</a>

<a name="myfootnote2">2</a>. Akhbari, P. et al. Differences in the composition of hip and knee synovial fluid in osteoarthritis: a nuclear magnetic resonance (NMR) spectroscopy study of metabolic profiles. Osteoarthritis and Cartilage, Volume 27, Issue 12, 1768 - 1777. https://doi.org/10.1016/j.joca.2019.07.017 </a>

<a name="myfootnote3">3</a>. Graça, G., Desterro, J., Sousa, J. et al. Identification of putative biomarkers for leptomeningeal invasion in B-cell non-Hodgkin lymphoma by NMR metabolomics. Metabolomics 13, 136 (2017). https://doi.org/10.1007/s11306-017-1269-9 </a>
