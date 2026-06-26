\name{mtbls1}
\alias{mtbls1}
\docType{data}
\title{Processed human urine spectra from MTBLS1 dataset}
\description{
As example data set this tutorial makes use of the publicly available dataset from Salek et al. 2007
\href{https://www.ebi.ac.uk/metabolights/MTBLS1/}{"A metabolomic study of urinary changes in type 2 diabetes in human compared to the control group"},
available through the \href{https://www.ebi.ac.uk/metabolights/}{MetaboLights repository}.
}
\details{
The dataset is composed of 132 1H-NMR spectra from human urine of healthy and diabetic subjects.
The spectra have been aligned using the R package \ref{https://cran.r-project.org/web/packages/speaq/index.html}{speaq} and normalised using the PQN method.
}

\format{A matrix with ppm as the first column and each spectrum in subsequent columns and 28493 rows.}
\usage{data("mtbls1")}
\references{
Salek, R. M., Maguire, M. L., et. al. A metabolomic comparison of urinary changes in type 2 diabetes in mouse, rat, and human. 
\emph{Physiological Genomics} 2007 29:2, 99-108. \ref{https://doi.org/10.1152/physiolgenomics.00194.2006}
}
