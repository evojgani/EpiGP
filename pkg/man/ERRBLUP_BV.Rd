% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ERRBLUP_bivar.R
\name{ERRBLUP_BV}
\alias{ERRBLUP_BV}
\title{Bivariate ERRBLUP Phenotype Prediction Function}
\usage{
ERRBLUP_BV(M, y1, y2, iters = 20, tolparinv = 1e-06, cores = 1)
}
\arguments{
\item{M}{The original marker matrix of \code{{-1,0,1}}, \code{{0,1}}, \code{{0,1,2}}, \code{{0,2}} or character coded markers with named individuals in the rows and the markers in the columns}

\item{y1}{A numeric vector of named phenotypes for the target triat to be predicted in bivariate ERRBLUP}

\item{y2}{A numeric vector of named phenotypes for the additional trait to be estimated in bivariate ERRBLUP}

\item{iters}{Maximum number of iterations allowed with the default value of \code{20} (The parameter of mmer function from sommer package)}

\item{tolparinv}{Tolerance parameter for matrix inverse used when singularities are encountered in the estimation procedure with the default value of \code{1e-06} (The parameter of mmer function from sommer package)}

\item{cores}{The number of cores with the default value of \code{1}}
}
\value{
A list of three components:

\describe{
  \item{Recodedmarkers}{A \code{{0,1,2}} or \code{{0,2}} coded marker matrix}
  \item{Relationshipmatrix}{A list of two components: ERRBLUP relationship matrix (G) and a vector of all genotype combinations frequencies in the population (P)}
  \item{Predictions}{A dataframe of both phenotype predictions of target triat test set, phenotype estimations of target triat training set and phenotype estimations of additional trait full set in bivarite ERRBLUP method}
}
}
\description{
Function to do phenotype prediction based on all pairwise SNP interactions in bivariate model
}
\examples{
library(BGLR)
data(wheat)
y1 <- Phenotype_BV[,1]
N <- length(y1)
n <- 60
test <- sample(1:N,n)
y1[test] <- NA
y2 <- Phenotype_BV[,2]
M <- wheat.X
rownames(M) <- names(y1)
ERRBLUP_Bivar <- ERRBLUP_BV(M, y1, y2, iters=20, tolparinv= 1e-06, cores=15)

}
