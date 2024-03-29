% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ERRBLUP_Bivar_Stepwise.R
\name{ERRBLUP_BV_Stepwise}
\alias{ERRBLUP_BV_Stepwise}
\title{Bivariate ERRBLUP Phenotype Prediction Function Relying On The Out Put Of ERRBLUP Relationship Matrix Function}
\usage{
ERRBLUP_BV_Stepwise(y1, y2, G_ERRBLUP, iters = 20, tolparinv = 1e-06)
}
\arguments{
\item{y1}{A numeric vector of named phenotypes for the target triat to be predicted in bivariate ERRBLUP}

\item{y2}{A numeric vector of named phenotypes for the additional trait to be estimated in bivariate ERRBLUP}

\item{G_ERRBLUP}{ERRBLUP relationship matrix with the row names and column names of all individuals}

\item{iters}{Maximum number of iterations allowed with the default value of \code{20} (The parameter of mmer function from sommer package)}

\item{tolparinv}{Tolerance parameter for matrix inverse used when singularities are encountered in the estimation procedure with the default value of \code{1e-06} (The parameter of mmer function from sommer package)}
}
\value{
A dataframe of both phenotype predictions of target triat test set, phenotype estimations of target triat training set and phenotype estimations of additional trait full set in bivarite ERRBLUP method
}
\description{
Function to do phenotype prediction based on all pairwise SNP interactions in bivariate model
}
\examples{
library(BGLR)
data(wheat)
m <- Recodemarkers(wheat.X)
y1 <- Phenotype_BV[,1]
N <- length(y1)
n <- 60
test <- sample(1:N,n)
y1[test] <- NA
y2 <- Phenotype_BV[,2]
G_ERRBLUP <- Gall(m, cores=15)
G <- G_ERRBLUP$G
rownames(G) <- names(y1)
colnames(G) <- names(y1)
ERRBLUP_Bivar <- ERRBLUP_BV_Stepwise(y1, y2, G, iters=20, tolparinv= 1e-06)

}
