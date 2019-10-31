
#' @title ERRBLUP Phenotype Prediction Function
#'
#' @description Function to do phenotype prediction based on all pairwise SNP interactions
#'
#' @param pheno_train A subset of one phenotype vector as a training set with names for each phenotypic value
#' @param G_ERRBLUP ERRBLUP relationship matrix with row names and column names of all the individuals
#'
#' @return A vector of both phenotype estimations of training set and phenotype predictions of test set based on ERRBLUP method
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' pheno <- wheat.Y[,1]
#' m <- Recodemarker(wheat.X)
#' rownames(m) <- names(pheno)
#' G_ERRBLUP <- Gall(m)
#' pheno_train <- pheno[1:round(4*length(pheno)/5)]
#' Pred_ERRBLUP <- ERRBLUP(pheno_train, G_ERRBLUP)
#'
#' @export
#'


ERRBLUP <- function(pheno_train , G_ERRBLUP) {



  Pheno <- pheno_train[stats::complete.cases(pheno_train)]
  Phenosid <- data.frame(ID = names(Pheno), observation = Pheno)

  n <- dim(G_ERRBLUP)[1]
  Zz <- diag(n)
  Xx <- matrix(1,n,ncol=1)

  y <- Phenosid[,2]
  ntrain <- length(y)
  X <- Xx[1:ntrain]
  Z <- Zz[1:ntrain,]

  Ginv <- MASS::ginv(G_ERRBLUP)

  Gtrain <- G_ERRBLUP[rownames(G_ERRBLUP) %in% Phenosid[,1],colnames(G_ERRBLUP) %in% Phenosid[,1]]

  abc <- EMMREML::emmreml(y=y, X=cbind(matrix(1, nrow = ntrain, ncol=1)), Z=diag(ntrain), K=Gtrain)
  vare <- abc$Ve
  varg <- abc$Vu
  lambda <- vare/varg

  LHS <- rbind(cbind(crossprod(X),crossprod(X,Z)),cbind(crossprod(Z,X),crossprod(Z)+Ginv*lambda))
  RHS <- rbind(crossprod(X,y),crossprod(Z,y))
  MME <- MASS::ginv(LHS)%*%RHS
  prediction <- MME[-1,]+MME[1,1]

  return(prediction)

}





