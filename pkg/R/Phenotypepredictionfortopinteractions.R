
#' @title sERRBLUP Phenotype Prediction Function
#'
#' @description Function to do phenotype prediction based on the desired proportion of pairwise SNP interactions
#'
#' @param pheno_train A subset of one numeric phenotype vector as a training set with names for each phenotypic value
#' @param G_ERRBLUP ERRBLUP relationship matrix with row names and column names of all the individuals
#' @param Gtop sERRBLUP Relationship matrix for the k percent of pairwise SNP interactions with row names and column names of all the individuals
#'
#' @return A numeric vector of both phenotype estimations of training set and phenotype predictions of test set based on sERRBLUP method
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' pheno <- wheat.Y[,1]
#' pheno_train <- pheno[1:round(4*length(pheno)/5)]
#' m <- Recodemarkers(wheat.X)
#' rownames(m) <- names(pheno)
#' G_ERRBLUP <- Gall(m)
#' t_hat <- SNP_effect(m, pheno_train, G_ERRBLUP)
#' sigma_hat <- SNP_var(m, pheno_train, t_hat)
#' k <- 10
#' Gtop_effect <- Gtop(m, pheno, t_hat, k)
#' Gtop_var <- Gtop(m, pheno, sigma_hat, k)
#' sERRBLUP_effect <- sERRBLUP(pheno_train, G_ERRBLUP, Gtop_effect)
#' sERRBLUP_var <- sERRBLUP(pheno_train, G_ERRBLUP, Gtop_var)
#'
#' @export
#'


sERRBLUP <- function(pheno_train, G_ERRBLUP, Gtop) {

  if(is.null(names(pheno_train))){

    stop("The individuals are not named")

  } else {

  Pheno <- pheno_train[stats::complete.cases(pheno_train)]
  Phenosid <- data.frame(ID = names(pheno_train), observation = Pheno)

  n <- dim(Gtop)[1]
  Zz <- diag(n)
  Xx <- matrix(1,n,ncol=1)

  y <- Phenosid[,2]
  ntrain <- length(y)
  X <- Xx[1:ntrain]
  Z <- Zz[1:ntrain,]

  Ginv <- MASS::ginv(Gtop)

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
}



