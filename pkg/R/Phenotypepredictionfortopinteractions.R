
#' @title sERRBLUP Phenotype Prediction Function
#'
#' @description Function to do phenotype prediction based on the top desired pairwise SNP interaction
#'
#' @param pheno_train A subset of one phenotype vector as a training set with names for each phenotypic value
#' @param Gall ERRBLUP relationship matrix with row names and column names of all the individuals
#' @param Gtop sERRBLUP Relationship matrix for the top k percent of pairwise SNP interactions with row names and column names of all the individuals
#'
#' @return A vector of both phenotype estimations of training set and phenotype predictions of test set based on sERRBLUP method
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' pheno <- wheat.Y[,1]
#' pheno_train <- pheno[1:round(4*length(pheno)/5)]
#' M <- Recodemarker(wheat.X)
#' rownames(M) <- names(pheno)
#' G_all <- Gall(M)
#' u_hat <- SNP_effect(M, pheno_train, G_all)
#' u_hat_var <- SNP_var(M, pheno_train, u_hat)
#' k <- 10
#' Gtop_effect <- Gtop(M, pheno , u_hat, k)
#' Gtop_var <- Gtop(M, pheno , u_hat_var, k)
#' sERRBLUP_effect <- sERRBLUP(pheno_train, G_all, Gtop_effect)
#' sERRBLUP_var <- sERRBLUP(pheno_train, G_all, Gtop_var)
#'
#' @export
#'


sERRBLUP <- function(pheno_train, Gall, Gtop) {


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

    Gtrain <- Gall[rownames(Gall) %in% Phenosid[,1],colnames(Gall) %in% Phenosid[,1]]

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




