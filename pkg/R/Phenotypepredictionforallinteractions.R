
#' @title ERRBLUP Phenotype Prediction Function
#'
#' @description Function to do phenotype prediction based on all pairwise SNP interactions
#'
#' @param Pheno_train A subset of one numeric phenotype vector as a training set with names for each phenotypic value
#' @param G_ERRBLUP ERRBLUP relationship matrix with row names and column names of all the individuals
#' @param Trainset A vector of individuals which are in the training set
#'
#' @return A numeric vector of both phenotype estimations of training set and phenotype predictions of test set based on ERRBLUP method
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' geno <- wheat.X
#' t1 <- sample(1:ncol(geno), 20)
#' t2 <- sample(1:ncol(geno), 20)
#' y1 <- rowSums((geno[,t1]==2) * (geno[,t2]==2))
#' t1 <- sample(1:ncol(geno), 20)
#' t2 <- sample(1:ncol(geno), 20)
#' y2 <- rowSums((geno[,t1]==2) * (geno[,t2]==0))
#' t1 <- sample(1:ncol(geno), 20)
#' t2 <- sample(1:ncol(geno), 20)
#' y3 <- rowSums((geno[,t1]==0) * (geno[,t2]==2))
#' t1 <- sample(1:ncol(geno), 20)
#' t2 <- sample(1:ncol(geno), 20)
#' y4 <- rowSums((geno[,t1]==0) * (geno[,t2]==0))
#' y <- y1+y2+y3+y4
#' pheno <- scale(y)
#' names(pheno) <- names(wheat.Y[,1])
#' m <- Recodemarkers(wheat.X)
#' rownames(m) <- names(pheno)
#' G_ERRBLUP <- Gall(m, cores=15)
#' G <- G_ERRBLUP$G
#' N <- length(pheno)
#' n <- 60
#' test <- sample(1:N,n)
#' training <- (1:N)[-test]
#' pheno_train <- pheno[training]
#' Pred_ERRBLUP <- ERRBLUP(pheno_train, G, training)
#'
#' @export
#'


ERRBLUP <- function(Pheno_train , G_ERRBLUP, Trainset) {

  if(is.null(names(Pheno_train))){

    stop("The individuals are not named")

  } else {

    Pheno <- Pheno_train[stats::complete.cases(Pheno_train)]
    Phenosid <- data.frame(ID = names(Pheno), observation = Pheno)

    n <- dim(G_ERRBLUP)[1]
    Zz <- diag(n)
    Xx <- matrix(1,n,ncol=1)

    y <- Phenosid[,2]
    ntrain <- length(y)
    X <- Xx[Trainset]
    Z <- Zz[Trainset,]

    Ginv <- MASS::ginv(G_ERRBLUP)

    Gtrain <- G_ERRBLUP[Trainset,Trainset]

    abc <- EMMREML::emmreml(y=y, X=cbind(matrix(1, nrow = ntrain, ncol=1)), Z=diag(ntrain), K=Gtrain)
    vare <- abc$Ve
    varg <- abc$Vu
    lambda <- vare/varg

    if(lambda<0.001){
      lambda = 0.001
    }
    if(lambda>1000){
      lambda = 1000
    }

    LHS <- rbind(cbind(crossprod(X),crossprod(X,Z)),cbind(crossprod(Z,X),crossprod(Z)+Ginv*lambda))
    RHS <- rbind(crossprod(X,y),crossprod(Z,y))
    MME <- MASS::ginv(LHS)%*%RHS
    prediction <- MME[-1,]+MME[1,1]


    return(prediction)

  }

}




