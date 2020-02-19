
#' @title ERRBLUP Phenotype Prediction Function Relying On The Out Put Of ERRBLUP Relationship Matrix Function
#'
#' @description Function to do phenotype prediction based on all pairwise SNP interactions
#'
#' @param Pheno A numeric vector of phenotypes
#' @param G_ERRBLUP ERRBLUP relationship matrix
#'
#' @return A numeric vector of both phenotype estimations of training set and phenotype predictions of test set based on ERRBLUP method
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' m <- Recodemarkers(wheat.X)
#' G_ERRBLUP <- Gall(m, cores=15)
#' G <- G_ERRBLUP$G
#' N <- length(Phenotype)
#' n <- 60
#' test <- sample(1:N,n)
#' Phenotype[test] <- NA
#' Pred_ERRBLUP <- ERRBLUP_Stepwise(Phenotype, G)
#'
#' @export
#'


ERRBLUP_Stepwise <- function(Pheno , G_ERRBLUP) {


    Y <- data.frame(ID = 1:length(Pheno), observation = Pheno)
    phenosid <- Y[stats::complete.cases(Y[,2]),]
    Trainset <- phenosid[,1]

    n <- dim(G_ERRBLUP)[1]
    Zz <- diag(n)
    Xx <- matrix(1,n,ncol=1)

    y <- phenosid[,2]
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




