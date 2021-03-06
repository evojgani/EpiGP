
#' @title sERRBLUP Phenotype Prediction Function Relying On The Out Put Of sERRBLUP Relationship Matrix Function
#'
#' @description Function to do phenotype prediction based on the desired proportion of pairwise SNP interactions
#'
#' @param Pheno A numeric vector of phenotypes
#' @param Gtop sERRBLUP Relationship matrix for the k percent of pairwise SNP interactions
#'
#' @return A numeric vector of both phenotype estimations of training set and phenotype predictions of test set based on sERRBLUP method
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' N <- length(Phenotype)
#' n <- 60
#' test <- sample(1:N,n)
#' pheno <- Phenotype
#' pheno[test] <- NA
#' m <- Recodemarkers(wheat.X)
#' G_ERRBLUP <- Gall(m, cores=15)
#' G <- G_ERRBLUP$G
#' P <- G_ERRBLUP$P
#' Estimation <- SNP_Effect_Var(m, pheno, G, P, cores=15)
#' t_hat <- Estimation$Effect
#' sigma_hat <- Estimation$Effect.Var
#' k <- 10
#' Gtop_effect <- Gtop(m, t_hat, k, cores=15)
#' Gtop_var <- Gtop(m, sigma_hat, k, cores=15)
#' sERRBLUP_effect <- sERRBLUP_Stepwise(pheno, Gtop_effect)
#' sERRBLUP_var <- sERRBLUP_Stepwise(pheno, Gtop_var)
#'
#' @export
#'


sERRBLUP_Stepwise <- function(Pheno, Gtop) {


    Y <- data.frame(ID = 1:length(Pheno), observation = Pheno)
    phenosid <- Y[stats::complete.cases(Y[,2]),]
    Trainset <- phenosid[,1]


    n <- dim(Gtop)[1]
    Zz <- diag(n)
    Xx <- matrix(1,n,ncol=1)

    y <- phenosid[,2]
    ntrain <- length(y)
    X <- Xx[Trainset]
    Z <- Zz[Trainset,]

    Ginv <- MASS::ginv(Gtop)

    Gtrain <- Gtop[Trainset,Trainset]

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





