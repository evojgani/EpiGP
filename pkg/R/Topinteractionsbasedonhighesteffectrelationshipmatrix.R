
#' @title sERRBLUP Relationship Matrix Function
#'
#' @description Function to generate relationship matrix based on the desired proportion of pairwise SNP interactions
#'
#' @param m {0,1,2} or {0,2} coded marker matrix with individuals in the rows and the markers in the columns
#' @param pheno A numeric phenotype vector with names for each phenotypic value
#' @param estimations A numeric vector of all estimated pairwise SNP interaction effects or all estimated pairwise SNP interaction effects variances
#' @param k Desired proportion of SNP interactions to be included in the model
#'
#' @return sERRBLUP Relationship matrix for the k percent of pairwise SNP interactions with row names and column names of all the individuals
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' pheno <- wheat.Y[1:100,1]
#' pheno_train <- pheno[1:round(4*length(pheno)/5)]
#' m <- Recodemarkers(wheat.X[1:100,])
#' rownames(m) <- names(pheno)
#' G_ERRBLUP <- Gall(m)
#' t_hat <- SNP_effect(m, pheno_train, G_ERRBLUP)
#' sigma_hat <- SNP_var(m, pheno_train, t_hat)
#' k <- 10
#' Gtop_effect <- Gtop(m, pheno , t_hat, k)
#' Gtop_var <- Gtop(m, pheno , sigma_hat, k)
#'
#' @export
#'


Gtop <- function(m, pheno, estimations, k){

  if(is.null(row.names(m))|is.null(names(pheno))){

    stop("The individuals are not named")

  } else {

  m <- m[rownames(m) %in% names(pheno), ]
  Z <- t(m)


  nsnp <- nrow(Z)
  nindi <- ncol(Z)


    G <- matrix(0, ncol=nindi, nrow=nindi)

  storage.mode(Z) = "integer"
  attr(Z, "dimnames") = NULL

  Z0 <- (Z==0)*2L
  Z1 <- (Z==1)*2L
  Z2 <- (Z==2)*2L


  if(sum(Z1==0)== nsnp*nindi){

    p_i <- numeric(nsnp*nsnp*4)

    include <- integer(nsnp*nsnp*4)+1L
    include[abs(estimations)< stats::quantile(abs(estimations),(1-(k/100)))] <- 0L

    rm(estimations)

    Z_share = matrix(0L, ncol=nindi, nrow=nsnp*4)
    check = prod(include)

    for(index in 1:nsnp){
      if(index %% 1000 == 0)print(index)

      Z_share[1:nsnp,] <- matrix(Z[index,]==0, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z0
      Z_share[1:nsnp+nsnp,] <- matrix(Z[index,]==2, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z0
      Z_share[1:nsnp+2*nsnp,] <- matrix(Z[index,]==0, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z2
      Z_share[1:nsnp+3*nsnp,] <- matrix(Z[index,]==2, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z2
      if(check!=1){
        Z_share <- matrix(include[((index-1)*nsnp*4+1):((index)*nsnp*4)], ncol=nindi, nrow=nsnp*4, byrow=TRUE) * Z_share
      }

      p_i[((index-1)*nsnp*4+1):((index)*nsnp*4)] <- rowSums(Z_share)/ncol(Z_share)/2

      if (requireNamespace("miraculix", quietly = TRUE)) {

        Z_miraculix <- miraculix::genomicmatrix(Z_share)
        G <- G +  miraculix::relationshipMatrix(Z_miraculix, centered=TRUE, normalized=FALSE)

      } else{

        A <- Z_share - 2*p_i[((index-1)*nsnp*4+1):((index)*nsnp*4)]
        G <- G + crossprod(A)
      }
    }
  } else {

    p_i <- numeric(nsnp*nsnp*9)

    include <- integer(nsnp*nsnp*9)+1L
    include[abs(estimations)< stats::quantile(abs(estimations),(1-(k/100)))] <- 0L

    rm(estimations)

    Z_share = matrix(0L, ncol=nindi, nrow=nsnp*9)
    check = prod(include)

    for(index in 1:nsnp){
      if(index %% 1000 == 0)print(index)


      Z_share[1:nsnp,] <- matrix(Z[index,]==0, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z0
      Z_share[1:nsnp+nsnp,] <- matrix(Z[index,]==2, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z0
      Z_share[1:nsnp+2*nsnp,] <- matrix(Z[index,]==0, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z2
      Z_share[1:nsnp+3*nsnp,] <- matrix(Z[index,]==2, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z2
      Z_share[1:nsnp+4*nsnp,] <- matrix(Z[index,]==0, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z1
      Z_share[1:nsnp+5*nsnp,] <- matrix(Z[index,]==1, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z0
      Z_share[1:nsnp+6*nsnp,] <- matrix(Z[index,]==1, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z1
      Z_share[1:nsnp+7*nsnp,] <- matrix(Z[index,]==1, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z2
      Z_share[1:nsnp+8*nsnp,] <- matrix(Z[index,]==2, ncol=nindi, nrow=nsnp, byrow=TRUE) * Z1

      if(check!=1){
        Z_share <- matrix(include[((index-1)*nsnp*9+1):((index)*nsnp*9)], ncol=nindi, nrow=nsnp*9, byrow=TRUE) * Z_share
      }

      p_i[((index-1)*nsnp*9+1):((index)*nsnp*9)] <- rowSums(Z_share)/ncol(Z_share)/2

      if (requireNamespace("miraculix", quietly = TRUE)) {

        Z_miraculix <- miraculix::genomicmatrix(Z_share)
        G <- G +  miraculix::relationshipMatrix(Z_miraculix, centered=TRUE, normalized=FALSE)

      } else{

        A <- Z_share - 2*p_i[((index-1)*nsnp*4+1):((index)*nsnp*4)]
        G <- G + crossprod(A)
      }
    }
  }



  G_k <- G / (2 * sum(p_i*(1-p_i)))

  rownames(G_k) <- rownames(m)
  colnames(G_k) <- rownames(m)

  return(G_k)

}

}






