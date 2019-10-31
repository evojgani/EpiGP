
#' @title Pairwise SNP Interaction Effects Variances Function
#'
#' @description Function to calculate all pairwise SNP interaction effects varinacs
#'
#' @param m {0,1,2} or {0,2} coding Marker matrix with individuals in the rows and the markers in the columns
#' @param pheno_train A subset of one phenotype vector as a training set with names for each phenotypic value
#' @param t_hat A vector of all estimated pairwise SNP interaction effects
#'
#' @return A vector of all pairwise SNP interaction effects variances
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' pheno <- wheat.Y[1:100,1]
#' pheno_train <- pheno[1:round(4*length(pheno)/5)]
#' m <- Recodemarker(wheat.X[1:100,])
#' rownames(m) <- names(pheno)
#' G_ERRBLUP <- Gall(m)
#' t_hat <- SNP_effect(m, pheno_train, G_ERRBLUP)
#' sigma_hat <- SNP_var(m, pheno_train, t_hat)
#'
#' @export
#'


SNP_var <- function(m, pheno_train, t_hat){

  m <- m[rownames(m) %in% names(pheno_train), ] # names(y_real)=Genotype
  Z <- t(m)


  nsnp <- nrow(Z)
  nindi <- ncol(Z)


  # this is cheating since i am assuming the heritability to be known but will not cost a lot of computation time

  Z0 <- (Z==0)*2L
  Z1 <- (Z==1)*2L
  Z2 <- (Z==2)*2L


  if(sum(Z1==0)== nsnp*nindi){

    sigma_hat <- numeric(nsnp*nsnp*4)
    p_i <- numeric(nsnp*nsnp*4)
    include <- integer(nsnp*nsnp*4)+1L

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
      sigma_hat[((index-1)*nsnp*4+1):((index)*nsnp*4)] <- (t_hat[((index-1)*nsnp*4+1):((index)*nsnp*4)])^2*2*p_i[((index-1)*nsnp*4+1):((index)*nsnp*4)]*(1-p_i[((index-1)*nsnp*4+1):((index)*nsnp*4)])

    }

  } else {

    sigma_hat <- numeric(nsnp*nsnp*9)
    p_i <- numeric(nsnp*nsnp*9)
    include <- integer(nsnp*nsnp*9)+1L

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
      sigma_hat[((index-1)*nsnp*4+1):((index)*nsnp*4)] <- (t_hat[((index-1)*nsnp*4+1):((index)*nsnp*4)])^2*2*p_i[((index-1)*nsnp*4+1):((index)*nsnp*4)]*(1-p_i[((index-1)*nsnp*4+1):((index)*nsnp*4)])

    }
  }

  return(sigma_hat)

}




