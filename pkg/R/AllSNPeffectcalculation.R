
#' @title Pairwise SNP Interaction Effects and variances Estimation Function
#'
#' @description Function to calculate all pairwise SNP interaction effects
#'
#' @param m {0,1,2} or {0,2} coded marker matrix with individuals in the rows and the markers in the columns
#' @param Pheno A numeric vector of phenotypes
#' @param G_ERRBLUP EERRBLUP relationship matrix
#' @param P A vector of all genotype combinations frequencies in the population
#' @param cores The number of cores with the default value of 1
#'
#'@return A list of two components effect and effectvar
#'
#' \describe{
#'   \item{effect}{A numeric vector of all estimated pairwise SNP interaction effects}
#'   \item{effectvar}{A numeric vector of all estimated pairwise SNP interaction effects variances}
#' }
#'
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' N <- length(Phenotype)
#' n <- 60
#' test <- sample(1:N,n)
#' Phenotype[test] <- NA
#' m <- Recodemarkers(wheat.X[,1:10])
#' G_ERRBLUP <- Gall(m, cores=15)
#' G <- G_ERRBLUP$G
#' P <- G_ERRBLUP$P
#' Estimation <- SNP_effect_var(m, Phenotype, G, P, cores=15)
#' t_hat <- Estimation$effect
#' sigma_hat <- Estimation$effectvar
#'
#' @export
#'


SNP_effect_var <- function(m, Pheno, G_ERRBLUP, P, cores=1){


  Y <- data.frame(ID = 1:length(Pheno), observation = Pheno)
  phenosid <- Y[stats::complete.cases(Y[,2]),]

  y <- phenosid[, 2]
  ntrain <- length(y)

  Trainset <- phenosid[,1]

  Gall_train <- G_ERRBLUP[Trainset, Trainset]

  R <- diag(ntrain)
  multi <- y - mean(y)

  abc <- EMMREML::emmreml(y, X=cbind(matrix(rep(1, ntrain), ncol=1)), Z=diag(ntrain), K=Gall_train)
  vare <- abc$Ve
  varg <- abc$Vu
  lambda <- vare/varg

  if(lambda<0.001){
    lambda = 0.001
  }
  if(lambda>1000){
    lambda = 1000
  }

  Rest_term <- (chol2inv(chol(Gall_train + R*lambda)) %*% multi)

  Z <- t(m)
  nsnp <- nrow(Z)
  nindi <- ncol(Z)

  Z0 <- (Z==0)*2L
  Z1 <- (Z==1)*2L
  Z2 <- (Z==2)*2L


  if(sum(Z1==0)== nsnp*nindi){

    u_hat <- rep(0L, nsnp*(nsnp+1)*2)
    Z_share = matrix(0L, ncol=nindi, nrow=(nsnp+1)*4)


    for(index in 1:ceiling(nsnp/2)){

      print(index)

      temp1 = matrix(Z[index,]==0, ncol=nindi, nrow=nsnp-index+1, byrow=TRUE)
      temp2 = Z0[index:nsnp,,drop=FALSE]
      temp3 = Z2[index:nsnp,,drop=FALSE]

      Z_share[1:(nsnp-index+1),] <- temp1 * temp2
      Z_share[1:(nsnp-index+1) + (nsnp+1),] <- (!temp1) * temp2
      Z_share[1:(nsnp-index+1) + 2*(nsnp+1),] <- temp1 * temp3
      Z_share[1:(nsnp-index+1) + 3*(nsnp+1),] <- (!temp1) * temp3

      if(index <= (nsnp/2)){
        temp1 = matrix(Z[(nsnp-index+1),]==0, ncol=nindi, nrow=index, byrow=TRUE)
        temp2 = Z0[(nsnp-index+1):nsnp,,drop=FALSE]
        temp3 = Z2[(nsnp-index+1):nsnp,,drop=FALSE]
      } else{
        temp1 = matrix(0L, ncol=nindi, nrow=index, byrow=TRUE)
      }

      Z_share[(nsnp-index+2):(nsnp+1),] <- temp1 * temp2
      Z_share[(nsnp-index+2):(nsnp+1) + (nsnp+1),] <- (!temp1) * temp2
      Z_share[(nsnp-index+2):(nsnp+1) + 2*(nsnp+1),] <- temp1 * temp3
      Z_share[(nsnp-index+2):(nsnp+1) + 3*(nsnp+1),] <- (!temp1) * temp3

      if (requireNamespace("miraculix", quietly = TRUE)) {

        RandomFieldsUtils::RFoptions(cores=cores)


        Z_new <- Z_share[c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1)),]
        Z_miraculix <- miraculix::genomicmatrix(Z_new[,Trainset])
        u_hat[1:(4*(nsnp-index+1)) + (index-1)*4*nsnp - 2 * (index-1) * (index-2)] <- miraculix::genoVector(Z_miraculix, Rest_term)

        if(index<= (nsnp/2)){
          index2 = nsnp - index + 1
          Z_new <- Z_share[-c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1)),]
          Z_miraculix <- miraculix::genomicmatrix(Z_new[,Trainset])
          u_hat[1:(4*(nsnp-index2+1)) + (index2-1)*4*nsnp - 2 * (index2-1) * (index2-2)] <- miraculix::genoVector(Z_miraculix, Rest_term)
        }

      } else{

        Z_new <- Z_share[c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1)),]
        A <- Z_new[,Trainset]
        u_hat[1:(4*(nsnp-index+1)) + (index-1)*4*nsnp - 2 * (index-1) * (index-2)] <- A %*% Rest_term

        if(index<= (nsnp/2)){
          index2 = nsnp - index + 1
          Z_new <- Z_share[-c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1)),]
          A <- Z_new[,Trainset]
          u_hat[1:(4*(nsnp-index2+1)) + (index2-1)*4*nsnp - 2 * (index2-1) * (index2-2)] <- A %*% Rest_term
        }
      }
    }

  } else {

    u_hat <- rep(NA, nsnp*(nsnp+1)*9/2)
    Z_share = matrix(0L, ncol=nindi, nrow=(nsnp+1)*9)

    for(index in 1:ceiling(nsnp/2)){

      print(index)

      temp0 = matrix(Z[index,]==0, ncol=nindi, nrow=nsnp-index+1, byrow=TRUE)
      temp1 = matrix(Z[index,]==1, ncol=nindi, nrow=nsnp-index+1, byrow=TRUE)
      temp2 = Z0[index:nsnp,,drop=FALSE]
      temp3 = Z2[index:nsnp,,drop=FALSE]
      temp4 = Z1[index:nsnp,,drop=FALSE]

      Z_share[1:(nsnp-index+1),] <- temp0 * temp2
      Z_share[1:(nsnp-index+1) + (nsnp+1),] <- ((!temp0)&(!temp1)) * temp2
      Z_share[1:(nsnp-index+1) + 2*(nsnp+1),] <- temp0 * temp3
      Z_share[1:(nsnp-index+1) + 3*(nsnp+1),] <- ((!temp0)&(!temp1)) * temp3
      Z_share[1:(nsnp-index+1) + 4*(nsnp+1),] <- temp1 * temp3
      Z_share[1:(nsnp-index+1) + 5*(nsnp+1),] <- ((!temp0)&(!temp1)) * temp4
      Z_share[1:(nsnp-index+1) + 6*(nsnp+1),] <- temp1 * temp2
      Z_share[1:(nsnp-index+1) + 7*(nsnp+1),] <- temp0 * temp4
      Z_share[1:(nsnp-index+1) + 8*(nsnp+1),] <- temp1 * temp4


      if(index <= (nsnp/2)){
        temp0 = matrix(Z[(nsnp-index+1),]==0, ncol=nindi, nrow=index, byrow=TRUE)
        temp1 = matrix(Z[(nsnp-index+1),]==1, ncol=nindi, nrow=index, byrow=TRUE)
        temp2 = Z0[(nsnp-index+1):nsnp,,drop=FALSE]
        temp3 = Z2[(nsnp-index+1):nsnp,,drop=FALSE]
        temp4 = Z1[(nsnp-index+1):nsnp,,drop=FALSE]

      } else{
        temp0 = matrix(0L, ncol=nindi, nrow=index, byrow=TRUE)
        temp1 = matrix(0L, ncol=nindi, nrow=index, byrow=TRUE)
      }

      Z_share[(nsnp-index+2):(nsnp+1),] <- temp0 * temp2
      Z_share[(nsnp-index+2):(nsnp+1) + (nsnp+1),] <- ((!temp0)&(!temp1)) * temp2
      Z_share[(nsnp-index+2):(nsnp+1) + 2*(nsnp+1),] <- temp0 * temp3
      Z_share[(nsnp-index+2):(nsnp+1) + 3*(nsnp+1),] <- ((!temp0)&(!temp1)) * temp3
      Z_share[(nsnp-index+2):(nsnp+1) + 4*(nsnp+1),] <- temp1 * temp3
      Z_share[(nsnp-index+2):(nsnp+1) + 5*(nsnp+1),] <- ((!temp0)&(!temp1)) * temp4
      Z_share[(nsnp-index+2):(nsnp+1) + 6*(nsnp+1),] <- temp1 * temp2
      Z_share[(nsnp-index+2):(nsnp+1) + 7*(nsnp+1),] <- temp0 * temp4
      Z_share[(nsnp-index+2):(nsnp+1) + 8*(nsnp+1),] <- temp1 * temp4

      if (requireNamespace("miraculix", quietly = TRUE)) {


        Z_new <- Z_share[c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1),
                           1:(nsnp-index+1) + 4*(nsnp+1), 1:(nsnp-index+1) + 5*(nsnp+1), 1:(nsnp-index+1) + 6*(nsnp+1),
                           1:(nsnp-index+1) + 7*(nsnp+1), 1:(nsnp-index+1) + 8*(nsnp+1)),]
        Z_miraculix <- miraculix::genomicmatrix(Z_new[,Trainset])
        u_hat[1:(9*(nsnp-index+1)) + (index-1)*9*nsnp - 2 * (index-1) * (index-2)] <- miraculix::genoVector(Z_miraculix, Rest_term)

        if(index<= (nsnp/2)){
          index2 = nsnp - index + 1
          Z_new <- Z_share[-c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1),
                              1:(nsnp-index+1) + 4*(nsnp+1), 1:(nsnp-index+1) + 5*(nsnp+1), 1:(nsnp-index+1) + 6*(nsnp+1),
                              1:(nsnp-index+1) + 7*(nsnp+1), 1:(nsnp-index+1) + 8*(nsnp+1)),]
          Z_miraculix <- miraculix::genomicmatrix(Z_new[,Trainset])
          u_hat[1:(9*(nsnp-index2+1)) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2)] <- miraculix::genoVector(Z_miraculix, Rest_term)
        }


      } else{

        Z_new <- Z_share[c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1),
                           1:(nsnp-index+1) + 4*(nsnp+1), 1:(nsnp-index+1) + 5*(nsnp+1), 1:(nsnp-index+1) + 6*(nsnp+1),
                           1:(nsnp-index+1) + 7*(nsnp+1), 1:(nsnp-index+1) + 8*(nsnp+1)),]
        A <- Z_new[,Trainset]
        u_hat[1:(9*(nsnp-index+1)) + (index-1)*9*nsnp - 2 * (index-1) * (index-2)] <- A %*% Rest_term

        if(index<= (nsnp/2)){
          index2 = nsnp - index + 1
          Z_new <- Z_share[-c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1),
                              1:(nsnp-index+1) + 4*(nsnp+1), 1:(nsnp-index+1) + 5*(nsnp+1), 1:(nsnp-index+1) + 6*(nsnp+1),
                              1:(nsnp-index+1) + 7*(nsnp+1), 1:(nsnp-index+1) + 8*(nsnp+1)),]
          A <- Z_new[,Trainset]
          u_hat[1:(9*(nsnp-index2+1)) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2)] <- A %*% Rest_term
        }
      }
    }
  }

  u_hat <- u_hat  * 1/ 2 / sum(P*(1-P))
  sigma_hat <- (u_hat^2)*2*P*(1-P)

  out <- list(effect = u_hat, effectvar = sigma_hat)
  return(out)


}





