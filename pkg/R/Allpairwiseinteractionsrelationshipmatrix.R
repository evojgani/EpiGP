
#' @title ERRBLUP Relationship Matrix Function
#'
#' @description Function to generate relationship matrix based on all pairwise SNP interactions
#'
#' @param m The \code{{0,1,2}} or \code{{0,2}} coded marker matrix with individuals in the rows and the markers in the columns
#' @param cores The number of cores with the default value of \code{1}
#'
#'@return A list of two components G and P
#'
#' \describe{
#'   \item{G}{ERRBLUP relationship matrix}
#'   \item{P}{A vector of all genotype combinations frequencies in the population}
#' }
#'
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' m <- Recodemarkers(wheat.X)
#' G_ERRBLUP <- Gall(m, cores=15)
#' G <- G_ERRBLUP$G
#' P <- G_ERRBLUP$P
#'
#' @export
#'

Gall <- function(m, cores=1){

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

       p_i <- rep(NA, nsnp*(nsnp+1)*2)
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

         Z_miraculix <- miraculix::genomicmatrix(Z_share)

         pi1 <- miraculix::allele_freq(Z_miraculix)


         p_i[1:(4*(nsnp-index+1)) + (index-1)*4*nsnp - 2 * (index-1) * (index-2)] <-
           pi1[c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1))]

         if(index<= (nsnp/2)){
           index2 = nsnp - index + 1
           p_i[1:(4*(nsnp-index2+1)) + (index2-1)*4*nsnp - 2 * (index2-1) * (index2-2)] <-
             pi1[-c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1))]
         }

         if(index > (nsnp/2)){
           Z_share <- Z_share[-(c((nsnp-index+2):(nsnp+1),
                                  (nsnp-index+2):(nsnp+1) + 1*(nsnp+1),
                                  (nsnp-index+2):(nsnp+1) + 2*(nsnp+1),
                                  (nsnp-index+2):(nsnp+1) + 3*(nsnp+1))),]
           Z_miraculix <- miraculix::genomicmatrix(Z_share)
         }

         G <- G +  miraculix::relationshipMatrix(Z_miraculix, centered=FALSE, normalized=FALSE)

       } else{

         pi1 <- rowSums(Z_share)/ncol(Z_share)/2

         p_i[1:(4*(nsnp-index+1)) + (index-1)*4*nsnp - 2 * (index-1) * (index-2)] <-
           pi1[c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1))]

         if(index<= (nsnp/2)){
           index2 = nsnp - index + 1
           p_i[1:(4*(nsnp-index2+1)) + (index2-1)*4*nsnp - 2 * (index2-1) * (index2-2)] <-
             pi1[-c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1))]
         }

          if(index > (nsnp/2)){
            Z_share <- Z_share[-(c((nsnp-index+2):(nsnp+1),
                                   (nsnp-index+2):(nsnp+1) + 1*(nsnp+1),
                                   (nsnp-index+2):(nsnp+1) + 2*(nsnp+1),
                                   (nsnp-index+2):(nsnp+1) + 3*(nsnp+1))),]
          }

          G <- G + crossprod(Z_share)
       }
      }
    } else {

    p_i <- rep(NA, nsnp*(nsnp+1)*9/2)
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

        RandomFieldsUtils::RFoptions(cores=cores)

        Z_miraculix <- miraculix::genomicmatrix(Z_share)

        pi1 <- miraculix::allele_freq(Z_miraculix)


        p_i[1:(9*(nsnp-index+1)) + (index-1)*9*nsnp - 2 * (index-1) * (index-2)] <-
          pi1[c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1),
                1:(nsnp-index+1) + 4*(nsnp+1), 1:(nsnp-index+1) + 5*(nsnp+1), 1:(nsnp-index+1) + 6*(nsnp+1),
                1:(nsnp-index+1) + 7*(nsnp+1), 1:(nsnp-index+1) + 8*(nsnp+1))]

        if(index<= (nsnp/2)){
          index2 = nsnp - index + 1
          p_i[1:(9*(nsnp-index2+1)) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2)] <-
            pi1[-c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1),
                   1:(nsnp-index+1) + 4*(nsnp+1), 1:(nsnp-index+1) + 5*(nsnp+1), 1:(nsnp-index+1) + 6*(nsnp+1),
                   1:(nsnp-index+1) + 7*(nsnp+1), 1:(nsnp-index+1) + 8*(nsnp+1))]
        }

        if(index > (nsnp/2)){
          Z_share <- Z_share[-(c((nsnp-index+2):(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 1*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 2*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 3*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 4*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 5*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 6*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 7*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 8*(nsnp+1))),]
          Z_miraculix <- miraculix::genomicmatrix(Z_share)
        }

        G <- G +  miraculix::relationshipMatrix(Z_miraculix, centered=FALSE, normalized=FALSE)

      } else{

        pi1 <- rowSums(Z_share)/ncol(Z_share)/2

        p_i[1:(9*(nsnp-index+1)) + (index-1)*9*nsnp - 2 * (index-1) * (index-2)] <-
          pi1[c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1),
                1:(nsnp-index+1) + 4*(nsnp+1), 1:(nsnp-index+1) + 5*(nsnp+1), 1:(nsnp-index+1) + 6*(nsnp+1),
                1:(nsnp-index+1) + 7*(nsnp+1), 1:(nsnp-index+1) + 8*(nsnp+1))]

        if(index<= (nsnp/2)){
          index2 = nsnp - index + 1
          p_i[1:(9*(nsnp-index2+1)) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2)] <-
            pi1[-c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp+1), 1:(nsnp-index+1)+2*(nsnp+1), 1:(nsnp-index+1) + 3*(nsnp+1),
                   1:(nsnp-index+1) + 4*(nsnp+1), 1:(nsnp-index+1) + 5*(nsnp+1), 1:(nsnp-index+1) + 6*(nsnp+1),
                   1:(nsnp-index+1) + 7*(nsnp+1), 1:(nsnp-index+1) + 8*(nsnp+1))]
        }

        if(index > (nsnp/2)){
          Z_share <- Z_share[-(c((nsnp-index+2):(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 1*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 2*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 3*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 4*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 5*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 6*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 7*(nsnp+1),
                                 (nsnp-index+2):(nsnp+1) + 8*(nsnp+1))),]
        }

        G <- G + crossprod(Z_share)
      }
     }
    }


  G_all <- G / (2 * sum(p_i*(1-p_i)))


  out <- list(G = G_all, P = p_i)

  return(out)

}




