
#' @title Bivariate ERRBLUP Phenotype Prediction Function
#'
#' @description Function to do phenotype prediction based on all pairwise SNP interactions in bivariate model
#'
#' @param M The original marker matrix of \code{{-1,0,1}}, \code{{0,1}}, \code{{0,1,2}}, \code{{0,2}} or character coded markers with named individuals in the rows and the markers in the columns
#' @param y1 A numeric vector of named phenotypes for the target triat to be predicted in bivariate ERRBLUP
#' @param y2 A numeric vector of named phenotypes for the additional trait to be estimated in bivariate ERRBLUP
#' @param iters Maximum number of iterations allowed with the default value of \code{20} (The parameter of mmer function from sommer package)
#' @param tolparinv Tolerance parameter for matrix inverse used when singularities are encountered in the estimation procedure with the default value of \code{1e-06} (The parameter of mmer function from sommer package)
#' @param cores The number of cores with the default value of \code{1}
#'
#'@return A list of three components:
#'
#' \describe{
#'   \item{Recodedmarkers}{A \code{{0,1,2}} or \code{{0,2}} coded marker matrix}
#'   \item{Relationshipmatrix}{A list of two components: ERRBLUP relationship matrix (G) and a vector of all genotype combinations frequencies in the population (P)}
#'   \item{Predictions}{A dataframe of both phenotype predictions of target triat test set, phenotype estimations of target triat training set and phenotype estimations of additional trait full set in bivarite ERRBLUP method}
#' }
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' y1 <- Phenotype_BV[,1]
#' N <- length(y1)
#' n <- 60
#' test <- sample(1:N,n)
#' y1[test] <- NA
#' y2 <- Phenotype_BV[,2]
#' M <- wheat.X
#' rownames(M) <- names(y1)
#' ERRBLUP_Bivar <- ERRBLUP_BV(M, y1, y2, iters=20, tolparinv= 1e-06, cores=15)
#'
#' @export
#'


ERRBLUP_BV <- function(M, y1, y2, iters=20, tolparinv= 1e-06, cores=1){

  Recodemarkers <- function(M){

    if(sum(c(0,1,2) %in% M)==3 | sum(c(0,2) %in% M)==2){

      stop("The marker matrix has proper coding for EpiGP package")

    } else{

      m <- M

      if(is.character(m)){

        m <- nchar(m, type = "chars")

        if(sum(m==2)>1){

          second <- (M[m==1][1]) != (M[m==1])
          code <- c(print(c((M[m==1][1]),("is recoded to 0"))),
                    print(c((M[m==1][which(second==TRUE)[1]] ),("is recoded to 1"))),
                    print(c((M[m==2][1]),(("is recoded to 2")))))

          a <- M[m==1][1]
          m[M==a] <- 0


        } else{

          second <- (M[m==1][1]) != (M[m==1])
          code <- c(print(c((M[m==1][1]),("is recoded to 0"))),
                    print(c((M[m==1][which(second==TRUE)[1]] ),("is recoded to 2"))))

          a <- M[m==1][1]
          m[M==a] <- 0
          m[M!=a] <-2
        }
      } else {

        if(sum(c(-1,0,1) %in% M)==3){

          code <- c(print(c((M[M==1][1]),("is recoded to 2"))),
                    print(c((M[M==0][1]),("is recoded to 1"))),
                    print(c((M[M==-1][1]),(("is recoded to 0")))))

          m[M==1] <- 2
          m[M==0] <- 1
          m[M==-1] <- 0

        } else{

          code <- c(print(c((M[M==1][1]),("is recoded to 2"))),
                    print(c((M[M==0][1]),("is not recoded"))))

          m[M==1] <- 2
          m[M==0] <- 0
        }
      }
    }

    return(m)
  }
  m <- Recodemarkers(M)

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
  G_ERRBLUP <- Gall(m, cores)
  G <- G_ERRBLUP$G
  rownames(G) <- rownames(M)
  colnames(G) <- rownames(M)


  ERRBLUP_Bivar_Stepwise <- function(y1 , y2, G_ERRBLUP, iters=20, tolparinv= 1e-06) {

    Y_pred <- data.frame(Name = names(y1), y1 = y1)
    Y_add <- data.frame(Name = names(y2), y2 = y2)
    Y <- base::merge(Y_pred, Y_add, by="Name", all=TRUE)
    Y$Name <- as.factor(Y$Name)

    G_ERRBLUP <- G_ERRBLUP[rownames(G_ERRBLUP) %in% Y[,1], colnames(G_ERRBLUP) %in% Y[,1]]


    Sommer_function <- sommer::mmer(cbind(y1,y2) ~ 1, tolparinv= tolparinv,
                                    random=~ sommer::vs(Name, Gu=G_ERRBLUP, Gtc=sommer::unsm(2)),
                                    rcov=~ sommer::vs(units, Gtc=sommer::unsm(2)),
                                    data=Y, iters=iters)


    Pred_random <- Sommer_function$U
    Pred_fix <- Sommer_function$Beta

    Y1_pred <- data.frame(Name = names(Pred_random$`u:Name`$y1), y1 = Pred_random$`u:Name`$y1+Pred_fix[1,3])

    Y2_est <- data.frame(Name = names(Pred_random$`u:Name`$y2),y2 = Pred_random$`u:Name`$y2+Pred_fix[2,3])


    prediction <- base::merge(Y1_pred, Y2_est, by = "Name")


    return(prediction)

  }

  prediction <- ERRBLUP_Bivar_Stepwise(y1, y2, G, iters=iters, tolparinv= tolparinv)

  out <- list(Recodedmarkers = m, Relationshipmatrix = G_ERRBLUP, Predictions = prediction)

  return(out)

}




