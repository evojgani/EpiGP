
#' @title Bivariate sERRBLUP Phenotype Prediction Function
#'
#' @description Function to do phenotype prediction based on desired proportion of pairwise SNP interactions in bivariate model
#'
#' @param M The original marker matrix of \code{{-1,0,1}}, \code{{0,1}}, \code{{0,1,2}}, \code{{0,2}} or character coded markers with named individuals in the rows and the markers in the columns
#' @param y1 A numeric vector of named phenotypes for the target triat to be predicted in bivariate sERRBLUP
#' @param y2 A numeric vector of named phenotypes for the additional trait to be estimated in bivariate sERRBLUP
#' @param k A desired proportion of SNP interactions to be included in the model
#' @param iters Maximum number of iterations allowed with the default value of \code{20} (The parameter of mmer function from sommer package)
#' @param tolparinv Tolerance parameter for matrix inverse used when singularities are encountered in the estimation procedure with the default value of \code{1e-06} (The parameter of mmer function from sommer package)
#' @param cores The number of cores with the default value of \code{1}
#'
#'@return A list of eight components:
#'
#' \describe{
#'   \item{Recodedmarkers}{A \code{{0,1,2}} or \code{{0,2}} coded marker matrix}
#'   \item{Relationshipmatrix}{A list of two components: ERRBLUP relationship matrix (G) and a vector of all genotype combinations frequencies in the population (P)}
#'   \item{Effect.Relationshipmatrix}{sERRBLUP Relationship matrix when \code{k} percent of pairwise SNP interactions maintained in the model based on absolute effect size estimations}
#'   \item{Var.Relationshipmatrix}{sERRBLUP Relationship matrix when \code{k} percent of pairwise SNP interactions maintained in the model based on effect variance estimations}
#'   \item{Effect}{A numeric vector of all estimated pairwise SNP interaction effect sizes}
#'   \item{Effect.Variance}{A numeric vector of all estimated pairwise SNP interaction effect variances}
#'   \item{Effect.Predictions}{A dataframe of both phenotype predictions of target triat test set, phenotype estimations of target triat training set and phenotype estimations of additional trait full set in bivarite sERRBLUP method based on effect sizes}
#'   \item{Var.Prediction}{A dataframe of both phenotype predictions of target triat test set, phenotype estimations of target triat training set and phenotype estimations of additional trait full set in bivarite sERRBLUP method based on effect variances}
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
#' sERRBLUP_Bivar <- sERRBLUP_BV(M, y1, y2, 10, iters=20, tolparinv= 1e-06, cores=15)
#'
#' @export
#'


sERRBLUP_BV <- function(M, y1, y2, k, iters=20, tolparinv= 1e-06, cores=1){

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
  G1 <- G_ERRBLUP$G
  P <- G_ERRBLUP$P

  SNP_Effect_Var <- function(m, Pheno, G_ERRBLUP, P, cores=1){


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

    out <- list(Effect = u_hat, Effect.Var = sigma_hat)
    return(out)


  }
  estimations <- SNP_Effect_Var(m, y2, G1, P, cores)
  t_hat <- estimations$Effect
  sigma_hat <- estimations$Effect.Var


  Gtop <- function(m, Estimations, k, cores=1){


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
      include <- integer((nsnp*(nsnp+1)*2))+1L

      include[abs(Estimations)< stats::quantile(abs(Estimations),(1-(k/100)))] <- 0L

      rm(Estimations)

      check = prod(include)


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


        if(check!=1){
          index2 = nsnp - index + 1

          activ = c(1:((nsnp-index+1)) + (index-1)*4*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + (index2-1)*4*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + (nsnp-index+1) + (index-1)*4*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + (nsnp-index2+1) + (index2-1)*4*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + 2*(nsnp-index+1) + (index-1)*4*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + 2*(nsnp-index2+1) +(index2-1)*4*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + 3*(nsnp-index+1) + (index-1)*4*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + 3*(nsnp-index2+1) + (index2-1)*4*nsnp - 2 * (index2-1) * (index2-2))



          Z_share <- matrix(include[activ], ncol=nindi, nrow=(nsnp+1)*4, byrow=FALSE) * Z_share
        }


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

          p_i[((index-1)*nsnp*4+1):((index)*nsnp*4)][c(1:(nsnp-index+1), 1:(nsnp-index+1)+(nsnp), 1:(nsnp-index+1)+2*(nsnp), 1:(nsnp-index+1) + 3*(nsnp))] <-
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

      include <- integer((nsnp*(nsnp+1)*9/2))+1L

      include[abs(Estimations)< stats::quantile(abs(Estimations),(1-(k/100)))] <- 0L

      rm(Estimations)

      check = prod(include)


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


        if(check!=1){
          index2 = nsnp - index + 1

          activ = c(1:((nsnp-index+1)) + (index-1)*9*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + (nsnp-index+1) + (index-1)*9*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + (nsnp-index2+1) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + 2*(nsnp-index+1) + (index-1)*9*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + 2*(nsnp-index2+1) +(index2-1)*9*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + 3*(nsnp-index+1) + (index-1)*9*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + 3*(nsnp-index2+1) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + 4*(nsnp-index+1) + (index-1)*9*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + 4*(nsnp-index2+1) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + 5*(nsnp-index+1) + (index-1)*9*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + 5*(nsnp-index2+1) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + 6*(nsnp-index+1) + (index-1)*9*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + 6*(nsnp-index2+1) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + 7*(nsnp-index+1) + (index-1)*9*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + 7*(nsnp-index2+1) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2),
                    1:((nsnp-index+1)) + 8*(nsnp-index+1) + (index-1)*9*nsnp - 2 * (index-1) * (index-2),
                    1:((nsnp-index2+1)) + 8*(nsnp-index2+1) + (index2-1)*9*nsnp - 2 * (index2-1) * (index2-2))



          Z_share <- matrix(include[activ], ncol=nindi, nrow=(nsnp+1)*9, byrow=FALSE) * Z_share
        }


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


    p_i <- p_i[stats::complete.cases(p_i)]

    G_k <- G / (2 * sum(p_i*(1-p_i)))

    return(G_k)

  }
  G2 <- Gtop(m, t_hat, k, cores)
  G3 <- Gtop(m, sigma_hat, k, cores)

  rownames(G2) <- rownames(m)
  colnames(G2) <- rownames(m)

  rownames(G3) <- rownames(m)
  colnames(G3) <- rownames(m)


  sERRBLUP_Bivar_Stepwise <- function(y1, y2, Gtop, iters=20, tolparinv= 1e-06) {


    Y_pred <- data.frame(Name = names(y1), y1 = y1)
    Y_add <- data.frame(Name = names(y2), y2 = y2)
    Y <- base::merge(Y_pred, Y_add, by="Name", all=TRUE)
    Y$Name <- as.factor(Y$Name)

    Gtop <- Gtop[rownames(Gtop) %in% Y[,1], colnames(Gtop) %in% Y[,1]]


    Sommer_function <- sommer::mmer(cbind(y1,y2) ~ 1, tolparinv= tolparinv,
                                    random=~ sommer::vs(Name, Gu=Gtop, Gtc=sommer::unsm(2)),
                                    rcov=~ sommer::vs(units, Gtc=sommer::unsm(2)),
                                    data=Y, iters=iters)


    Pred_random <- Sommer_function$U
    Pred_fix <- Sommer_function$Beta

    Y1_pred <- data.frame(Name = names(Pred_random$`u:Name`$y1), y1 = Pred_random$`u:Name`$y1+Pred_fix[1,3])

    Y2_est <- data.frame(Name = names(Pred_random$`u:Name`$y2),y2 = Pred_random$`u:Name`$y2+Pred_fix[2,3])


    prediction <- base::merge(Y1_pred, Y2_est, by = "Name")

  }
  prediction_effect <- sERRBLUP_Bivar_Stepwise(y1, y2 , G2)
  prediction_var <- sERRBLUP_Bivar_Stepwise(y1, y2 , G3)

  out <- list(Recodedmarkers = m, Relationshipmatrix = G_ERRBLUP,
              Effect.Relationshipmatrix = G2, Var.Relationshipmatrix = G3,
              Effect = t_hat, Effect.Variance = sigma_hat,
              Effect.Predictions = prediction_effect, Var.Prediction = prediction_var)

  return(out)

}
