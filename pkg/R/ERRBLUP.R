
#' @title ERRBLUP Phenotype Prediction Function
#'
#' @description Function to do phenotype prediction based on all pairwise SNP interactions
#'
#' @param M The original marker matrix of \code{{-1,0,1}}, \code{{0,1}}, \code{{0,1,2}}, \code{{0,2}} or character coded markers
#' @param Pheno A numeric vector of phenotypes
#' @param cores The number of cores with the default value of \code{1}
#'
#'@return A list of three components:
#'
#' \describe{
#'   \item{Recodedmarkers}{A \code{{0,1,2}} or \code{{0,2}} coded marker matrix}
#'   \item{Relationshipmatrix}{A list of two components: ERRBLUP relationship matrix (G) and a vector of all genotype combinations frequencies in the population (P)}
#'   \item{Predictions}{A numeric vector of both phenotype estimations of training set and phenotype predictions of test set based on ERRBLUP method}
#' }
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' M <- wheat.X
#' N <- length(Phenotype)
#' n <- 60
#' test <- sample(1:N,n)
#' pheno <- Phenotype
#' pheno[test] <- NA
#' ERRBLUP <- ERRBLUP(M, pheno)
#'
#' @export
#'


ERRBLUP <- function(M, Pheno, cores=1){

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

  ERRBLUP <- function(Pheno , G_ERRBLUP) {


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
  prediction <- ERRBLUP(Pheno , G)

  out <- list(Recodedmarkers = m, Relationshipmatrix = G_ERRBLUP, Predictions = prediction)

  return(out)

}




