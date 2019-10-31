
#' @title Marker Recoding Function
#'
#' @description Function to recode marker data from {-1, 0, 1} to {0, 1, 2} or from {0, 1} to {0, 2} marker coding
#'
#' @param M A {-1, 0, 1} or {0, 1} coded Marker matrix
#'
#' @return A {0, 1, 2} or {0, 2} coded marker matrix
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' m <- Recodemarker(wheat.X)
#'
#' @export
#'


Recodemarker <- function(M){



  if(sum(c(-1,0,1) %in% M)==3){


    M[M==1] <- 2
    M[M==0] <- 1
    M[M==-1] <- 0

  }  else{


    M[M==1] <- 2
    M[M==0] <- 0

  }

  return(M)
}




