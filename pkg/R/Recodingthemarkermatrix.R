
#' @title Marker Recoding Function
#'
#' @description Function to recode marker matrix to {0, 2} or {0, 1, 2} coded marker matrix
#'
#' @param M A {-1, 0, 1} or {0, 1} or character coded marker matrix
#'
#' @return A {0, 1, 2} or {0, 2} coded marker matrix. The function reports how the recoding was implemented in the input marker matrix
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' m <- Recodemarkers(wheat.X)
#'
#' @export
#'


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




