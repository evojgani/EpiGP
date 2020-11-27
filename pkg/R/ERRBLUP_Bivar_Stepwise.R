
#' @title Bivariate ERRBLUP Phenotype Prediction Function Relying On The Out Put Of ERRBLUP Relationship Matrix Function
#'
#' @description Function to do phenotype prediction based on all pairwise SNP interactions in bivariate model
#'
#' @param y1 A numeric vector of named phenotypes for the target triat to be predicted in bivariate ERRBLUP
#' @param y2 A numeric vector of named phenotypes for the additional trait to be estimated in bivariate ERRBLUP
#' @param G_ERRBLUP ERRBLUP relationship matrix with the row names and column names of all individuals
#' @param iters Maximum number of iterations allowed with the default value of \code{20} (The parameter of mmer function from sommer package)
#' @param tolparinv Tolerance parameter for matrix inverse used when singularities are encountered in the estimation procedure with the default value of \code{1e-06} (The parameter of mmer function from sommer package)

#' @return A dataframe of both phenotype predictions of target triat test set, phenotype estimations of target triat training set and phenotype estimations of additional trait full set in bivarite ERRBLUP method
#'
#' @examples
#' library(BGLR)
#' data(wheat)
#' m <- Recodemarkers(wheat.X)
#' y1 <- Phenotype_BV[,1]
#' N <- length(y1)
#' n <- 60
#' test <- sample(1:N,n)
#' y1[test] <- NA
#' y2 <- Phenotype_BV[,2]
#' G_ERRBLUP <- Gall(m, cores=15)
#' G <- G_ERRBLUP$G
#' rownames(G) <- names(y1)
#' colnames(G) <- names(y1)
#' ERRBLUP_Bivar <- ERRBLUP_BV_Stepwise(y1, y2, G, iters=20, tolparinv= 1e-06)
#'
#' @export
#'


ERRBLUP_BV_Stepwise <- function(y1 , y2, G_ERRBLUP, iters=20, tolparinv= 1e-06) {


  Y_pred <- data.frame(Name = names(y1), y1 = y1)
  Y_add <- data.frame(Name = names(y2), y2 = y2)
  Y <- base::merge(Y_pred, Y_add, by = "Name",all = TRUE)
  Y$Name <- as.factor(Y$Name)

  G_ERRBLUP <- G_ERRBLUP[rownames(G_ERRBLUP) %in% Y[,1], colnames(G_ERRBLUP) %in% Y[,1]]


  Sommer_function <- sommer::mmer(cbind(y1,y2) ~ 1, tolparinv= tolparinv,
                                  random=~ sommer::vs(Name, Gu=G_ERRBLUP, Gtc=sommer::unsm(2)),
                                  rcov=~ sommer::vs(units, Gtc=sommer::unsm(2)),
                                  data=Y,iters=iters)


 Pred_random <- Sommer_function$U
 Pred_fix <- Sommer_function$Beta

 Y1_pred <- data.frame(Name = names(Pred_random$`u:Name`$y1), y1 = Pred_random$`u:Name`$y1+Pred_fix[1,3])

 Y2_est <- data.frame(Name = names(Pred_random$`u:Name`$y2),y2 = Pred_random$`u:Name`$y2+Pred_fix[2,3])


 prediction <- base::merge(Y1_pred, Y2_est, by = "Name")


 return(prediction)

}




