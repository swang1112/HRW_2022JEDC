#' relative mean squared error
#'
#' @param B Matrix of dimension (KxK), the true structural decomposition
#' @param Bhat Matrix of dimension (KxK), the estimated decomposition
#' @param xcol Numeric, if provided, only the specific column in matrix Bhat is evaluated
#'
#' @return Numeric, relative mean squared error
#' @export eva.RMSE
#'
#' @examples NULL
eva.RMSE <- function(B, Bhat, xcol = NULL){

    if (is.null(xcol)){
      erg <- sqrt(((B - Bhat) / B)^2)
    } else {
      erg <- sqrt(((B[,xcol] - Bhat[,xcol]) / B[,xcol])^2 )
    }

    erg
}

#' check signs / labeling ratio
#'
#' @param Bhat Matrix of dimension (KxK), the estimated decomposition
#' @param TestSign Matrix of dimension (KxK), containing the true Sign pattern in DGP 1: positive sign; -1: negative sign
#' @param Permut Logical, whether apply function to alternative column orderings
#'
#' @return Logical
#' @export eva.sign
#' @importFrom combinat permn
#'
#' @examples NULL
eva.sign <- function(Bhat, TestSign, Permut = FALSE){

  if (!Permut){
    erg   <- all(sign(Bhat) == TestSign, na.rm = T)
  } else {

    K     <- nrow(Bhat)
    perms <- permn(K)
    erg   <- FALSE
    i     <- 1
    while (erg == 0 & i <= length(perms)) {
      erg <- all(sign(resign.mat(Bhat[,perms[[i]]], K = K)) == TestSign, na.rm = T)
      i   <- i + 1
    }
  }
  return(erg)
}

#' check unique monetary policy shock (works only for K = 3) based on sign pattern (-1,-1,1)
#'
#' @param B Matrix of dimension (KxK), the true structural decomposition
#' @param Reorder Logical, whether column ordering of Bhat is unidentified
#' @param Bhat Matrix of dimension (KxK), the estimated decomposition
#'
#' @return Numeric, vector of dimension (1x4): c(indicator, RMSEs for demand, supply and mp shock)
#' @export eva.sign.UMP
#'
#' @examples NULL
eva.sign.UMP <- function(Bhat, B, Reorder = FALSE){
  if (Reorder){
    mps <- rep(NA, 3)
    for (k in 1:3) {
      ifelse(all(Bhat[2,k] < 0, Bhat[3,k] > 0) | all(Bhat[2,k] > 0, Bhat[3,k] < 0),
             mps[k] <- 1, mps[k] <- 0)
    }

    if (sum(mps) == 0) {
      return(c(0, rep(NA, 3))) # no monetary shock identified!
    } else if (sum(mps) > 1){
      return(c(0, rep(NA, 3))) # not unique!
    } else {
      xcol <- which(mps == 1)
      if (Bhat[3, xcol] < 0) {Bhat[,xcol] <- Bhat[,xcol] * -1}
      rmse <- eva.RMSE(B = B, Bhat = Bhat, xcol = xcol)
      return(c(1, rmse))
    }
  } else {
    mps <- 0
    ifelse(all(Bhat[2,3] < 0, Bhat[3,3] > 0) | all(Bhat[2,3] > 0, Bhat[3,3] < 0),
           mps <- 1, mps <- 0)
    if (mps == 0) {
      return(c(0, rep(NA, 3)))
    } else {
      if (Bhat[3,3] < 0) {Bhat[,3] <- Bhat[,3] * -1}
      rmse <- eva.RMSE(B = B, Bhat = Bhat, xcol = 3)
      return(c(1, rmse))
    }
  }
}
