#' Test for weak instruments in Lunsford 2016
#'
#' @param Tob Numeric, length of simulated time series
#' @param u Vector, a (Tx1) vector containing residuals from reduced form VAR
#' @param instruments Vector of dimension (Tx1) containing SINGLE instrument for ONE of structural shocks
#' with the same length as the time series
#'
#' @return Numeric, test statistic. For critical value: look at Lunsford 2016: p.16 Table 2
#' @export test.weak.IV
#'
#' @examples NULL
test.weak.IV <- function(Tob, u, instruments){
  u      <- u
  Pi     <- solve(crossprod(u)/Tob) %*% (crossprod(u, instruments)/Tob)
  K      <- ncol(u)

  # Test statistic Lunsford 2016: p.15 equation 29
  F_stat <- ((Tob - K)/K) * (crossprod(instruments) - crossprod((instruments - u %*% Pi))) /  crossprod((instruments - u %*% Pi))
  return(F_stat)
}


