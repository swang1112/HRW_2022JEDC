#' Independence-based identification based on distance covariance
#'
#' @param Model  A class of 'varest' or 'vec2var'
#'
#' @return list object of class "my.id"
#' @importFrom steadyICA steadyICA
#' @export get.id.dcov
#'
#' @examples NULL
get.id.dcov <- function(Model){
  p       <- Model$p
  K       <- Model$K
  u       <- resid(Model)
  Tob     <- Model$obs
  Varname <- names(Model$varresult)
  Beta    <- Model %>% vars::Bcoef() %>% t() # (Kp + mu) x K parameter matrix
  Z       <- Model$datamat[,-c(1:K)] %>% as.matrix() # (T-p) x (Kp + mu) covariates matrix
  Type    <- Model$type
  Covmat  <- crossprod(u) / (Tob - K*p - 1)

  # Choleski decomp
  B_chol  <- t(chol(Covmat))
  u_test <- t(solve(B_chol)%*%t(u))

  rotat <- suppressMessages(steadyICA(u_test, symmetric=TRUE))
  Bmat <- B_chol %*% rotat$W

  rownames(Bmat) <- Varname

  erg <- list("B" = Bmat,
              "Varname" = Varname,
              "dat" = Model$y,
              "p" = p,
              "Tob" = Tob,
              "u" = u,
              "SigmaU" = Covmat,
              "Beta" = Beta,
              "Z" = Z,
              "type" = Type,
              "method" = "dcov")

  class(erg) <- "my.id"

  return(erg)
}
