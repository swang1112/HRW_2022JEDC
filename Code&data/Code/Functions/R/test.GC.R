#' Test for Granger Causality
#'
#' @param x Class of "varest"
#' @param C1 a KxK matrix contains the restrictions
#'
#' @return test statistics and p value
#' @export test.GC
#'
#' @examples NULL
test.GC <- function(x, C1){
  K <- x$K
  p <- x$p
  Tob <- x$obs
  As <- Bcoef(x)[,1:(K*p)]
  nu <- Bcoef(x)[,-(1:(K*p))]
  beta_hat <- matrix(c(c(nu),c(As)), nrow = K^2*p+length(nu), ncol = 1)


  Cmat <- matrix(0, nrow = nrow(Bcoef(x)), ncol = ncol(Bcoef(x)))
  for (i in 1:p) {
    Cmat[, (K*(i-1)+1):(K*i)] <- C1
  }
  As_C <- Cmat[,1:(K*p)]
  nu_C <- Cmat[,-(1:(K*p))]
  Cvec <- matrix(c(c(nu_C),c(As_C)), nrow = K^2*p+length(nu), ncol = 1)

  Js <- which(Cvec == 1)
  N <- length(Js)
  Cs <- matrix(0, nrow = nrow(Cvec), ncol = N)
  for (j in 1:N) {
    Cs[Js[j], j] <- 1
  }

  C <- t(Cs)

  Sigma_u <- crossprod(resid(x)) / (Tob - K*p - 1)

  Z_temp <- x$datamat[,-c(1:K)]
  Z <-cbind(Z_temp[-c(1:(K*p))], Z_temp[,1:(K*p)])
  Z <- t(Z)

  Gamma = Z %*% t(Z) / Tob

  lambda_W <- Tob * t(C %*% beta_hat) %*% solve(C %*% kronecker(solve(Gamma), Sigma_u) %*% t(C)) %*% (C %*% beta_hat)
  lambda_F <- lambda_W / N
  p_W <- 1 - pchisq(lambda_W, df = N)
  p_F <- 1 - pf(lambda_F, df1 = N, df2 = K*(Tob - K*p - 1))

  erg <- list("TestStat" = data.frame("lambda_W" = lambda_W, "lambda_F" = lambda_F),
              "p_val" = data.frame("p_W" = p_W, "p_F"= p_F))

  return(erg)

}
