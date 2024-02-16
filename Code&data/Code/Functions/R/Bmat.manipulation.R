#' Resign a matrix such that diagonal elements are all positive
#'
#' @param Bmat Matrix, a (KxK) input matix
#' @param K Integer, dimension of matrix
#'
#' @return Matrix, a (KxK) output matix
#' @export resign.mat
#'
#' @examples NULL
resign.mat <- function(Bmat, K){
  for (k in 1:K) {
    if (Bmat[k,k] < 0){
      Bmat[,k] <- Bmat[,k] * -1
    }
  }
  Bmat
}

#' Reordering columns by minimizing the sqaured Frob. distance
#'
#' @param B.bench KxK benchmark decomposition Matrix
#' @param B.nonBench KxK decomposition Matrix which this function is applied to
#'
#' @return reordered KxK benchmark decomposition Matrix
#' @export Reorder.Frob
#'
#' @examples NULL
Reorder.Frob <- function(B.bench, B.nonBench){
  K <- ncol(B.bench)
  perms <- combinat::permn(K)
  score <- rep(NA, length(perms))
  for (i in 1:length(perms)) {
    d <- resign.mat(B.nonBench[,perms[[i]]], K) - B.bench
    score[i] <- sum(diag(t(d) %*% d))
  }
  erg <- resign.mat(B.nonBench[,perms[[which.min(score)]]], K)
  colnames(erg)   <- colnames(B.bench)
  row.names(erg)  <- row.names(B.bench)
  return(erg)
}
