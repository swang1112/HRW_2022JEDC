#' Calculate a space of potential candidates of decomposition matrix
#'
#' @param B Matrix, (KxK) dimensional original decomposition matrix, usually the cholesky factor
#' @param J Integer, number of grid points of rotation angles [0, pi/2]. Function will generates
#'  J + 1 candidates of decomposition matrix
#' @param Rand Integer, number of simulated rotation angels from a uniform distribution in the
#'  interval [0, pi].
#'
#' @return An array containing potential candidates of decomposition matrix.
#' @export B.space
#'
#' @examples
B.space <- function(B, J = NULL, Rand = NULL){
  if (ncol(B) != nrow(B)) {
    stop("The decomposition matrix (B) needs to be a square matrix.")
  }
  K <- ncol(B)
  M <- choose(K, 2)
  if (is.null(Rand)){
    if (is.null(J)) {
      stop("Please specify the number of grid points of rotation angles!")
    }

    theta.grid <- list(seq(0, pi/2, length = (J+1)))
    theta.expand <- expand.grid(rep(theta.grid, M))
    N <- (J + 1)^M
    erg <- array(NA, dim = c(K, K, N))
    for (n in 1:N) {
      erg[,,n] <-  givens.Q.svars(theta = theta.expand[n,], B)
    }
  } else {
    erg <- array(NA, dim = c(K, K, Rand))
    for (n in 1:Rand) {
      theta.rand <- runif(n = M, min = 0, max = pi)
      erg[,,n] <- givens.Q.svars(theta = theta.rand, B)
    }
  }
  erg
}


