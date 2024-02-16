givens.rotation <- function(i, j, theta, K){
  erg <- diag(K)
  erg[i,i] <- cos(theta) %>% as.numeric()
  erg[i,j] <- -sin(theta)%>% as.numeric()
  erg[j,i] <- sin(theta) %>% as.numeric()
  erg[j,j] <- cos(theta) %>% as.numeric()
  erg
}

#' Calculate (a sequence of) Givens Rotation matrices
#'
#' @param theta Vector, containing \code{choose(K,2)} rotation angles
#' @param K Integer, dimension of decomposition matrix
#'
#' @return a (KxK) rotation matrix
#' @export givens.Q
#'
#' @examples NULL
givens.Q <- function(theta, K){
  if (K < 2) {
    stop("The smallest dimension needs to be 2.")
  }
  erg <- diag(K)
  c <- 1
  for (i in 1:(K-1)) {
    for (j in (i + 1):K) {
      erg <- erg %*% givens.rotation(i, j, theta[c], K)
      c <- c + 1
    }
  }
  erg
}

#' From svars package: B matrix after givens rotation at giben angles,
#' faster than \code{givens.Q(theta, K)}
#'
#' @param theta Vector, containing \code{choose(K,2)} rotation angles
#' @param B Matrix, (KxK) dimensional original decomposition matrix, usually the cholesky factor
#'
#' @return a (KxK) rotated B matrix
#' @export givens.Q.svars
#'
#' @examples NULL
givens.Q.svars <- function(theta, B){
  K      <- (1 + sqrt(1 + 8 * length(theta))) / 2
  combns <- combn(K, 2, simplify = FALSE)
  erg    <- diag(K)

  for (i in seq_along(theta)) {

    temp <- diag(K)
    temp[combns[[i]][1], combns[[i]][1]] <- cos(theta[i])
    temp[combns[[i]][2], combns[[i]][2]] <- cos(theta[i])
    temp[combns[[i]][1], combns[[i]][2]] <- - sin(theta[i])
    temp[combns[[i]][2], combns[[i]][1]] <- sin(theta[i])
    erg <- erg %*% temp

  }
  tcrossprod(B, erg)
}
