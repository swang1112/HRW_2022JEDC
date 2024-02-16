#' Return structural impulse response functions (Theta)
#'
#' @param Model A class of 'varest', 'svarest' or 'vec2var'
#' @param Step Integer, IRF horizons
#' @param Bmat A (KxK) decomposition matrix of covariance matrix of reduced form errors, such that Bmat*Bmat'= Sigma_U / Identification matrix
#' @param Integrate logical, If Integrate = TRUE, the IRFs will be integrated
#' @param Plot logical, If Plot = TRUE, a plot with generic plot function will be generated
#' @param Acoef A (KxKxp) array. If this is provided, IRF will be calculated based on the given VAR Coefficients.
#' @param Varname A (1xK) vector of strings, Set variable names manuelly
#' @param Normalize logical, if Normalize = True, size of the shocks will be normalized
#' @param Bcoef A (Kp + mu) x K dimensional VAR parameter vectors, where mu is the number of deterministic terms
#' @param Simple Logical, if Simple = TRUE, run only kernel function. This is particularly useful when function is called within other functions, e.g. bootstraps
#'
#' @return An array collects the structural impulse response functions (Theta)
#' @importFrom vars Phi
#' @export get.my.irf
#'
#' @examples NULL
#'
get.my.irf <- function(Model = NULL, Acoef = NULL, Bcoef = NULL, Simple = FALSE, Lag = NULL, Step,
                       Normalize = FALSE, Varname = NULL, Bmat = NULL, Integrate = FALSE, Plot = FALSE){

  if (Simple){
    if (is.null(Bmat)){Bmat <- diag(K)}
    if (Normalize){Bmat <- Bmat %*% diag(1/diag(Bmat))}

    # dimension of Bcoef (Kp + mu) x K
    K <- ncol(Bcoef)
    if (ncol(Bmat) == 1){Bmat <- matrix(Bmat, nrow = K, ncol = K)}
    p <- Lag
    Phi <- array(0, c(K, K, Step + 1))
    Theta <- Phi
    # recover A form
    A <- array(0, c(K, K, Step + 1))
    for (i in 1:p) {
      A[,,i] <- t(Bcoef[(1 + (i-1)*K):(i*K), ])
    }
    # recursive fomula for Phis
    Phi[,,1] <- diag(K)
    Theta[,,1] <- Bmat


    for (h in 2:(Step + 1)) {
      temp <- matrix(0, nrow = K, ncol = K)
      for (j in 1:(h-1)) {
        temp <- temp + Phi[,,h-j] %*% A[, , j]
      }
      Phi[,,h] <- temp
      if (Integrate){
        Theta[,,h] <- Theta[,,h-1] + temp %*% Bmat
      } else {
        Theta[,,h] <- temp %*% Bmat
      }
    }

    colnames(Theta) <- row.names(Theta) <- colnames(Bcoef)

    Theta

  } else {
    if (is.null(Model)) {

      if (is.null(Acoef)) {
        stop("Please put in VAR parameters !")
      }

      K <- dim(Acoef)[1]
      p <- dim(Acoef)[3]

      # reduced form MA parameter
      Phi      <- array(0, c(K, K, Step + 1))
      Phi[,,1] <- diag(K)

      # structural MA parameter
      Theta     <- array(0, c(K, K, Step + 1))
      if (is.null(Bmat)){Bmat <- diag(K)}
      if (Normalize){Bmat <- Bmat %*% diag(1/diag(Bmat))}
      Theta[,,1] <- Bmat

      # Acoef for iteration, A[,,i] = 0 for i > p
      A <- array(0,c(K, K, Step))
      for (i in 1:p) {
        A[,,i] <- Acoef[,,i]
      }

      for (i in 2:(Step + 1)) {
        summ <- matrix(0, nrow = K, ncol = K)
        for (j in 1:(i-1)) {
          summ <- summ + Phi[,,i-j] %*% A[,,j]
        }
        Phi[,,i] <- summ
        Theta[,,i] <- Phi[,,i] %*% Bmat
      }

      if (Integrate == TRUE){
        for (i in 2:(Step + 1)) {
          Theta[,,i] <- Theta[,,i-1] + Theta[,,i]
        }
      }
    } else {
      myclass <- c("varest", "vec2var")
      if (is(Model, myclass) == FALSE) {
        stop("Please put in a model from class 'varest' or 'vec2var' !")
      }
      if (class(Model) == "varest"){
        Varname <- names(Model$varresult)
        K <- Model$K
      } else if (class(Model) == "vec2var"){
        Varname <- colnames(Model$y)
        K <- Model$K
      }

      # reduced form MA parameter
      Phi <- vars::Phi(x = Model, Step)
      # structural MA parameter
      Theta <- array(0, c(K, K, Step + 1))
      if (is.null(Bmat)){Bmat <- diag(K)}
      if (Normalize){Bmat <- Bmat %*% diag(1/diag(Bmat))}
      Theta[,,1] <- Bmat
      if (Integrate == TRUE){
        for (i in 2:(Step + 1)) {
          Theta[,,i] <- Theta[,,i-1] + Phi[,,i] %*% Bmat
        }
      } else {
        for (i in 2:(Step + 1)) {
          Theta[,,i] <- Phi[,,i] %*% Bmat
        }
      }
    }

    # rename
    row.names(Theta) <- Varname
    # plot
    if (Plot == TRUE) {
      par(mfrow = c(K,K))
      for (i in 1:K) {
        for (j in 1:K) {
          plot(Theta[i,j,], xlab = "Horizont", ylab = "Response",
               main = bquote(epsilon[.(j)]~" -> " ~ .(Varname[i])), type = "l")
          abline(h = 0, lty = 2)}}
      par(mfrow = c(1,1))}

    return(Theta)
  }

}


