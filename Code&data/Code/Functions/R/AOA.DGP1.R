#' Calculate Area Of Acceptance based on DGP1 simulation which is then used to evaluate the identified IFRs.
#'
#' @param Simu Integer, number of simulations
#' @param Step Integer, IRF horizont
#' @param Ci Double, between 0 and 1, confidence interval
#'
#' @return A list contains lower and upper bounds arrays of dimension (KxKxh) respectively
#'
#' @export AOA.DGP1
#'
#' @examples NULL
AOA.DGP1 <- function(Simu, Step, Ci) {
  beta    <- rep(0.99, Simu)
  alpha   <- runif(Simu, 0.40, 0.60)
  k       <- runif(Simu, 0.03, 0.07)
  gamma   <- runif(Simu, 0.40, 0.60)
  delta_x <- runif(Simu, 0.05, 0.15)
  tau_x   <- runif(Simu, 0.30, 0.70)
  tau_pi  <- runif(Simu, 1.60, 2.00)
  tau_r   <- runif(Simu, 0.40, 0.80)
  rho_x   <- runif(Simu, 0.40, 0.60)
  rho_pi  <- runif(Simu, 0.40, 0.60)
  rho_r   <- runif(Simu, 0.40, 0.60)

  IRFs    <- array(0, c(3, 3, Step + 1, Simu))

  start.time <- Sys.time()
  cat("\r", "...calculating finish time...")

  for (i in 1:Simu) {
    VAR.temp <-
     DSGE2VAR.DGP1(alpha   = alpha[i],
                   beta    = beta[i],
                   k       = k[i],
                   gamma   = gamma[i],
                   delta_x = delta_x[i],
                   tau_x   = tau_x[i],
                   tau_pi  = tau_pi[i],
                   tau_r   = tau_r[i],
                   rho_x   = rho_x[i],
                   rho_pi  = rho_pi[i],
                   rho_r   = rho_r[i])

    A_array      <- array(NA, c(3, 3, 2))
    A_array[,,1] <- VAR.temp$A_1
    A_array[,,2] <- VAR.temp$A_2

    # impulse response
    IRFs[,,,i]   <- get.my.irf(Acoef = A_array, Step = Step, Normalize = TRUE, Varname = c("x", "pi", "r"), Bmat = VAR.temp$B)

    # print progress
    progress(i, max = Simu, start.time = start.time)
  }

  L <- U <- array(NA, c(3, 3, Step + 1))
  for (h in 1:(Step + 1)) {
    for (i in 1:3) {
      for (j in 1:3) {
        L[i,j,h] <- quantile(IRFs[i,j,h,], probs = c(0.5 * (1 - Ci), 0.5 * (1 + Ci)))[1]
        U[i,j,h] <- quantile(IRFs[i,j,h,], probs = c(0.5 * (1 - Ci), 0.5 * (1 + Ci)))[2]
      }
    }
  }

  erg <- list()
  erg[["L"]] <- L
  erg[["U"]] <- U
  return(erg)
}







