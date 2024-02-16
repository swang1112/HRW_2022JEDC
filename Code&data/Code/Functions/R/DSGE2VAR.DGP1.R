# function to solve
optim.fct1 <- function(Par, alpha, beta, k, gamma, delta_x, tau_x, tau_pi, tau_r, rho_x, rho_pi, rho_r){
  # DSGE Model: Gamma_0 %*% z_t = Cmat %*% E_t(z_t+1) + Gamma_1 %*% z_t-1 + Hmat %*% w_t
  Gamma_0 <- matrix(c(1, 0, delta_x,
                      -k, 1, 0,
                      (tau_r - 1)*tau_x, (tau_r - 1)*tau_pi, 1), nrow = 3, ncol = 3, byrow = T)
  Cmat    <- matrix(c(gamma, delta_x, 0,
                      0,  beta*(1 + alpha * beta)^-1, 0,
                      0, 0, 0), nrow = 3, ncol = 3,  byrow = T)
  Gamma_1 <- matrix(c(1-gamma, 0,0,
                      0, alpha * (1 + alpha * beta)^-1, 0,
                      0, 0, tau_r), nrow = 3, ncol = 3,  byrow = T)
  Hmat    <- diag(3)

  # VAR Model: z_t = PHImat %*% z_t+1 + Bmat %*% w_t
  # E_t(z_t+1) = PHImat %*% z_t + Bmat %*% Fmat %*% w_t
  Fmat   <-  diag(c(rho_x, rho_pi, rho_r))
  PHImat <- matrix(Par[1 : 9], nrow = 3, ncol = 3)
  Bmat   <- matrix(Par[10:18], nrow = 3, ncol = 3)

  obj.eq1 <- c(Gamma_1) - kronecker(diag(3), Gamma_0) %*% c(PHImat) + kronecker(t(PHImat) %*% t(PHImat), diag(3)) %*% c(Cmat)
  obj.eq2 <- kronecker(diag(3), Gamma_0) %*% c(Bmat) - kronecker(t(Bmat) %*% t(PHImat), diag(3)) %*%c(Cmat) - kronecker(t(Fmat) %*% t(Bmat), diag(3)) %*% c(Cmat) - c(diag(3))
  erg     <- rbind(obj.eq1, obj.eq2) %>% as.numeric()

  return(erg)
}

#' Solve a 3 dimensional DSGE model and convert it into VAR(2)
#'
#' @param alpha
#' @param beta
#' @param k
#' @param gamma
#' @param delta_x
#' @param tau_x
#' @param tau_pi
#' @param tau_r
#' @param rho_x
#' @param rho_pi
#' @param rho_r
#'
#' @return A list containing VAR and structural parameters.
#'
#' @import dplyr
#' @importFrom nleqslv nleqslv
#' @export DSGE2VAR.DGP1
#'
#' @examples NULL
DSGE2VAR.DGP1 <- function(alpha, beta, k, gamma, delta_x, tau_x, tau_pi, tau_r, rho_x, rho_pi, rho_r){
  # parameter initialization
  Par <- rep(0, 18)
  Par_erg <- nleqslv(Par, optim.fct1,
                              alpha = alpha,
                              beta = beta,
                              k = k,
                              gamma = gamma,
                              delta_x = delta_x,
                              tau_x = tau_x,
                              tau_pi = tau_pi,
                              tau_r = tau_r,
                              rho_x = rho_x,
                              rho_pi = rho_pi,
                              rho_r = rho_r)$x %>% round(2)
  # VAR parameters
  PHImat <- matrix(Par_erg[1 : 9], nrow = 3, ncol = 3)
  Bmat   <- matrix(Par_erg[10:18], nrow = 3, ncol = 3)
  Fmat   <- diag(c(rho_x, rho_pi, rho_r))

  # VAR(2) representation of DSGE
  A_1       <-   PHImat + Bmat %*% Fmat %*% solve(Bmat)
  A_2       <- - Bmat %*% Fmat %*% solve(Bmat) %*% PHImat

  erg          <- list()
  erg[["Phi"]] <- PHImat
  erg[["B"]]   <- Bmat
  erg[["A_1"]] <- A_1
  erg[["A_2"]] <- A_2
  return(erg)
}



