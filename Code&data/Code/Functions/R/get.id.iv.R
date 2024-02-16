#' Identification based on Proxy SVAR
#'
#' @param Model An object of class 'vars', 'vec2var', 'nlVar'. Estimated VAR object
#' @param instruments Vector of dimension (Tx1) containing SINGLE instrument for ONE of structural shocks
#' with the same length as the time series
#' @param Synchro Logical, if Synchro = TRUE, instruments and VAR-residuals are synchronized. This only works if
#' the given instruments is an object of class "ts".
#'
#' @return  list object of class "my.id"
#' @export get.id.iv
#' @importFrom dplyr left_join
#' @examples NULL
get.id.iv <- function(Model, instruments = NULL, Synchro = FALSE){
  if (is.null(instruments)) {
    stop("Please provide a valid instrument")
  }

  p       <- Model$p
  K       <- Model$K
  u       <- resid(Model)
  Tob     <- Model$obs
  Varname <- names(Model$varresult)
  Beta    <- Model %>% vars::Bcoef() %>% t() # (Kp + mu) x K parameter matrix
  Z       <- Model$datamat[,-c(1:K)] %>% as.matrix() # (T-p) x (Kp + mu) covariates matrix
  Type    <- Model$type
  Covmat  <- crossprod(u) / (Tob - K*p - 1)

  iv      <- instruments

  if (Synchro){
    if(class(instruments) != 'ts'){
      stop("Automatic synchronization only works for time series object (ts). Please convert instruments into ts object!")
    }

    u_time  <- time(Model$y)
    u_time  <- u_time[-(1:p)]
    iv_time <- time(instruments)

    u.ts    <- data.frame("time" = u_time, "u" = u)
    iv.ts   <- data.frame("time" = iv_time, "iv" = instruments)
    syn.ts  <- suppressWarnings(left_join(u.ts, iv.ts, by = "time") %>% tidyr::drop_na())

    iv <- as.matrix(syn.ts$iv)
    u <- as.matrix(syn.ts[,c(2:(K+1))])
    Tob <- nrow(u)
  }

  if (Tob != length(iv)) {
    stop("length of instrument variable is unequal to data length")
  }

  Pi        <- solve(crossprod(u)/Tob) %*% (crossprod(u, iv)/Tob)
  phi2      <- (crossprod(iv, u)/Tob) %*% solve(crossprod(u)/Tob) %*% (crossprod(u, iv)/Tob)
  B_k       <- (crossprod(u, iv)/Tob) %*% (1/sqrt(phi2))
  s_errors  <- u %*% Pi %*% sqrt(phi2)^(-1)
  F_test    <- test.weak.IV(Tob = Tob, u = u, instruments = iv)

  row.names(B_k) <- Varname

  erg <- list("B" = B_k,
              "Varname" = Varname,
              "dat" = Model$y,
              "Instrument" = instruments,
              "F_test" = F_test,
              "p" = p,
              "Tob" = Model$obs,
              "u" = resid(Model),
              "epsilon" = s_errors,
              "SigmaU" = Covmat,
              "Beta" = Beta,
              "Z" = Z,
              "type" = Type,
              "Synchro" = Synchro,
              "method" = "iv")
  if (Synchro) {
    erg[["eps.ts"]] <- data.frame("time" = syn.ts$time, "epsilon" = s_errors)
  }

  class(erg) <- "my.id"

  return(erg)
}
