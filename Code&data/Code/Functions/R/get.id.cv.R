#' Cloned from svars: Identification of SVAR models based on Changes in volatility (CV)
#'
#' @param Model An object of class 'vars', 'vec2var', 'nlVar'. Estimated VAR object
#' @param SB Integer, vector or date character. The structural break is specified either by an integer (number of observations in the pre-break period),
#'                    a vector of ts() frequencies if a ts object is used in the VAR or a date character. If a date character is provided, either a date vector containing the whole time line
#'                    in the corresponding format (see examples) or common time parameters need to be provided
#' @param dateVector Vector. Vector of time periods containing SB in corresponding format
#' @param start Character. Start of the time series (only if dateVector is empty)
#' @param end Character. End of the time series (only if dateVector is empty)
#' @param frequency Character. Frequency of the time series (only if dateVector is empty)
#' @param format Character. Date format (only if dateVector is empty)
#' @param restriction_matrix Matrix. A matrix containing presupposed entries for matrix B, NA if no restriction is imposed (entries to be estimated). Alternatively, a K^2*K^2 matrix can be passed, where ones on the diagonal designate unrestricted and zeros restricted coefficients. (as suggested in Luetkepohl, 2017, section 5.2.1).
#'
#' @return
#' @export
#'
#' @examples
get.id.cv <- function(Model, SB, start = NULL, end = NULL, frequency = NULL,
                      format = NULL, dateVector = NULL, max.iter = 50, crit = 0.001,
                      restriction_matrix = NULL){
  p       <- Model$p
  K       <- Model$K
  u       <- resid(Model)
  Tob     <- Model$obs
  Varname <- names(Model$varresult)
  Beta    <- Model %>% vars::Bcoef() %>% t() # (Kp + mu) x K parameter matrix
  Z       <- Model$datamat[,-c(1:K)] %>% as.matrix() # (T-p) x (Kp + mu) covariates matrix
  Type    <- Model$type
  Covmat  <- crossprod(u) / (Tob - K*p - 1)

  B_chol  <- t(chol(Covmat))

  erg_svars <- svars::id.cv(x = Model, SB = SB, start = start, end = end,
                           frequency = frequency, format = format,
                           dateVector = dateVector, max.iter = max.iter, crit = crit)

  Bmat <- erg_svars$B
  rownames(Bmat) <- Varname

  erg <- list("B" = Bmat,
              "Varname" = Varname,
              "dat" = Model$y,
              "svars" = erg_svars,
              "p" = p,
              "Tob" = Tob,
              "u" = u,
              "SigmaU" = Covmat,
              "Beta" = Beta,
              "Z" = Z,
              "type" = Type,
              "method" = "cv")

  class(erg) <- "my.id"

  return(erg)
}
