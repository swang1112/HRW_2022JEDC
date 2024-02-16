#' Remove Autocorrelation in IV
#'
#' @param iv Instruments, ts object
#' @param lag
#'
#' @return
#' @export Remove_AC
#'
#' @examples NULL
Remove_AC <- function(iv, lag){
  time <- time(iv) %>% as.numeric()
  iv_mat <- matrix(NA, nrow = length(iv), ncol = lag+1)
  for (i in 1:lag) {
    iv_mat[,i] <- lag(as.numeric(iv),i)
  }
  iv_mat[,lag+1] <- iv
  erg <- lm(iv_mat[,lag+1]~iv_mat[,-(lag+1)]-1) %>% resid
  for (i in 1:lag) {
    print(Box.test(erg, lag = i, type = "Ljung-Box"))
  }
  erg <- ts(erg, frequency = 4, start = c(floor(time[1+lag]), (time[1+lag] %% 1) / 0.25 + 1 ))
  return(erg)

}
