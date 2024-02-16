#' a small functions that calculates and print the progress
#'
#' @param x integer, current iteration step
#' @param max integer, total iteration steps
#' @param start.time double, starting time
#'
#' @return rest time
#' @export progress
#'
#' @examples NULL
progress <- function (x, max = 100, start.time) {
  percent <- x / max * 100
  unit <- "secs"
  used.time <- difftime(Sys.time(), start.time, units = unit)
  finish.time <- (100 - percent) * used.time / percent
  if (finish.time < 121) {
    finish.time <- finish.time
  }
  else if (finish.time < 7201) {
    unit <- "mins"
    used.time <- difftime(Sys.time(), start.time, units = unit)
    finish.time <- (100 - percent) * used.time / percent
  }
  else {
    unit <- "hours"
    used.time <- difftime(Sys.time(), start.time, units = unit)
    finish.time <- (100 - percent) * used.time / percent
  }
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent)), " finish in: ", paste0(round(finish.time, digits = 0)), " ", unit)
  if (x == max)
    cat('\n')
}
