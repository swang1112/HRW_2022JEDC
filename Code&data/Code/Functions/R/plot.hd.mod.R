#' plot function for historical decompostion
#'
#' @param x historical decomposition
#'
#' @return plot of historical decomposition
#' @export plot.hd.mod
#'
#' @examples NULL
plot.hd.mod <- function(x){
  reshape2::melt(x, id = "t") %>% ggplot(aes(x = t, y = value)) +
    geom_line() +
    facet_wrap(~variable, ncol = 1) +
    xlab("Time") + theme_bw()
}
