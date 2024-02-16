
#' Title
#'
#' @param x cf.mod
#' @param ...
#'
#' @return NULL
#' @export plot.cf.mod
#'
#' @examples NULL

plot.cf.mod <- function(x, ...){

  y_all_plot <- reshape2::melt(as.data.frame(x), id = "time")

  ggplot(y_all_plot, aes(x = time, y = value, group = variable)) +
    geom_line(aes(linetype = variable)) +
    scale_linetype_manual(values =  c("solid", "longdash")) +
    xlab("") + ylab(" ") + theme_bw()  + theme(legend.position = "none")
}
