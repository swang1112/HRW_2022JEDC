#' Forecast error variance decomposition for SVAR Models
#'
#' @param My.id Object of class "my.id", identified SVAR
#' @param Step Calculation steps
#' @param Plot Logical, plot is generated if Plot = TRUE
#' @param Simple Logical, if Simple = TRUE, run only kernel function based on the provided IRF.
#' @param IRF Array, of dimension KxKx(Step+1)
#'
#' @return Kx1 list containing FEVD
#' @export get.my.fevd
#' @import ggplot2
#'
#' @examples NULL
#'
get.my.fevd <- function(My.id = NULL, Step = NULL, Plot = FALSE, Simple = F, IRF = NULL){
  if (Simple) {
    if (is.null(IRF)) {
      stop("Missing IRF input!")
    }

    K       <- dim(IRF)[1]
    Step    <- dim(IRF)[3]-1
    Varname <- row.names(IRF)

    Theta.temp <- IRF

  } else {
    if(class(My.id) != "my.id"){
      stop("Input should be an identified SVAR model of class 'my.id'!")
    }

    Varname <- My.id$Varname
    Method  <- My.id$method
    Tob     <- My.id$Tob
    K       <- length(Varname)
    p       <- My.id$p
    Beta    <- My.id$Beta
    B       <- My.id$B
    Z       <- My.id$Z

    if (Method[1] == "iv"){
      tryCatch(stop("Not supported for partial identification by Proxy SVAR."))
    }

    if (Method[1] != "sr"){
      stop("For other id.methods, please use fevd from svars package!")
    }

    Theta.temp <- get.my.irf(Simple = T, Bcoef = Beta, Bmat = My.id$IRF.M[,,1], Step = Step, Lag = p)
  }

  fe <- list()
  for(i in 1:K){
    fe[[i]] <- as.data.frame(t(Theta.temp[i,,]))
    colnames(fe[[i]]) <- Varname
  }
  names(fe) <- Varname
  fe2 <- fe

  for(i in 1:K){
    for(j in 1:(Step + 1)){
      fe2[[i]][j,] <- (colSums(fe[[i]][j:1,]^2)/sum(fe[[i]][j:1,]^2))*100
    }
  }

  if (Plot){
    fe.plot <- data.frame(V1 = seq(1, Step+1), Variables = c(sapply(Varname, rep, Step+1)), sapply(fe2, unlist))
    for(i in 3:ncol(fe.plot)){
      colnames(fe.plot)[i] <- paste("FEVD of", colnames(fe.plot)[i])
    }
    fe.plot <- melt(fe.plot, id = c('V1', 'Variables'))
    plot(ggplot(fe.plot, aes(x = V1, y = value, fill = Variables)) + geom_bar(stat="identity", position='stack') +
      facet_wrap(~variable, ncol = 1) +
      xlab("Horizon") + ylab("Contribution to FEV [in %]") + scale_fill_grey() +
      theme_bw() + theme(legend.title=element_blank()))
  }
  return(fe2)
}
