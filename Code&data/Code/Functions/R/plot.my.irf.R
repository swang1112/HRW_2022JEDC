
#' Plot structural impulse response functions
#'
#' @param Step Calculation steps
#' @param My.id Object of class "my.id" or "my.id.boot", identified SVAR
#' @param Partial Integer, indexing the partially identified shocks. NULL by default.
#' @param Hide.Legend Logical, if Hide.Legend = TRUE, legend is removed
#' @param Color.manual Charactuer, set colors manually
#' @param Lty.manual Charactuer, set linetype manually
#' @param Epsname Characters, Name the structural shocks
#' @param Normalize Numeric, normalization of IRFs, IRFs(new) = IRFs / Normalize
#'
#' @return Creat a ggplot object of class "gg"-"ggplot"
#' @importFrom reshape2 melt
#' @importFrom dplyr left_join mutate_all
#' @importFrom urca ca.jo
#' @import ggplot2
#' @export plot.my.irf
#'
#' @examples NULL
plot.my.irf <- function(My.id, Step = NULL, Partial = NULL, Hide.Legend = FALSE, Color.manual = NULL, Lty.manual = NULL, Epsname = NULL, Normalize = NULL){

  myclass <- c("my.id", "my.id.boot")
  if(is(My.id, myclass) == FALSE){
    stop("Input should be an identified SVAR model of class 'my.id' or 'my.id.boot'!")
  }


  if(class(My.id) == "my.id"){
    Varname <- My.id$Varname
    Method <- My.id$method
    K <- length(Varname)
    p <- My.id$p

    # Sign restrictions
    if (Method[1] == "sr"){
      Response <- My.id$IRF.plot

      Step <- dim(My.id$IRF.M)[3] - 1
      # Partial identification
      plot.col <- K
      if(!is.null(Partial)){
        Response.partial <- data.frame()
        for (k in 1:(length(Method)*K)) {
          Response.temp <- Response[((Step+1)*K*(k-1) + 1): ((Step+1)*K*k),]
          Response.partial <- bind_rows(Response.partial,Response.temp[((Step+1)*(Partial-1)+1):((Step+1)*Partial),])
        }
        Response <- Response.partial
        plot.col <- 1
      }

      erg <- ggplot(data = Response, aes(x = h, y = value, group = Label)) +
        geom_line(aes(linetype = Label, color = Label), alpha = 0.9) +
        geom_hline(yintercept = 0, color = 'red') +
        facet_wrap(~variable, scales = "free_y", labeller = label_parsed, ncol = plot.col) +
        #geom_line(aes(y = L), linetype = "dashed", alpha = 0.9) +
        #geom_line(aes(y = U), linetype = "dashed", alpha = 0.9) +
        geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.2) +
        xlab("Horizont") + ylab("Response") +
        theme_bw()

      if (Hide.Legend){
        erg <- erg + theme(legend.position = "none")
      }

      if (!is.null(Color.manual)){
        erg <- erg + scale_color_manual(values = Color.manual)
      }

      if (!is.null(Lty.manual)){
        erg <- erg + scale_linetype_manual(values = Lty.manual)
      }

      return(erg)
    }

    Bcoef <- My.id$Beta
    B <- My.id$B

    Theta <- get.my.irf(Simple = T, Bcoef = Bcoef, Bmat = B, Step = Step, Lag = p)

    Response <- matrix(0, ncol = K^2 + 1, nrow = Step + 1)
    colnames(Response) <- rep("h", K^2 + 1)
    c <- 1
    Response[,c] <- 0:Step
    for(i in 1:K){
      for(j in 1:K){
        c <- c + 1
        Response[,c] <- Theta[i,j,]
        colnames(Response)[c] <- paste("epsilon[", Varname[j],"]", "%->%", Varname[i])
      }
    }
    Response <- melt(as.data.frame(Response), id = 'h')

    # Partial identification
    plot.col <- K
    if(!is.null(Partial)){
      Response.partial <- data.frame()
      for (k in 1:K) {
        Response.temp <- Response[((Step+1)*K*(k-1) + 1): ((Step+1)*K*k),]
        Response.partial <- bind_rows(Response.partial,Response.temp[((Step+1)*(Partial-1)+1):((Step+1)*Partial),])
      }
      Response <- Response.partial
      plot.col <- 1
    }

      erg <- ggplot(data = Response, aes(x = h, y = value)) + geom_line() + geom_hline(yintercept = 0, color = 'red') +
        facet_wrap(~variable, scales = "free_y", labeller = label_parsed, ncol = plot.col) +
        xlab("Horizon") + ylab("Response") +
        theme_bw()

  } else if (class(My.id) == "my.id.boot") {
    irf.start <- My.id$IRF$point
    Varname   <- row.names(irf.start)
    if(is.null(Epsname)){Epsname <- Varname}
    K         <- dim(irf.start)[1]
    Step      <- dim(irf.start)[3] - 1

    Response           <- matrix(0, ncol = K^2 + 1, nrow = Step + 1)
    colnames(Response) <- rep("h", K^2 + 1)
    c                  <- 1
    Response[,c]       <- 0:Step

    # Bootstrap CI
    Response_L <- Response_U <- Response

    for(i in 1:K){
      for(j in 1:K){
        c <- c + 1
        Response[,c] <- irf.start[i,j,]

        Response_L[,c] <- My.id$IRF$L[i,j,]
        Response_U[,c] <- My.id$IRF$U[i,j,]
        colnames(Response_L)[c] <- colnames(Response_U)[c] <- colnames(Response)[c] <-
          paste("epsilon[", Epsname[j],"]", "%->%", Varname[i])
      }
    }
    Response_point <- melt(as.data.frame(Response), id = 'h')

    Response_L <- melt(as.data.frame(Response_L), id = 'h', value.name = "L")
    Response_U <- melt(as.data.frame(Response_U), id = 'h', value.name = "U")
    Response_CI   <- left_join(Response_L, Response_U, by = c("variable", "h"))
    Response <- left_join(Response_point, Response_CI, by = c("variable", "h"))

    if(!is.null(Normalize)){
      Response %<>% mutate(value = value / Normalize, L = L / Normalize, U = U  / Normalize)
    }

    # Partial identification
    plot.col <- K
    if(!is.null(Partial)){
      Response.partial <- data.frame()
      for (k in 1:K) {
        Response.temp <- Response[((Step+1)*K*(k-1) + 1): ((Step+1)*K*k),]
        Response.partial <- bind_rows(Response.partial,Response.temp[((Step+1)*(Partial-1)+1):((Step+1)*Partial),])
      }
      Response <- Response.partial
      plot.col <- 1
    }

    erg <- ggplot(data = Response, aes(x = h, y = value)) + geom_line() + geom_hline(yintercept = 0, color = 'red') +
      facet_wrap(~variable, scales = "free_y", labeller = label_parsed, ncol = plot.col) +
      geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
      xlab("Horizon") + ylab("Response") +
      theme_bw()

  }

  return(erg)
}
