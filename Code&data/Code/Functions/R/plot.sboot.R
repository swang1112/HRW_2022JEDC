
#' modified plot svar.boot
#'
#' @param x sboot object
#' @param scales NA
#' @param lowerq NA
#' @param upperq NA
#' @param Partial Parial identification
#' @param percentile 'standard', 'hall' , 'bonferroni'
#' @param ...
#' @param base NA
#' @param Normalize Numeric, normalization of IRFs, IRFs(new) = IRFs / Normalize
#'
#' @return irf plot
#' @export plot.sboot.mod
#'
#' @examples NA
plot.sboot.mod <- function(x, scales = "free_y", lowerq = 0.16, upperq = 0.84, Partial = NULL, Normalize = NULL,  percentile = 'standard', ..., base){

  # define
  probs <- NULL
  V1 <- NULL
  value <- NULL

  n.ahead <- nrow(x$true$irf)
  kk <- ncol(x$true$irf)
  K <- sqrt(kk-1)
  bootstrap <- x$bootstrap
  nboot <- length(bootstrap)
  rest <- x$rest_mat

  n.probs <- length(lowerq)
  if(length(lowerq) != length(upperq)){
    stop("Vectors 'lowerq' and 'upperq' must be of same length!")
  }

  intervals <- array(0, c(n.ahead, kk, nboot))
  for(i in 1:nboot){
    intervals[,,i] <- as.matrix(bootstrap[[i]]$irf)
  }

  # find quantiles for lower and upper bounds
  lower <- array(0, dim = c(n.ahead, kk, n.probs))
  upper <- array(0, dim = c(n.ahead, kk, n.probs))
  if(percentile == 'standard' | percentile == 'hall'){
    for(i in 1:n.ahead){
      for(j in 1:kk){
        lower[i,j, ] <- quantile(intervals[i,j, ], probs = lowerq)
        upper[i,j, ] <- quantile(intervals[i,j, ], probs = upperq)
      }
    }

    if(percentile == 'hall'){
      for(p in 1:n.probs){
        lower[ ,-1, p] <- as.matrix(2*x$true$irf)[ ,-1] - lower[ ,-1, p]
        upper[ ,-1, p] <- as.matrix(2*x$true$irf)[ ,-1] - upper[ ,-1, p]
      }
    }

  }else if(percentile == 'bonferroni'){
    rest <- matrix(t(rest), nrow = 1)
    rest[is.na(rest)] <- 1
    rest <- c(1, rest)
    for(i in 1:n.ahead){
      for(j in 1:kk){
        if(rest[j] == 0){
          lower[i,j, ] <- quantile(intervals[i,j, ], probs = (lowerq/n.ahead))
          upper[i,j, ] <- quantile(intervals[i,j, ], probs = 1 + (((upperq - 1)/n.ahead)))
        }else{
          lower[i,j, ] <- quantile(intervals[i,j, ], probs = (lowerq/(n.ahead + 1)))
          upper[i,j, ] <- quantile(intervals[i,j, ], probs = 1 + (((upperq - 1)/(n.ahead + 1))))
        }
      }
    }
  }else{
    stop("Invalid choice of percentile; choose between standard, hall and bonferroni")
  }

  # plot IRF with confidence bands
  alp <- 0.7 * (1+log(n.probs, 10))/n.probs
  irf <- reshape2::melt(x$true$irf, id = 'V1')
  cbs <- data.frame(V1 = rep(irf$V1, times=n.probs),
                    variable = rep(irf$variable, times=n.probs),
                    probs = rep(1:n.probs, each=(kk-1)*n.ahead),
                    lower = c(lower[,-1,]),
                    upper = c(upper[,-1,]))
  Response <- dplyr::left_join(cbs, irf, by = c("V1", "variable"))

  if(!is.null(Normalize)){
    Response %<>% mutate(value = value / Normalize, upper = upper / Normalize, lower = lower  / Normalize)
  }

  # Partial identification
  plot.col <- K
  if(!is.null(Partial)){
    Response.partial <- data.frame()

    for (k in 1:K) {
      Response.temp <- Response[((n.ahead)*K*(k-1) + 1): ((n.ahead)*K*k),]
      Response.partial <- dplyr::bind_rows(Response.partial,Response.temp[((n.ahead)*(Partial-1)+1):((n.ahead)*Partial),])
    }
    Response <- Response.partial
    plot.col <- 1
  }

  erg <- ggplot(data=Response) +
    geom_ribbon(data=Response, aes(x=V1, ymin=lower, ymax=upper, group=probs), alpha=alp, fill='darkgrey') +
    geom_line(data=Response, aes(x=V1, y=value)) +
    geom_hline(yintercept = 0, color = 'red') +
    facet_wrap(~variable, scales = scales, labeller = label_parsed,  ncol = plot.col) +
    xlab("Horizon") + ylab("Response") +
    theme_bw()

  return(erg)
}
