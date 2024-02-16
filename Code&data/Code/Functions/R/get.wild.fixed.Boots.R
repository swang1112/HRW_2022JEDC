#' Wild bootstrap with fixed design
#'
#' @param My.id An object of of class "my.id", identified SVAR
#' @param Dist Distribution, from which independent random noises are drawn.
#' @param nboot Integer, number of bootstrap replicates.
#' @param Step Integer, Integer, IRF horizons
#' @param Ci Numeric, quantile of boostraped IRFs
#'
#' @return
#' @export get.wild.fixed.Boots
#'
#' @examples
get.wild.fixed.Boots <- function(My.id, Dist = "rademacher", Step, nboot, Ci){
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

  if (Method != "iv"){
    stop("For other id.methods, please use wild.boot from svars package!")
  }

  cat("\r", "...calculating finish time...")

  u <- My.id$u - colMeans(My.id$u)
  instruments <- My.id$Instrument


  if (My.id$Synchro){
    u_time  <- time(My.id$dat)
    u_time  <- u_time[-(1:p)]
    iv_time <- time(instruments)
    u.ts    <- data.frame("time" = u_time, "u" = u)
    iv.ts   <- data.frame("time" = iv_time, "iv" = instruments)
    syn.ts  <- suppressWarnings(dplyr::left_join(u.ts, iv.ts, by = "time"))

    colnames(syn.ts) <- c("time", Varname, "iv")
    time <- seq(min(u_time), max(u_time, iv_time), by = 1/frequency(iv_time))

    USTAR <- list()
    for (i in 1:nboot) {
      if (Dist == "rademacher"){
        etas <- rnorm(n = length(time))
        etas <- (etas > 0) - (etas < 0)
      } else if (Dist == "mammen"){
        cu <- (sqrt(5)+1)/(2*sqrt(5))
        etas <- rep(1,length(time))*(-(sqrt(5)-1)/2)
        uni <- runif(n = length(time), min = 0, max = 1)
        etas[uni > cu] <- (sqrt(5)+1)/2
      } else if(Dist == "gaussian"){
        etas <- rnorm(n = length(time))
      } else {
        stop("Undefined distribution! Dist = c('rademacher', 'mammen', 'gaussian').")
      }

      eta   <- cbind("time" = time, "eta"=etas) %>% as.data.frame()
      temp  <- dplyr::left_join(syn.ts, eta, by = "time") %>%
        mutate_at(.funs = list(~.*eta), .vars = c(Varname, "iv"))

      USTAR[[i]] <- list("ustar"  = temp %>% dplyr::select(Varname) %>% tidyr::drop_na(),
                         "ivstar" = temp %>% dplyr::select(c("time", "iv")) %>% tidyr::drop_na())
    }
  } else {
    USTAR <- list()
    for (i in 1:nboot) {
      if (Dist == "rademacher"){
        etas <- rnorm(n = Tob)
        etas <- (etas > 0) - (etas < 0)
      } else if (Dist == "mammen"){
        cu <- (sqrt(5)+1)/(2*sqrt(5))
        etas <- rep(1,Tob)*(-(sqrt(5)-1)/2)
        uni <- runif(n = Tob, min = 0, max = 1)
        etas[uni > cu] <- (sqrt(5)+1)/2
      } else if(Dist == "gaussian"){
        etas <- rnorm(n = Tob)
      } else {
        stop("Undefined distribution! Dist = c('rademacher', 'mammen', 'gaussian').")
      }
      USTAR[[i]] <- list("ustar" = u * etas, "ivstar" = instruments * etas)
    }
  }

  # kernel function
  kernel.boots.fun <- function(ustars){

      Y.star <- Z %*% Beta + ustars$ustar
      Beta.B <- solve(t(Z) %*% Z) %*% t(Z) %*% as.matrix(Y.star)
      U.star <- as.matrix(Y.star) - Z %*% Beta.B
      #Sigma_u <- crossprod(U.star)/(Tob - K * p - 1)
      #colnames(Sigma_u) <- row.names(Sigma_u) <- Varnames

      #var.temp <- vars::VAR(Y.star, p = p, type = Type)

    if (Method == "iv"){

      if (My.id$Synchro){
        U.star.ts <- data.frame("time" = u_time, "u" = U.star)
        Syn.star.ts <- suppressWarnings(dplyr::left_join(U.star.ts, ustars$ivstar, by = "time") %>% tidyr::drop_na())

        iv <- as.matrix(Syn.star.ts$iv)
        U.star <- as.matrix(Syn.star.ts[,c(2:(K+1))])
        Tob <- nrow(U.star)
      }

      Pi          <- solve(crossprod(U.star)/Tob) %*% (crossprod(U.star, iv)/Tob)
      phi2        <- (crossprod(iv, U.star)/Tob) %*% solve(crossprod(U.star)/Tob) %*% (crossprod(U.star, iv)/Tob)
      B_k         <- (crossprod(U.star, iv)/Tob) %*% (1/sqrt(phi2))
    } else {
      stop("For other id.methods, please use wild.boot from svars package!")
    }

    Theta.temp <- get.my.irf(Simple = T, Bcoef = Beta.B, Bmat = B_k, Step = Step, Lag = p)

    temp.erg <- list("B.boot" = B_k, "irf.B"= Theta.temp)

    return(temp.erg)

  }

  BOOTS <- pbapply::pblapply(USTAR, FUN = kernel.boots.fun, cl = 1)

  # summerize results
  erg.Bmats <- matrix(NA, nrow = nboot, ncol = K^2)
  erg.irf <- matrix(NA, nrow = nboot, ncol = K^2*(Step + 1))

  for (i in 1:nboot) {
    erg.Bmats[i,] <- c(BOOTS[[i]][["B.boot"]])
    erg.irf[i,]   <- c(BOOTS[[i]][["irf.B"]])
  }

  # Bmats
  B.boot.mean <- array(apply(erg.Bmats, 2, mean), dim = c(K, K))
  B.boot.sd   <- array(apply(erg.Bmats, 2, sd), dim = c(K, K))
  # IRF
  IRF.boot.mean <- array(apply(erg.irf, 2, mean), dim = c(K, K, (Step + 1)))
  IRF.boot.ci <- apply(erg.irf, 2, quantile, probs = c((1-Ci)/2, (1+Ci)/2))
  IRF.boot.L <- array(IRF.boot.ci[1,], dim = c(K, K, (Step + 1)))
  IRF.boot.U <- array(IRF.boot.ci[2,], dim = c(K, K, (Step + 1)))

  # rename
  colnames(B.boot.mean) <- colnames(B.boot.sd) <- colnames(IRF.boot.mean) <-
    colnames(IRF.boot.L) <- colnames(IRF.boot.U) <- Varname
  row.names(B.boot.mean) <- row.names(B.boot.sd) <- row.names(IRF.boot.mean) <-
    row.names(IRF.boot.L) <- row.names(IRF.boot.U) <- Varname

  erg <- list("Bmat" = list("mean" = B.boot.mean, "sd" = B.boot.sd),
              "IRF" = list("point" = get.my.irf(Simple = T, Bcoef = Beta, Bmat = B, Step = Step, Lag = p),
                           "mean" = IRF.boot.mean, "L" = IRF.boot.L, "U" = IRF.boot.U))
  erg[["Method"]] <- c("wild", Dist, "Fixed", Method)
  erg[["nBoots"]] <- nboot
  class(erg) <- "my.id.boot"

  return(erg)
}

#' Check sboot sign
#'
#' @param x object "sboot"
#'
#' @return
#' @export signcheck.sboot
#'
#' @examples NULL
signcheck.sboot <- function(x){
  irfs <- lapply(x$bootstrap, '[[', 'irf')
  B0 <- irfs[[1]] %>% filter(V1 == 1)
  K <- sqrt(ncol(B0)-1)
  for (i in 2:length(irfs)) {
    B0 <- rbind(B0, irfs[[i]] %>% filter(V1 == 1))
  }

  erg <- list("negative"= apply(B0 %>% select(-1), 2, function(x){mean(x<0)})  %>% matrix(nrow = K, ncol = K, byrow = T),
              "positive"= apply(B0 %>% select(-1), 2, function(x){mean(x>0)})  %>% matrix(nrow = K, ncol = K, byrow = T))
  return(erg)
}

