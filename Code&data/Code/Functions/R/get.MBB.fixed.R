#' Moving Block Boostraps with fixed design
#'
#' @param My.id An object of of class "my.id", identified SVAR
#' @param Step Integer, Integer, IRF horizons
#' @param nboot Integer, number of bootstrap replicates.
#' @param Ci Numeric, quantile of boostraped IRFs
#' @param b.length Integer, length of each block
#' @param design character. If design="fixed", a fixed design bootstrap is performed. If design="recursive", a recursive design bootstrap is performed.
#' @param Integrate logical, If Integrate = TRUE, the IRFs will be integrated
#'
#' @return
#' @export get.MBB.fixed
#'
#' @references Jentsch, Carsten, and Kurt G. Lunsford. 2019. "The Dynamic Effects of Personal and Corporate Income Tax Changes in the United States: Comment." American Economic Review, 109 (7): 2655-78.
#'
#' @examples NULL
get.MBB.fixed <- function(My.id, design = "fixed", b.length = 15, Step, nboot, Ci, Integrate = F){
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
  y       <- My.id$dat
  A       <- t(Beta)

  if (Method != "iv"){
    stop("For other id.methods, please use mb.boot from svars package!")
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


    USTAR <- list()

    # creating blocks
    N <- Tob/b.length
    U_blocks    <- array(NA, c(b.length, K, Tob - b.length + 1))
    IV_blocks   <-  array(NA, c(b.length, Tob - b.length + 1))
    #time_blocks <-

    for (i in 0:(Tob - b.length)) {
      U_blocks[, , (i + 1)] <- as.matrix(syn.ts[(i + 1):(i + b.length),Varname])
      IV_blocks[,(i + 1)] <- as.matrix(syn.ts[(i + 1):(i + b.length),(K+2)])
      #time_blocks[,(i + 1)] <- as.matrix(syn.ts[(i + 1):(i + b.length),1])
    }


    for (i in 1:nboot) {

      epsilon.star <- list()
      iv.star <- matrix(0, nrow = (ceiling(N) + 1)*b.length, ncol = 1)
      #time.star <-
      # stacking randomly selected blocks at each other
      for (kk in 1:(ceiling(N) + 1)) {
        eta <- floor(runif(1, 1, Tob - b.length + 2))
        epsilon.star[[kk]] <- U_blocks[, , eta]
        iv.star[(b.length*(kk-1)+1):(b.length*kk),]     <- IV_blocks[, eta]
        #time.star[(b.length*(kk-1)+1):(b.length*kk),]   <- time_blocks[, eta]
      }
      epsilon.star  <- do.call('rbind', epsilon.star)

      # centering
      for(s in 1:b.length){
        b.mean <- colSums(epsilon.star[1 : (s + (Tob - b.length)),])/(Tob - b.length + 1)
        z.mean <- sum(iv.star[1 : (s + (Tob - b.length)),], na.rm = T)/(Tob - b.length + 1)
        for(j in 0:floor(N)){
          epsilon.star[j * b.length + s,] <- epsilon.star[j * b.length + s,] - b.mean
          iv.star[j * b.length + s,] <- iv.star[j * b.length + s,] - z.mean
        }
      }


      # cutting of unnecessary observations
      epsilon.star <- epsilon.star[1:Tob, ]
      iv.star     <- iv.star[1:Tob, ]

      USTAR[[i]] <- list("ustar" = epsilon.star,"ivstar" = iv.star)
    }
  } else {
    stop("Under development!")
  }

  # kernel function
  kernel.boots.fun <- function(ustars){

    if (design == "recursive"){

      stop("Under development!")


      Ystar <- matrix(0, nrow(y), K)
      # adding pre sample values
      Ystar[1:p,] <- y[1:p,]

      if (My.id$type == 'const' | My.id$type == 'trend') {
        for (i in (p + 1):nrow(y)) {
          for (j in 1:K) {
            Ystar[i, j] <- A[j, 1] + A[j, -1] %*% c(t(Ystar[(i - 1):(i - p), ])) + ustars$ustar[(i - p),j]
          }
        }
      } else {
        stop("Under development!")
      }

      # Delete pre sample values
      #Ystar <- Ystar[-c(1:p), ]

      varb <- suppressWarnings(VAR(Ystar, p = My.id$p, type = My.id$type))
      Beta.B <- t(Bcoef(varb))
      U.star <- residuals(varb)

    } else if (design == "fixed"){
      Y.star <- Z %*% Beta + ustars$ustar
      Beta.B <- solve(t(Z) %*% Z) %*% t(Z) %*% as.matrix(Y.star)
      U.star <- as.matrix(Y.star) - Z %*% Beta.B
    }
    Sigma_u <- crossprod(U.star)/(Tob - K * p - 1)
    #colnames(Sigma_u) <- row.names(Sigma_u) <- Varnames

    #var.temp <- vars::VAR(Y.star, p = p, type = Type)

    if (Method == "iv"){

      iv <- iv <- ifelse(is.na(ustars$ivstar), 0, ustars$ivstar)

      Pi          <- solve(crossprod(U.star)/Tob) %*% (crossprod(U.star, iv)/Tob)
      phi2        <- (crossprod(iv, U.star)/Tob) %*% solve(crossprod(U.star)/Tob) %*% (crossprod(U.star, iv)/Tob)
      B_k         <- (crossprod(U.star, iv)/Tob) %*% (1/sqrt(phi2))
    } else {
      stop("For other id.methods, please use wild.boot from svars package!")
    }

    Theta.temp <- get.my.irf(Simple = T, Bcoef = Beta.B, Bmat = B_k, Step = Step, Lag = p, Integrate = Integrate)

    temp.erg <- list("B.boot" = B_k, "irf.B"= Theta.temp, "Sigma_u" = Sigma_u)

    return(temp.erg)
  }
  BOOTS <- pbapply::pblapply(USTAR, FUN = kernel.boots.fun, cl = 1)

  # summerize results
  erg.Bmats <- matrix(NA, nrow = nboot, ncol = K^2)
  erg.irf <- matrix(NA, nrow = nboot, ncol = K^2*(Step + 1))

  arry.Bmats <- array(NA, dim = c(K, K, nboot))
  arry.Sigma <- array(NA, dim = c(K, K, nboot))

  for (i in 1:nboot) {
    erg.Bmats[i,] <- c(BOOTS[[i]][["B.boot"]])
    erg.irf[i,]   <- c(BOOTS[[i]][["irf.B"]])

    arry.Bmats[,,i] = c(BOOTS[[i]][["B.boot"]])
    arry.Sigma[,,i] = c(BOOTS[[i]][["Sigma_u"]])
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
              "IRF" = list("point" = get.my.irf(Simple = T, Bcoef = Beta, Bmat = B, Step = Step, Lag = p, Integrate = Integrate),
                           "mean" = IRF.boot.mean, "L" = IRF.boot.L, "U" = IRF.boot.U))
  erg[["Method"]] <- c("MBB", "Fixed", Method)
  erg[["boot_B"]] <- arry.Bmats
  erg[["boot_Sigma"]] <- arry.Sigma
  erg[["nBoots"]] <- nboot
  class(erg) <- "my.id.boot"

  return(erg)
}






