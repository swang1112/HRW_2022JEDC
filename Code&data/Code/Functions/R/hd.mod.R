#' modifited historical decomposition
#'
#' @param x svar or my.id
#' @param series which series to be decomposed
#' @param Partial partial identification, by default Partial = NULL
#' @param Epsname Vector of shock names (string)
#'
#' @return data.frame
#' @export hd.mod
#'
#' @examples NULL
hd.mod <- function(x, series = 1, Partial = NULL, Epsname = NULL){
  # Function to calculate matrix potence
  "%^%" <- function(A, n){
    if(n == 1){
      A
    }else{
      A %*% (A %^% (n-1))
    }
  }


  # function to calculate impulse response
  IrF <- function(A_hat, B_hat, horizon){
    k <- nrow(A_hat)
    p <- ncol(A_hat)/k
    if(p == 1){
      irfa <- array(0, c(k, k, horizon))
      irfa[,,1] <- B_hat
      for(i in 1:horizon){
        irfa[,,i] <- (A_hat%^%i)%*%B_hat
      }
      return(irfa)
    }else{
      irfa <- array(0, c(k, k, horizon))
      irfa[,,1] <- B_hat
      Mm <- matrix(0, nrow = k*p, ncol = k*p)
      Mm[1:k, 1:(k*p)] <- A_hat
      Mm[(k+1):(k*p), 1 : ((p-1)*k)] <- diag(k*(p-1))
      Mm1 <- diag(k*p)
      for(i in 1:(horizon-1)){
        Mm1 <- Mm1%*%Mm
        irfa[,,(i+1)] <- Mm1[1:k, 1:k]%*%B_hat
      }
      return(irfa)
    }
  }

  if (class(x) == "my.id"){
    if (x$method[1] == "sr"){
      Varname <- x$Varname
      k <- length(x$Varname)
      A_hat <- t(x$Beta)[,  1 : (x$p * k)]

      B_hat <- x$IRF.MT[,,1] #median target
      horizon <- x$Tob

      # MA coeff
      IR <- IrF(A_hat, B_hat, horizon)
      if(is.null(Epsname)){Epsname <- Varname}
      impulse <- IR[series,,]

      ## struct shocks
      p <- x$p
      s.time <- time(x$dat)[-(1:p)]
      y <- x$dat[-c(1:p), series]

      if (is.null(Partial)){
        s.errors <- x$epsilon.MT
        y_hat <- matrix(NA, nrow = horizon, ncol = k)
        for (i in 1:horizon) {
          for (j in 1:k) {
            y_hat[i,j] <-  t(impulse[j,1:i]) %*% s.errors[j,i:1]
          }
        }

        y_hat_a <- rowSums(y_hat)

        yhat <- data.frame(s.time, (y - mean(y)), y_hat_a, y_hat)
        colnames(yhat)[1:3]<- c("t",
                                paste("Demeaned series ", Varname[series]),
                                paste("Constructed series ", Varname[series]))
        for(i in 4:ncol(yhat)){
          colnames(yhat)[i] <- paste("Cumulative effect of ", Epsname[i-3], "shock on ", Varname[series])
        }

      } else {
        s.errors <- x$epsilon.MT[Partial,]
        y_hat <- rep(NA, horizon)
        for (i in 1:horizon) {
          y_hat[i] <- t(impulse[Partial,1:i]) %*% s.errors[i:1]
        }
        yhat <- data.frame(s.time, y_hat)
        colnames(yhat) <- c("t",
                            paste("Cumulative effect of ", Epsname[Partial], "shock on ", Varname[series]))

      }


    } else if(x$method[1] == "sr.b"){
      Varname <- x$Varname
      k <- length(x$Varname)
      A_hat <- x$Beta

      B_hat <- x$B
      horizon <- x$Tob

      # MA coeff
      IR <- IrF(A_hat, B_hat, horizon)
      if(is.null(Epsname)){Epsname <- Varname}
      impulse <- IR[series,,]

      ## struct shocks
      p <- x$p
      s.time <- time(x$dat)[-(1:p)]
      y <- x$dat[-c(1:p), series]

      if (is.null(Partial)){
        s.errors <- x$epsilon.M
        y_hat <- matrix(NA, nrow = horizon, ncol = k)
        for (i in 1:horizon) {
          for (j in 1:k) {
            y_hat[i,j] <-  t(impulse[j,1:i]) %*% s.errors[j,i:1]
          }
        }

        y_hat_a <- rowSums(y_hat)

        yhat <- data.frame(s.time, (y - mean(y)), y_hat_a, y_hat)
        colnames(yhat)[1:3]<- c("t",
                                paste("Demeaned series ", Varname[series]),
                                paste("Constructed series ", Varname[series]))
        for(i in 4:ncol(yhat)){
          colnames(yhat)[i] <- paste("Cumulative effect of ", Epsname[i-3], "shock on ", Varname[series])
        }

      } else {
        s.errors <- x$epsilon.M[Partial,]
        y_hat <- rep(NA, horizon)
        for (i in 1:horizon) {
          y_hat[i] <- t(impulse[Partial,1:i]) %*% s.errors[i:1]
        }
        yhat <- data.frame(s.time, y_hat)
        colnames(yhat) <- c("t",
                            paste("Cumulative effect of ", Epsname[Partial], "shock on ", Varname[series]))

      }


    }  else if (x$method[1] == "iv"){


      Varname <- x$Varname
      k <- length(x$Varname)
      A_hat <- t(x$Beta)[,  1 : (x$p * k)]

      B_hat <- matrix(x$B, nrow = k, ncol = k)
      horizon <- length(x$epsilon)

      # MA coeff
      IR <- IrF(A_hat, B_hat, horizon)
      if(is.null(Epsname)){Epsname <- Varname}
      impulse <- IR[series,,]

      ## struct shocks
      p <- x$p
      s.time <- x$eps.ts$time
      y <-  ts2df(x$dat, Varname) %>%
        dplyr::filter(time >= s.time[1] & time <= s.time[length(s.time)]) %>%
        dplyr::select(-1)

      s.errors <- x$epsilon
      y_hat <- rep(NA, horizon)

      if (is.null(Partial)){
        stop("For proxy SVAR identification: only partial identification allowed!")
      }

      for (i in 1:horizon) {
        y_hat[i] <- t(impulse[Partial,1:i]) %*% s.errors[i:1]
      }
      yhat <- data.frame(s.time, y_hat)
      colnames(yhat) <- c("t",
                          paste("Cumulative effect of ", Epsname[Partial], "shock on ", Varname[series]))



    }

  } else if (class(x) == "svars"){
    Varname <- colnames(x$y)
    k <- x$K
    p <- x$p
    horizon <- x$n
    B_hat <- x$B

    if(x$type == "const"){
      A_hat <- x$A_hat[,-1]
    }else if(x$type == "trend"){
      A_hat <- x$A_hat[,-1]
    }else if(x$type == "both"){
      A_hat <- x$A_hat[,-c(1,2)]
    }else{
      A_hat <- x$A_hat
    }

    IR <- IrF(A_hat, B_hat, horizon)

    if(is.null(Epsname)){Epsname <- Varname}
    impulse <- IR[series,,]

    ## struct shocks

    p <- x$p
    s.time <- time(x$y)[-(1:p)]
    y <- x$y[-c(1:p), series]

    u <- t(resid(x$VAR))

    s.errors <- solve(B_hat)%*%u

    if (is.null(Partial)){
      y_hat <- matrix(NA, nrow = horizon, ncol = k)
      for (i in 1:horizon) {
        for (j in 1:k) {
          y_hat[i,j] <-  t(impulse[j,1:i]) %*% s.errors[j,i:1]
        }
      }

      y_hat_a <- rowSums(y_hat)

      yhat <- data.frame(s.time, (y - mean(y)), y_hat_a, y_hat)
      colnames(yhat)[1:3]<- c("t",
                              paste("Demeaned series ", Varname[series]),
                              paste("Constructed series ", Varname[series]))
      for(i in 4:ncol(yhat)){
        colnames(yhat)[i] <- paste("Cumulative effect of ", Epsname[i-3], "shock on ", Varname[series])
      }

    } else {
      s.errors <- s.errors[Partial,]
      y_hat <- rep(NA, horizon)
      for (i in 1:horizon) {
        y_hat[i] <- t(impulse[Partial,1:i]) %*% s.errors[i:1]
      }
      yhat <- data.frame(s.time, y_hat)
      colnames(yhat) <- c("t",
                          paste("Cumulative effect of ", Epsname[Partial], "shock on ", Varname[series]))
    }


  }

  return(yhat)



}
