#' Modified Historical Counterfactuals for SVAR Models
#'
#' Calculation of historical decomposition (for partial identification) for a proxy identified SVAR object 'my.id'
#'
#' @param x my.id
#' @param series Integer. indicating the series that should be decomposed.
#' @param transition Numeric. Value from [0, 1] indicating how many initial values should be discarded, i.e., 0.1 means that the first 10 per cent observations of the sample are considered as transient.
#' @param Epsname Characters, Name the structural shocks
#' @param Partial Integer. indicating which shock is identified
#'
#' @examples NULL
#' @export cf.mod

cf.mod <- function(x, series = 1, transition = 0, Partial = NULL, Epsname = NULL){

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


  if(class(x) == "my.id"){

    ## Step 1: Calculate MA coefficients

    k <- length(x$Varname)
    A_hat <- t(x$Beta)[,  1 : (x$p * k)]

  if (x$method[1] == "sr"){

    B_hat <- x$IRF.MT[,,1] #median target
    horizon <- x$Tob

    IR <- IrF(A_hat, B_hat, horizon)
    if(is.null(Epsname)){Epsname <- x$Varname}

    impulse <- matrix(0, ncol = dim(IR)[2]^2 + 1, nrow = dim(IR)[3])
    colnames(impulse) <- rep("V1", ncol(impulse))
    cc <- 1
    impulse[,1] <- seq(1, dim(IR)[3])
    for(i in 1:dim(IR)[2]){
      for(j in 1:dim(IR)[2]){
        cc <- cc + 1
        impulse[,cc] <- IR[i,j,]
        colnames(impulse)[cc] <- paste("epsilon[",Epsname[j],"]", "%->%", x$Varname[i])
      }
    }


    # Step 2: Calculate structural errors
    p <- x$p
    s.errors <- x$epsilon.M[Partial,]
    s.time <- time(x$dat)[-(1:p)]

    # Step 3: Match up structural shocks with appropriate impuslse response
    impulse <- impulse[,-1]
    y_hat <- rep(NA, length(s.errors))

    if (transition == 0) {
      for (i in 1:length(s.errors)) {

        y_hat[i] <- impulse[1:i, 1 + series] %*% t(s.errors)[i:1]

      }
    } else {
      for (i in (length(s.errors) - floor(length(s.errors) * (1 - transition))):length(s.errors)) {
        y_hat[i] <- impulse[(length(s.errors) - floor(length(s.errors) * (1 - transition))):i, 1 + series] %*% t(s.errors)[i:(length(s.errors) - floor(length(s.errors) * (1 - transition)))]

      }
      s.time <- s.time[-transition]
    }

    #y_hat <- ts(y_hat, frequency = frequency(x$dat), start = c(floor(x$eps.ts[1,1]), (x$eps.ts[1,1] %% 1) / 0.25 + 1 ))


    y_hat <- data.frame("time" = s.time, y_hat)

    y_true <- ts2df(x$dat[,series], Varname = "y_true")
    y_true[,2] <- y_true[,2] - mean(y_true[,2])

    yhat <- left_join(y_hat,  y_true, by = "time")

    y_cf <- yhat[,3] - yhat$y_hat

    y_all <- cbind(yhat, y_cf) #%>% select(-y_hat)


  } else if (x$method[1] == "iv"){



    B_hat <- matrix(x$B, k, k)

    horizon <- x$Tob

    IR <- IrF(A_hat, B_hat, horizon)

    if(is.null(Epsname)){Epsname <- x$Varname}

    impulse <- matrix(0, ncol = dim(IR)[2]^2 + 1, nrow = dim(IR)[3])
    colnames(impulse) <- rep("V1", ncol(impulse))
    cc <- 1
    impulse[,1] <- seq(1, dim(IR)[3])
    for(i in 1:dim(IR)[2]){
      for(j in 1:dim(IR)[2]){
        cc <- cc + 1
        impulse[,cc] <- IR[i,j,]
        colnames(impulse)[cc] <- paste("epsilon[",Epsname[j],"]", "%->%", x$Varname[i])
      }
    }

    # Step 2: Calculate structural errors


    s.errors <- t(x$epsilon)
    s.time <- x$eps.ts[,1]
    # Step 3: Match up structural shocks with appropriate impuslse response
    impulse <- impulse[,-1]
    y_hat <- rep(NA, length(s.errors))

    if (transition == 0) {
      for (i in 1:length(s.errors)) {

        y_hat[i] <- impulse[1:i, 1 + series] %*% t(s.errors)[i:1]

      }
    } else {
      for (i in (length(s.errors) - floor(length(s.errors) * (1 - transition))):length(s.errors)) {
        y_hat[i] <- impulse[(length(s.errors) - floor(length(s.errors) * (1 - transition))):i, 1 + series] %*% t(s.errors)[i:(length(s.errors) - floor(length(s.errors) * (1 - transition)))]

      }
      s.time <- s.time[-transition]
    }

    #y_hat <- ts(y_hat, frequency = frequency(x$dat), start = c(floor(x$eps.ts[1,1]), (x$eps.ts[1,1] %% 1) / 0.25 + 1 ))


    y_hat <- data.frame("time" = s.time, y_hat)

    y_true <- ts2df(x$dat[,series], Varname = "y_true")
    y_true[,2] <- y_true[,2] - mean(y_true[,2])

    yhat <- left_join(y_hat,  y_true, by = "time")

    y_cf <- yhat[,3] - yhat$y_hat

    y_all <- cbind(yhat, y_cf) #%>% select(-y_hat)
  }
    } else {
  p <- x$p
  k <- x$K
  obs <- x$VAR$obs

  ## Step 1: Calculate MA coefficients
  if(x$type == "const"){
    A_hat <- x$A_hat[,-1]
  }else if(x$type == "trend"){
    A_hat <- x$A_hat[,-1]
  }else if(x$type == "both"){
    A_hat <- x$A_hat[,-c(1,2)]
  }else{
    A_hat <- x$A_hat
  }


  IR <- IrF(A_hat, x$B, horizon)
  impulse <- matrix(0, ncol = dim(IR)[2]^2 + 1, nrow = dim(IR)[3])
  colnames(impulse) <- rep("V1", ncol(impulse))
  cc <- 1
  impulse[,1] <- seq(1, dim(IR)[3])
  for(i in 1:dim(IR)[2]){
    for(j in 1:dim(IR)[2]){
      cc <- cc + 1
      impulse[,cc] <- IR[i,j,]
      colnames(impulse)[cc] <- paste("epsilon[",Epsname[j],"]", "%->%", x$Varname[i])
    }
  }

  # Step 2: Calculate structural errors
  s.errors <- solve(x$B) %*% t(resid( x$VAR))
  s.time <- time(x$y)[-(1:p)]

  # Step 3: Match up structural shocks with appropriate impuslse response
  impulse <- impulse[,-1]
  y_hat <- rep(NA, length(s.errors))

  y_hat <- matrix(NA, nrow = obs, ncol = k)
  if (transition == 0) {
    for (i in 1:obs) {
      for (j in 1:k) {
        y_hat[i, j] <- impulse[1:i, j + series * k - k] %*% t(s.errors)[i:1, j]
      }
    }
  } else {
    for (i in (obs - floor(obs * (1 - transition))):obs) {
      for (j in 1:k) {
        y_hat[i, j] <- impulse[(obs - floor(obs * (1 - transition))):i, j + series * k - k] %*% t(s.errors)[i:(obs - floor(obs * (1 - transition))), j]
      }
    }
    s.time <- s.time[-transition]
    }

  y_hat <- data.frame("time" = s.time, y_hat)
  if(is.null(Epsname)){Epsname <- colnames( x$VAR$y)}
  colnames(y_hat)[-1] <- Epsname

  y_hat <- data.frame(y_hat, "y_hat" =  apply(y_hat[,-1], 1, sum))

  y_true <- ts2df(x$VAR$y[,series], Varname = "y_true")
  y_true[,2] <- y_true[,2] - mean(y_true[,2])

  yhat <- left_join(y_hat,  y_true, by = "time")

  y_cf <- matrix(NA, nrow = nrow(yhat), ncol = k)
  tempname <- rep("tempname", k)
  for (i in 1:k) {
    y_cf[,i] <- yhat$y_true - yhat[,i+1]
    tempname[i] <- paste0("cf_", Epsname[i])
  }
  colnames(y_cf) <- tempname


  y_all <- cbind(yhat, y_cf) #%>% select(-y_hat)

  }

  return(y_all)
}


ts2df <- function(x, Varname){
  erg <- data.frame(as.numeric(time(x)), x)
  colnames(erg) <- c("time", Varname)
  erg
}
