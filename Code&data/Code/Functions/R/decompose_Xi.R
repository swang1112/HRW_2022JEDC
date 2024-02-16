#' NULL
#'
#' @param x
#' @param series
#' @param Partial
#' @param Epsname
#' @param Start
#' @param End
#' @param Freq
#' @param wrong_shocks to replace epsilon3
#'
#' @return NULL
#' @export decompose_Xi
#'
#' @examples NULL
decompose_Xi <- function(x, series = 1, Partial = NULL, Epsname = NULL, Start, End, Freq = 4, wrong_shocks){

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

    ## Step 1: Calculate MA coefficients

    k <- length(x$Varname)
    A_hat <- t(x$Beta)[,  1 : (x$p * k)]


      B_hat <- matrix(x$B, k, k)

      horizon <- (End - Start)*Freq + 1

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

      colnames(wrong_shocks) <- c("time", "eps")
      s.errors <- wrong_shocks %>% filter(time >= Start & time <= End) %>% select(eps)
      s.impulse <- impulse[,(k*(series-1)+Partial +1)]

      erg <- 0
      for (i in 1:horizon) {
        erg <- erg + s.errors$eps[i]*sum(s.impulse[1:(horizon-(i-1))])
      }




  return(erg)

}

