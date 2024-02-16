#' Identification based on Sign restrictions
#'
#' @param Model A class of 'varest' or 'vec2var'
#' @param Iter Integer, number of successful draws
#' @param Signmat Matrix, a (KxK) matrix containing Sign restrictions. 1: positive sign; -1: negative sign; NA: unrestricted
#' @param RHorizont Integer, RHorizonts of impulse responses, for which the sign restrictions are imposed
#' @param PositiveDiag Logical, if PositiveDiag = TRUE, the diagonal of impact matrix will be forced to be positive
#' @param Step Integer, IRF horizont
#' @param Ci Double, between [0,1] sample quantiles for IRFs
#' @param Plot Logical, if Plot = TRUE, IRF plot will be generated
#' @param MT Logical, if MT = TRUE, Fry and Pagan’s (2011) median target (MT) method will be used
#' @param Itermax Integer, number of maximum draws
#' @param Hist Logical, if Hist = TRUE, histogram of distribution of unrestricted elements in the identified B matrices is provided
#' @param Normalize logical, if Normalize = True, size of the shocks will be normalized
#' @param Integrate logical, If Integrate = TRUE, the IRFs will be integrated
#' @param Epsname Characters, Name the structural shocks
#'
#' @return A list containing a (KxKxm) array (m-sets) of decomposition matrices which fulfill the sign restrictions, arrays of
#' impulse responses (median, lower and upper bound)
#' @importFrom reshape2 melt
#' @importFrom dplyr left_join
#' @export get.id.Sign
#'
#' @examples NULL
get.id.Sign <- function(Model, Iter, Itermax, Signmat, RHorizont, MT = FALSE, Integrate = FALSE, PositiveDiag = TRUE, Normalize = TRUE,
                        Step, Ci, Plot = FALSE, Hist = FALSE, Epsname = NULL){
  p       <- Model$p
  K       <- Model$K
  u       <- resid(Model)
  Tob     <- Model$obs
  Beta    <- Model %>% vars::Bcoef() %>% t() # (Kp + mu) x K parameter matrix
  Z       <- Model$datamat[,-c(1:K)] %>% as.matrix() # (T-p) x (Kp + mu) covariates matrix
  Varname <- names(Model$varresult)
  if(is.null(Epsname)){Epsname <- Varname}
  Covmat  <- crossprod(u) / (Tob - 1)
  B_chol  <- t(chol(Covmat))
  M       <- choose(K, 2)
  i <- im <- 1
  B.win <- list()
  # check btw. full sign restriction or argnostic sign restriction in the sense of Uhlig(2005)
  nosign <- which(is.na(Signmat))
  if (length(nosign) == 0) {
    agnostic <- 0
  } else if (length(nosign) == 1) {
    agnostic <- 1
  } else {
    agnostic <- 2
    if (Hist){
      warning("Histogram is only allowed for case where there exists only one unrestricted sign.")
      Hist <- FALSE
    }
  }


  while (i <= Iter & im <= Itermax) {

    theta.rand <- runif(n = M, min = 0, max = pi)
    B_temp     <- givens.Q.svars(theta = theta.rand, B_chol)

    # make sure impacts on the diagonal are all positive
    if (PositiveDiag) {
      for (j in 1:K) {
        if (B_temp[j,j] < 0) {
          B_temp[,j] <- B_temp[,j] * -1
        }
      }
    }

    if (RHorizont == 0) {

      # only the impact matrix is subject to sign check
      sign_test  <- sign(B_temp) == Signmat

      if (all(sign_test, na.rm = TRUE)){
        #B.win[,,i] <- B_temp
        B.win[[i]] <- B_temp
        i          <- i + 1
      }
    } else {

      # multi-step IRFs are subject to sign check
      irf.temp  <- get.my.irf(Model = Model, Step = RHorizont, Normalize = Normalize, Integrate = Integrate, Bmat = B_temp)
      if (check.irf.sign(irf.temp, Signmat, RHorizont)) {
        B.win[,,i] <- B_temp
        i          <- i + 1
      }
    }

    im <- im + 1

  }

  if (length(B.win) == 0) {
    # no impact matrix that fulfills the sign restrictions is found
    erg   <- NULL
    warning("No successful draw is found!")
    return(erg)
  } else {
    if (length(B.win) == 1) {
      # only one matrix is found
      B.win <- B.win[[1]]

      # IRFs
      IRFs  <- get.my.irf(Model = Model, Step = Step, Normalize = Normalize, Integrate = Integrate, Bmat = B.win)
      IRF.M <- IRF.L <- IRF.U <- IRFs

      # Fry and Pagan’s (2011) MT method
      if(MT) {
        IRF.MT <- IRFs
      }

      warning("Only one successful draw is found!")

    } else {
      # more matrices are found
      N     <- length(B.win)

      if (N < Iter) {
        warning(paste0("Only ", N, " successful draws are found!"))
      }

      # IRFs
      IRFs  <- array(0, c(K, K, Step + 1, N))
      for (n in 1:N) {
        IRFs[,,,n] <- get.my.irf(Model = Model, Step = Step, Normalize = Normalize, Integrate = Integrate, Bmat = B.win[[n]])
      }

      # Median and sample quantiles corresponding to the given probabilities
      IRF.M <- IRF.L <- IRF.U <- array(0, c(K, K, Step + 1))
      for (h in 1:(Step + 1)) {
        for (i in 1:K) {
          for (j in 1:K) {
            IRF.M[i,j,h] <- median(IRFs[i,j,h,])
            IRF.L[i,j,h] <- quantile(IRFs[i,j,h,], probs = c((1-Ci)/2, (1+Ci)/2))[1]
            IRF.U[i,j,h] <- quantile(IRFs[i,j,h,], probs = c((1-Ci)/2, (1+Ci)/2))[2]
          }
        }
      }

      # Fry and Pagan’s (2011) MT method
      if(MT) {
        IRF.MT <- MT.irf(IRFs = IRFs, IRF.M = IRF.M)
      }
    }
  }

    # plot IRFs
    if (Plot) {
      if(!MT) {
        IRF.MT <- IRF.M
      }
      Response.M           <- matrix(0, ncol = K^2 + 1, nrow = Step + 1)
      colnames(Response.M) <- rep("h", K^2 + 1)
      Response.M[,1]       <- 0:Step
      Response.MT <- Response.L <- Response.U <- Response.M
      c <- 1

      for(i in 1:K){
        for(j in 1:K){
          c <- c + 1
          Response.M[,c]  <- IRF.M[i,j,]
          Response.L[,c]  <- IRF.L[i,j,]
          Response.U[,c]  <- IRF.U[i,j,]
          Response.MT[,c] <- IRF.MT[i,j,]
          colnames(Response.M)[c] <- paste("epsilon[", Epsname[j],"]", "%->%", Varname[i])
        }
      }

      colnames(Response.MT) <- colnames(Response.L) <- colnames(Response.U) <- colnames(Response.M)

      Response.M    <- melt(as.data.frame(Response.M), id = 'h')
      Response.L    <- melt(as.data.frame(Response.L), id = 'h', value.name = "L")
      Response.U    <- melt(as.data.frame(Response.U), id = 'h', value.name = "U")
      Response.MT   <- melt(as.data.frame(Response.MT), id = 'h')
      Response.Ci   <- left_join(Response.L, Response.U, by = c("variable", "h"))
      Response.plot <- left_join(Response.M, Response.Ci, by = c("variable", "h")) %>% add_column(Label = "median")

      if (!MT) {
        IRF.plot <- ggplot(data = Response.plot, aes(x = h, y = value)) + geom_line() + geom_hline(yintercept = 0, color = 'red') +
          facet_wrap(~variable, scales = "free_y", labeller = label_parsed) +
          geom_line(aes(y = L), linetype = "dashed", alpha = 0.9) +
          geom_line(aes(y = U), linetype = "dashed", alpha = 0.9) +
          geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.1) +
          xlab("Horizont") + ylab("Response") +
          theme_bw()
      } else {
        Response.plot.MT <- left_join(Response.MT, Response.Ci, by = c("variable", "h")) %>% add_column(Label = "median target (MT)")
        Response.plot <- rbind(Response.plot, Response.plot.MT)
        IRF.plot <- ggplot(data = Response.plot, aes(x = h, y = value, group = Label)) +
          geom_line(aes(linetype = Label, color = Label), alpha = 0.9, size = 0.9) +
          geom_hline(yintercept = 0, color = 'red') +
          facet_wrap(~variable, scales = "free_y", labeller = label_parsed) +
          geom_line(aes(y = L), linetype = "dashed", alpha = 0.9) +
          geom_line(aes(y = U), linetype = "dashed", alpha = 0.9) +
          geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.1) +
          xlab("Horizont") + ylab("Response") +
          theme_bw()
      }

      plot(IRF.plot)
    }

    erg <- list()
    erg[["Bs"]] <-  array(unlist(B.win), c(K, K, N))
    erg[["Varname"]] <- Varname
    erg[["Tob"]] <- Tob
    erg[["Beta"]] <- Beta
    erg[["Z"]] <- Z
    erg[["p"]] <- p
    erg[["dat"]] <-  Model$y
    erg[["epsilon.M"]] <- solve(IRF.M[,,1]) %*% t(u)
    erg[["method"]] <- "sr"
    erg[["IRF.M"]] <- IRF.M
    if (MT) {
      erg[["IRF.MT"]] <- IRF.MT
      erg[["epsilon.MT"]] <- solve(IRF.MT[,,1]) %*% t(u)
      erg[["method"]] <- c("sr", "mt")
    }
    if (agnostic == 1) {
      total <- unlist(lapply(B.win, '[[', nosign))
      unclear <- sum(sign(quantile(total, c((1-Ci)/2, (1+Ci)/2)))) == 0
      erg[["unrest"]] <- list("total" = total, "Positive" = mean(total>0), "Negative" = mean(total<0), "Unclear" = unclear)
      if (Hist) {
        plot(ggplot(data = as.data.frame(total), aes(x = total)) + labs(x= "") +
          geom_histogram(bins = 50, colour = "black", fill = "white") +
          geom_vline(aes(xintercept = median(total)), color = "deepskyblue", linetype = "dashed", size = 1) +
          geom_vline(aes(xintercept = 0), color = "red", linetype = "solid", size = 0.8) +
          theme_bw())
      }
    }
    erg[["IRF.L"]] <- IRF.L
    erg[["IRF.U"]] <- IRF.U
    erg[["# of draws"]] <- im - 1
    if (Plot){
      erg[["IRF.plot"]] <- Response.plot
    }
    class(erg) <- "my.id"

  return(erg)
}

# Function to check sign of IRF of multiple horizonts
check.irf.sign <- function(IRF, Signmat, RHorizont){

  h <- 1
  while (h <= (RHorizont + 1)) {

    sign_test <- sign(IRF[,,h]) == Signmat

    if (all(sign_test, na.rm = TRUE)){
      erg <- TRUE
      h   <- h + 1
    } else {
      erg <- FALSE
      break
    }
  }
  return(erg)
}

# Fry and Pagan’s (2011) MT method
MT.irf <- function(IRFs, IRF.M){
  N    <- dim(IRFs)[4]
  K    <- dim(IRFs)[1]
  Step <- dim(IRFs)[3] - 1

  IRF.gap <- array(NA, dim = dim(IRFs))
  for (n in 1:N) {
    IRF.gap[,,,n] <- IRFs[,,,n] - IRF.M
  }

  IRF.sd  <- array(0, c(K, K, Step + 1))
  for (h in 2:(Step + 1)) {
    for (i in 1:K) {
      for (j in 1:K) {
        IRF.sd[i,j,h]   <- sd(IRFs[i,j,h,])
        IRF.gap[i,j,h,] <- IRF.gap[i,j,h,]/IRF.sd[i,j,h]
      }
    }
  }

  gap <- rep(NA, N)
  for (n in 1:N) {
    gap[n] <- crossprod(IRF.gap[,,,n])
  }

  erg <- IRFs[,,,which.min(gap)]

}
