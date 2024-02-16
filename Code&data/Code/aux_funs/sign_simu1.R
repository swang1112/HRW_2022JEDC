if (server_mode){
  Rcpp::sourceCpp('../../../aux_funs/aux_funs.cpp')
  # Uhlig's agnostic restrictions on the impact multiplier
  sign_id0 <- function(Model, MC,  Epsname = NULL, simu_modus = FALSE, cl = NULL){
    
    # test --------------------------------------------------------------------
    # Model = var4
    # MC = 2000
    # Epsname =  c("d", "s", "mp")
    
    # end test ----------------------------------------------------------------
    
    # internal functions ------------------------------------------------------
    .norm_irf_3 = function(dat){
      dat       = dat %>% as.data.frame() 
      # scale_fac = dat$V9[1]/0.25
      # dat$V9 =  dat$V9/scale_fac
      # dat$V7 =  dat$V7/scale_fac
      # dat$V8 =  dat$V8/scale_fac
      dat[,c('V7','V8','V9')]
    }
    
    
    .QR1 = function(x){
      foo = qr(x)
      Q = qr.Q(foo)
      R = qr.R(foo)
      for (kk in 1:ncol(x)) {
        if (R[kk,kk]<0){
          Q[,kk] = Q[,kk] * -1
        }
      }
      Q
    }
    
    .sign_check = function(x, sr){
      sign_test  <- sign(x) == sr
      all(sign_test, na.rm = TRUE)
    }
    
    .postiv_sign = function(x){
      for (kk in 1:ncol(x)) {
        if (x[kk,kk]<0) {
          x[,kk] = x[,kk] * -1
        }
      }
      x
    }
    
    .sr_b_mprf = function(x){
      
      C = t(chol(x$Sigma))
      pass = 0
      while (pass == 0) {
        G = matrix(rnorm(K*K), nrow = K, ncol = K)
        Q = .QR1(G)
        B = C %*% Q
        B = .postiv_sign(B)
        pass    = .sign_check(B, matrix(c(1,1,1,-1,1,1,NA,-1,1), nrow = 3, ncol = 3))
        
      }
      
      A_hat = matrix(x$b_vec, nrow = K, byrow = T)[,-1]
      list("IRFs" = IRF_fast(A_hat = A_hat, B_hat = B, horizon = 16), "Bmats" = B)
      
    }
    
    
    
    p       <- Model$p
    y       <- Model$y[-c(1:p),]
    K       <- Model$K
    u       <- resid(Model)
    Tob     <- Model$obs
    Varname <- names(Model$varresult)
    Beta    <- Model %>% vars::Bcoef() %>% t()
    if(is.null(Epsname)){Epsname <- Varname}
    if(is.null(cl)){cl <- parallel::detectCores()-1}
    Covmat  <- crossprod(u) / (Tob - 1)
    
    Z = Ylag(Model$y, p, 1)
    y_vec = matrix(y, ncol = 1)
    b_vec = var_est1(y_vec, t(Z), K)
    
    
    # Beta_vec = matrix(Beta, ncol = 1)
    # kronecker(diag(K), Z) %*% Beta_vec
    # matrix(kronecker(diag(K), Z) %*% Beta_vec, ncol = 3)
    # crossprod(y  -  matrix(kronecker(diag(K), Z) %*% Beta_vec, ncol = 3)) / (Tob-1)
    
    # matrix(b_vec, nrow = K, byrow = T)[,-1]
    
    MCMC = Gibbs_VAR_noninfo(y_vec, b_vec, Z, Covmat, MC = MC, p, K, Tob)
    
    #   
    # foo1 = MCMC %>% lapply("[[", 1)
    # foo2 = MCMC %>% lapply("[[", 2)
    # 
    # par(mfrow = c(3,3), mar = c(0,0,0,0))
    # for (i in 1:9) {
    #   idx =  sample(1:length(b_vec), 1)
    #   foo1 %>% lapply("[[", idx) %>%unlist %>% plot(type = 'l')
    #   abline(h = b_vec[idx], col = 2)
    # }
    # par(mfrow = c(1,1))
    # 
    # 
    # par(mfrow = c(3,3), mar = c(0,0,0,0))
    # for (i in 1:9) {
    #   foo2 %>% lapply("[[", i) %>% unlist %>% plot(type = 'l')
    #   abline(h = Covmat[i], col = 2)
    # }
    # par(mfrow = c(1,1))
    # 
    
    SR_model = pbapply::pblapply(MCMC, .sr_b_mprf, cl = cl)
    
    
    if (simu_modus){
      # IRFs (partial identification of one shock)
      IRFs = matrix(NA, nrow = length(SR_model), ncol = K^2*16)
      for (i in 1:length(SR_model)) {
        dat = SR_model[[i]]$IRFs %>% sapply(function(x) x) %>% t
        IRFs[i,] = dat %>% as.data.frame %>% unlist %>%  matrix(nrow = 1)
      }
      
      # point estimates of Inoue and Kilian (2022)
      loss_val = rep(NA, length(SR_model))
      loss_val = IK_point(loss_val, IRFs, 2)
      winner   = which.min(loss_val)
      
      
      
      
      IRF_list = IRF_fast(A_hat = matrix(MCMC[[winner]]$b_vec, nrow = K, byrow = T)[,-1], B_hat = SR_model[[winner]]$Bmats %*% diag(1/diag(SR_model[[winner]]$Bmats)), horizon = 16) 
      IRF.M =  array(0, c(K, K, 16))
      for (h in 1:16) {
        IRF.M[,,h] = IRF_list[[h]]%>% matrix(nrow = K, ncol = K)
      }
      
      IRF.U = 0
      IRF.L = 0
    } else {
      # IRFs (partial identification of one shock)
      IRFs = matrix(NA, nrow = length(SR_model), ncol = K*16)
      for (i in 1:length(SR_model)) {
        dat = SR_model[[i]]$IRFs %>% sapply(function(x) x) %>% t
        IRFs[i,] = dat %>% .norm_irf_3 %>% unlist %>%  matrix(nrow = 1)
      }
      
      # point estimates of Inoue and Kilian (2022)
      loss_val = rep(NA, length(SR_model))
      loss_val = IK_point(loss_val, IRFs, 2)
      winner   = which.min(loss_val)
      
      # Median and sample quantiles corresponding to the given probabilities
      IRF.M = IRFs[winner,]
      IRF.U = IRFs %>% apply(MARGIN = 2, quantile, probs = 0.84)  
      IRF.L = IRFs %>% apply(MARGIN = 2, quantile, probs = 0.16)  
    }
    
    
    
    
    erg <- list()
    erg[["B"]] <- SR_model[[winner]]$Bmats
    erg[["Varname"]] <- Varname
    erg[["Epsname"]] <- Epsname
    erg[["dat"]] <-  Model$y
    erg[["Tob"]] <- Tob
    erg[["Z"]] <- Z
    erg[["p"]] <- p
    erg[["Beta"]] <-  matrix(MCMC[[winner]]$b_vec, nrow = K, byrow = T)[,-1]
    erg[["dat"]] <-  Model$y
    erg[["epsilon.M"]] <- solve(SR_model[[winner]]$Bmats) %*% t(u)
    erg[["method"]] <- "sr.b"
    erg[["IRF.M"]] <- IRF.M
    erg[["IRF.L"]] <- IRF.L
    erg[["IRF.U"]] <- IRF.U
    
    
    
    total <- unlist(lapply(SR_model, '[[', 2) %>% lapply('[[', 7))
    unclear <- sum(sign(quantile(total, c(0.16, 0.84)))) == 0
    erg[["unrest"]] <- list("total" = total, "Positive" = mean(total>0), "Negative" = mean(total<0), "Unclear" = unclear)
    
    class(erg) <- "my.id"
    return(erg)
  }
  
  # Arias et al.'s restrictions on MP rule
  sign_id1 <- function(Model, MC,  Epsname = NULL, simu_modus = FALSE, cl = NULL){
    
    # test --------------------------------------------------------------------
    # Model = var4
    # MC = 2000
    # Epsname =  c("d", "s", "mp")
    
    # end test ----------------------------------------------------------------
    
    # internal functions ------------------------------------------------------
    .norm_irf_3 = function(dat){
      dat       = dat %>% as.data.frame() 
      # scale_fac = dat$V9[1]/0.25
      # dat$V9 =  dat$V9/scale_fac
      # dat$V7 =  dat$V7/scale_fac
      # dat$V8 =  dat$V8/scale_fac
      dat[,c('V7','V8','V9')]
    }
    
    
    .QR1 = function(x){
      foo = qr(x)
      Q = qr.Q(foo)
      R = qr.R(foo)
      for (kk in 1:ncol(x)) {
        if (R[kk,kk]<0){
          Q[,kk] = Q[,kk] * -1
        }
      }
      Q
    }
    
    .sign_check = function(x, sr){
      sign_test  <- sign(x) == sr
      all(sign_test, na.rm = TRUE)
    }
    
    .postiv_sign = function(x){
      for (kk in 1:ncol(x)) {
        if (x[kk,kk]<0) {
          x[,kk] = x[,kk] * -1
        }
      }
      x
    }
    
    .sr_b_mprf = function(x){
      
      C = t(chol(x$Sigma))
      pass = 0
      while (pass == 0) {
        G = matrix(rnorm(K*K), nrow = K, ncol = K)
        Q = .QR1(G)
        B = C %*% Q
        B = .postiv_sign(B)
        pass_b    = .sign_check(B, matrix(c(1,1,1,-1,1,1,NA,-1,1), nrow = 3, ncol = 3))
        if (pass_b == 0){
          next
        } else {
          B_inv = solve(B)
          pass  = all(c((B_inv[3,1]/B_inv[3,3]<0),(B_inv[3,2]/B_inv[3,3]<0)))
        }
        
      }
      
      A_hat = matrix(x$b_vec, nrow = K, byrow = T)[,-1]
      list("IRFs" = IRF_fast(A_hat = A_hat, B_hat = B, horizon = 16), "Bmats" = B)
      
    }
    
    
    
    p       <- Model$p
    y       <- Model$y[-c(1:p),]
    K       <- Model$K
    u       <- resid(Model)
    Tob     <- Model$obs
    Varname <- names(Model$varresult)
    Beta    <- Model %>% vars::Bcoef() %>% t()
    if(is.null(Epsname)){Epsname <- Varname}
    if(is.null(cl)){cl <- parallel::detectCores()-1}
    Covmat  <- crossprod(u) / (Tob - 1)
    
    Z = Ylag(Model$y, p, 1)
    y_vec = matrix(y, ncol = 1)
    b_vec = var_est1(y_vec, t(Z), K)
    
    
    # Beta_vec = matrix(Beta, ncol = 1)
    # kronecker(diag(K), Z) %*% Beta_vec
    # matrix(kronecker(diag(K), Z) %*% Beta_vec, ncol = 3)
    # crossprod(y  -  matrix(kronecker(diag(K), Z) %*% Beta_vec, ncol = 3)) / (Tob-1)
    
    # matrix(b_vec, nrow = K, byrow = T)[,-1]
    
    MCMC = Gibbs_VAR_noninfo(y_vec, b_vec, Z, Covmat, MC = MC, p, K, Tob)
    
    
    SR_model = pbapply::pblapply(MCMC, .sr_b_mprf, cl = cl)
    
    
    
    if (simu_modus){
      # IRFs (partial identification of one shock)
      IRFs = matrix(NA, nrow = length(SR_model), ncol = K^2*16)
      for (i in 1:length(SR_model)) {
        dat = SR_model[[i]]$IRFs %>% sapply(function(x) x) %>% t
        IRFs[i,] = dat %>% as.data.frame %>% unlist %>%  matrix(nrow = 1)
      }
      
      # point estimates of Inoue and Kilian (2022)
      loss_val = rep(NA, length(SR_model))
      loss_val = IK_point(loss_val, IRFs, 2)
      winner   = which.min(loss_val)
      
      IRF_list = IRF_fast(A_hat = matrix(MCMC[[winner]]$b_vec, nrow = K, byrow = T)[,-1], B_hat = SR_model[[winner]]$Bmats %*% diag(1/diag(SR_model[[winner]]$Bmats)), horizon = 16) 
      IRF.M =  array(0, c(K, K, 16))
      for (h in 1:16) {
        IRF.M[,,h] = IRF_list[[h]]%>% matrix(nrow = K, ncol = K)
      }
      
      IRF.U = 0
      IRF.L = 0
    } else {
      # IRFs (partial identification of one shock)
      IRFs = matrix(NA, nrow = length(SR_model), ncol = K*16)
      for (i in 1:length(SR_model)) {
        dat = SR_model[[i]]$IRFs %>% sapply(function(x) x) %>% t
        IRFs[i,] = dat %>% .norm_irf_3 %>% unlist %>%  matrix(nrow = 1)
      }
      
      # point estimates of Inoue and Kilian (2022)
      loss_val = rep(NA, length(SR_model))
      loss_val = IK_point(loss_val, IRFs, 2)
      winner   = which.min(loss_val)
      
      # Median and sample quantiles corresponding to the given probabilities
      IRF.M = IRFs[winner,]
      IRF.U = IRFs %>% apply(MARGIN = 2, quantile, probs = 0.84)  
      IRF.L = IRFs %>% apply(MARGIN = 2, quantile, probs = 0.16)  
    }
    
    
    
    
    erg <- list()
    erg[["B"]] <- SR_model[[winner]]$Bmats
    erg[["Varname"]] <- Varname
    erg[["Epsname"]] <- Epsname
    erg[["dat"]] <-  Model$y
    erg[["Tob"]] <- Tob
    erg[["Z"]] <- Z
    erg[["p"]] <- p
    erg[["Beta"]] <-  matrix(MCMC[[winner]]$b_vec, nrow = K, byrow = T)[,-1]
    erg[["dat"]] <-  Model$y
    erg[["epsilon.M"]] <- solve(SR_model[[winner]]$Bmats) %*% t(u)
    erg[["method"]] <- "sr.b"
    erg[["IRF.M"]] <- IRF.M
    erg[["IRF.L"]] <- IRF.L
    erg[["IRF.U"]] <- IRF.U
    
    total <- unlist(lapply(SR_model, '[[', 2) %>% lapply('[[', 7))
    unclear <- sum(sign(quantile(total, c(0.16, 0.84)))) == 0
    erg[["unrest"]] <- list("total" = total, "Positive" = mean(total>0), "Negative" = mean(total<0), "Unclear" = unclear)
    
    class(erg) <- "my.id"
    return(erg)
  }
  
  sign_plot1 = function(SR){
    IRF.M = SR$IRF.M %>% matrix(nrow = 16, ncol = 3) %>% as.data.frame %>% dplyr::mutate(h= 0:15) %>% reshape2::melt(id.vars = "h") %>% dplyr::mutate(Label = "M")
    IRF.L = SR$IRF.L %>% matrix(nrow = 16, ncol = 3) %>% as.data.frame %>% dplyr::mutate(h= 0:15) %>% reshape2::melt(id.vars = "h") %>% dplyr::mutate(Label = "L")
    IRF.U = SR$IRF.U %>% matrix(nrow = 16, ncol = 3) %>% as.data.frame %>% dplyr::mutate(h= 0:15) %>% reshape2::melt(id.vars = "h") %>% dplyr::mutate(Label = "U")
    
    scale_fac = dplyr::filter(IRF.M, h == 0 & variable == "V3")$value/0.25 
    IRF.M$value = IRF.M$value/scale_fac
    
    
    Response = cbind(IRF.M, IRF.L$value, IRF.U$value)
    colnames(Response)[5:6] = c("L", "U")
    
    my_labeller <- as_labeller(c(V1 = paste("epsilon[","mp","]", "%->%", SR$Varname[1]), 
                                 V2 = paste("epsilon[","mp","]", "%->%", SR$Varname[2]), 
                                 V3 = paste("epsilon[","mp","]", "%->%", SR$Varname[3])), 
                               default = label_parsed)
    
    ggplot(data = Response, aes(x = h, y = value)) +
      geom_line() + geom_hline(yintercept = 0, color = 'red') +
      facet_wrap(~variable, scales = "free_y", labeller = my_labeller, ncol = 1) +
      geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
      xlab("Horizon") + ylab("Response") +
      theme_bw()
    
  }
  
} else {
  source('../aux_funs/sign_id1.R')
}


simulate.sign1 <- function(Data, AOA, Subsample = NULL, Iter,  Signmat, TrueSign, 
                          MT = FALSE, Varname = c("x", "pi", "r"), Bmat,
                          Plot = FALSE, Hist = FALSE, Arias = FALSE, cl, cl_sub = NULL){
  

# test --------------------------------------------------------------------
  # Data = Y
  # Iter = 50
  # Subsample = NULL
  # Varname = c("x", "pi", "r")
  # Bmat = DSGE.Par$B
  # MT = F
  # Arias = T
  # Plot = FALSE
  # Hist = FALSE
  # cl_sub = NULL
  # cl = 8
# end test ----------------------------------------------------------------

    
  
  
  Size     <- length(Data)
  K        <- ncol(Data[[1]]$Y)
  p        <- 2
  Tob      <- nrow(Data[[1]]$Y)
  Step     <- dim(AOA$L)[3] - 1
  eva.temp <- array(NA, c(K, K, Step + 1, Size))
  mse.temp <- matrix(rep(NA, Size*3), nrow = Size, ncol = 3)
  
  # check whether to use Fry and Pagan’s (2011) Median Target method
  if (MT) {
    eva.temp.MT   <- eva.temp
  }
  
  # check btw. full sign restriction or argnostic sign restriction in the sense of Uhlig(2005)
  nosign <- which(is.na(Signmat))
  if (length(nosign) == 0) {
    agnostic <- FALSE
  } else if (length(nosign) == 1) {
    agnostic <- TRUE
    correct  <- rep(NA, Size)
    unclear  <- rep(NA, Size)
  } else {stop("For now, only one unrestricted sign is allowed.")}
  
  
  
  .simu_kernel = function(x){
    # Data-generation
    Y.temp <- x
    
    # Subsample
    if (!is.null(Subsample)){
      Y.temp$Y     <- Y.temp$Y[Subsample[1]:Subsample[2],]
      Y.temp$shock <- Y.temp$shock[Subsample[1]:Subsample[2],]
      Y.temp$iv    <- Y.temp$iv[Subsample[1]:Subsample[2],]
    }
    
    # Estimation
    var2 <- vars::VAR(Y.temp$Y, p = 2, type = "none")
    
    # Identification
    if (Arias){
      ID.temp = sign_id1(Model = var2, MC = Iter, Epsname =  c("d", "s", "mp"), simu_modus = TRUE, cl = cl_sub)
    } else {
      ID.temp = sign_id0(Model = var2, MC = Iter, Epsname =  c("d", "s", "mp"), simu_modus = TRUE, cl = cl_sub)
    }
    
    

    IRF.temp            <- ID.temp$IRF.M
    row.names(IRF.temp) <- Varname
    list('eva' = Evaluate.AOA(IRF.temp, AOA = AOA),
         'mse' = eva.RMSE(B = Bmat, Bhat = ID.temp$B, xcol = 3),
         'unclear' = ID.temp$unrest$Unclear,
         'correect' = ID.temp$unrest$Negative)
  }
  
  Out = pbapply::pblapply(Data, .simu_kernel, cl = cl)
  
  Out1 = Out %>% lapply('[[',1)
  Out2 = Out %>% lapply('[[',2)
  unclear = Out %>% lapply('[[',3) %>% unlist
  correct = Out %>% lapply('[[',4) %>% unlist
  rm('Out')
  for (i in 1:Size) {
    eva.temp[,,,i] = Out1[[i]]
    mse.temp[i,] = Out2[[i]]
  }
  
  # 
  # start.time <- Sys.time()
  # cat("\r", "...calculating finish time...")
  # 
  # for (i in 1:Size) {
  #   # Data-generation
  #   Y.temp <- Data[[i]]
  #   
  #   # Subsample
  #   if (!is.null(Subsample)){
  #     Y.temp$Y     <- Y.temp$Y[Subsample[1]:Subsample[2],]
  #     Y.temp$shock <- Y.temp$shock[Subsample[1]:Subsample[2],]
  #     Y.temp$iv    <- Y.temp$iv[Subsample[1]:Subsample[2],]
  #   }
  #   
  #   # Estimation
  #   var2 <- vars::VAR(Y.temp$Y, p = 2, type = "none")
  #   
  #   # Identification
  #   if (Arias){
  #     ID.temp = sign_id1(Model = var2, MC = Iter, Epsname =  c("d", "s", "mp"), simu_modus = TRUE)
  #   } else {
  #     ID.temp = sign_id0(Model = var2, MC = Iter, Epsname =  c("d", "s", "mp"), simu_modus = TRUE)
  #   }
  #   
  #   
  #     IRF.temp            <- ID.temp$IRF.M
  #     row.names(IRF.temp) <- Varname
  #     eva.temp[,,,i]      <- Evaluate.AOA(IRF.temp, AOA = AOA)
  #     mse.temp[i,]         <- eva.RMSE(B = Bmat, Bhat = ID.temp$B, xcol = 3)
  #     
  #     # document frequency of detection of the correct sign
  #     if (agnostic) {
  #       unclear[i] <- ID.temp$unrest$Unclear
  #       if (TrueSign[nosign] == 1){
  #         correct[i] <- ID.temp$unrest$Positive
  #       } else if (TrueSign[nosign] == -1) {
  #         correct[i] <- ID.temp$unrest$Negative
  #       } else {stop("Invalid argument in TrueSign!")}
  #     }
  #     
  #     # evaluate MT method
  # 
  #   
  #   
  #   # print progress
  #   progress(i, max = Size, start.time = start.time)
  # }
  
    
    
    
    erg.M <- array(NA, c(K, K, Step + 1))
    for (h in 1:(Step + 1)) {
      for (i in 1:K) {
        for (j in 1:K) {
          erg.M[i,j,h] <- mean(eva.temp[i,j,h,])
        }
      }
    }
    row.names(erg.M) <- Varname
    
    erg <- list("M" = erg.M, "ump" = NaN, "mse.mp" = colMeans(mse.temp, na.rm = T))
    
    if (Plot == TRUE){
      plot(plot.simu(erg.M))
    }
  
  
  if (agnostic) {
    entry <- which(is.na(Signmat), arr.ind = T)
    cat(paste0("\n", "Frequency of detecting the correct sign pattern at position ",
               "row " , entry[1], " column ", entry[2], ": ",
               round(mean(correct)*100, digits = 2), " percent.", "\n"))
    erg[["unrest"]] <- list("position" = entry, "freq" = correct, "freq_mean" = mean(correct), "Unclear" = mean(unclear))
    if (Hist) {
      plot(ggplot(data = as.data.frame(correct), aes(x = correct)) + labs(x= "") +
             geom_histogram(bins = 20, colour = "black", fill = "#FF6666", alpha = 0.1) +
             geom_vline(aes(xintercept = mean(correct)), color = "deepskyblue", linetype = "dashed", size = 1) +
             theme_bw())
    }
  }
  
  return(erg)
}
