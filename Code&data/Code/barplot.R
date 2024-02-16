
Start = 1979.75
End1 = 1982.5
End2 = 1983.5


Volcker_shock <- function(x, Punkt = 1979.75, End){
  series = 2
  Partial = 3
  Epsname = c("d", "s", "mp")
  Start = Start
  Freq = 4
  
  if(class(x) == "my.id"){
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
    
    
    
    if (x$method[1] == "sr"){
      A_hat <- t(x$Beta)[,  1 : (x$p * k)]
      B_hat <- x$IRF.MT[,,1] #median target
      horizon <- (End - Punkt)*Freq + 1
      
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
      
      
      s.errors <- data.frame("time" = s.time, "eps" = s.errors) %>% filter(time == Punkt) %>% select(eps)
      s.impulse <- impulse[,(k*(series-1)+Partial +1)]
      
      erg <- s.errors$eps*sum(s.impulse)
      
      
    } else if (x$method[1] == "sr.b"){
      A_hat <- x$Beta
      B_hat <- x$B #median target
      horizon <- (End - Punkt)*Freq + 1
      A_hat <- x$Beta
      
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
      
      
      s.errors <- data.frame("time" = s.time, "eps" = s.errors) %>% filter(time == Punkt) %>% select(eps)
      s.impulse <- impulse[,(k*(series-1)+Partial +1)]
      
      erg <- s.errors$eps*sum(s.impulse)
      
      
    } else if (x$method[1] == "iv"){
      
      A_hat <- t(x$Beta)[,  1 : (x$p * k)]
      B_hat <- matrix(x$B, k, k)
      
      horizon <- (End - Punkt)*Freq + 1
      
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
      
      
      s.errors <- x$epsilon
      s.time <- x$eps.ts[,1]
      
      s.errors <- data.frame("time" = s.time, "eps" = s.errors) %>% filter(time == Punkt) %>% select(eps)
      s.impulse <- impulse[,(k*(series-1)+Partial +1)]
      
      erg <- s.errors$eps*sum(s.impulse)
      
    }
  } else if (class(x) == "svars"){
    
    k <- x$K
    p <- x$p
    horizon <- (End - Punkt)*Freq + 1
    
    impulse <- irf(x, n.ahead = horizon)
    
    
    s.errors <- solve(x$B) %*% t(resid( x$VAR))
    s.errors <- s.errors[Partial,]
    s.time <- time(x$y)[-(1:p)]
    
    s.errors <- data.frame("time" = s.time, "eps" = s.errors) %>%  dplyr::filter(time == Punkt) %>% dplyr::select(eps)
    s.impulse <- impulse$irf[,(k*(series-1)+Partial +1)]
    
    erg <- s.errors$eps*sum(s.impulse)
    
    
  }
  
  return(erg)
  
}


# cumulated ---------------------------------------------------------------
df_pro <- data.frame("period" = c("1979Q4 - 1982Q3", "1979Q4 - 1983Q3"),
                     "SR+" = c(quasi.hd(SR, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                               quasi.hd(SR, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                     "dCov" = c(quasi.hd(dCov, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                                quasi.hd(dCov, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                     "CV" = c(quasi.hd(CV, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                              quasi.hd(CV, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                     "SW" = c(quasi.hd(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                              quasi.hd(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                     "SZ" = c(quasi.hd(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                              quasi.hd(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)),
                     "RR" = c(quasi.hd(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                              quasi.hd(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)),
                     "SZ_pro" = c(quasi.hd(IV.SZ_pro, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                                  quasi.hd(IV.SZ_pro, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                     "RR_pro" = c(quasi.hd(IV.RR_pro, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                                  quasi.hd(IV.RR_pro, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4))) %>% melt(id = c("period"))
df_pro %<>% add_column(lab = "old") %>% mutate(lab = ifelse(variable == "SZ_pro" | variable == "RR_pro", "new", "old"))
df_pro[13:14,]$variable <- "SZ"
df_pro[15:16,]$variable <- "RR"

df_pro$variable = factor(c("SR+", "SR+", "dCov", "dCov", "CV", 
                           "CV",   "SW",   "SW",   "SZ",   "SZ",
                           "RR",   "RR",   "SZ",   "SZ",   "RR",   "RR" ),
                         levels = c("SR+", "dCov", "CV", "SW", "SZ", "RR" ))




#df_pro %<>% mutate(value = value/4)
df_pro %<>% add_column(type = "agg")



# shock of 1979Q4 ---------------------------------------------------------

df_1979Q4 <- data.frame("period" = c("1979Q4 - 1982Q3", "1979Q4 - 1983Q3"),
                        "SR+" = c(Volcker_shock(SR, End = End1),
                                  Volcker_shock(SR, End = End2)),
                        "dCov" = c(Volcker_shock(dCov, End = End1),
                                   Volcker_shock(dCov, End = End2)),
                        "CV" = c(Volcker_shock(CV, End = End1),
                                 Volcker_shock(CV, End = End2)),
                        "SW" = c(Volcker_shock(IV.SW, End = End1),
                                 Volcker_shock(IV.SW, End = End2)),
                        "SZ" = c(Volcker_shock(IV.SZ, End = End1),
                                 Volcker_shock(IV.SZ, End = End2)),
                        "RR" = c(Volcker_shock(IV.RR, End = End1),
                                 Volcker_shock(IV.RR, End = End2)),
                        "SZ_pro" = c(Volcker_shock(IV.SZ_pro, End = End1),
                                     Volcker_shock(IV.SZ_pro, End = End2)),
                        "RR_pro" = c(Volcker_shock(IV.RR_pro, End = End1),
                                     Volcker_shock(IV.RR_pro, End = End2)))  %>% melt(id = c("period")) #%>% mutate(value = value * 0.25)
df_1979Q4 %<>% add_column(lab = "old") %>% mutate(lab = ifelse(variable == "SZ_pro" | variable == "RR_pro", "new", "old"))
df_1979Q4[13:14,]$variable <- "SZ"
df_1979Q4[15:16,]$variable <- "RR"

df_1979Q4$variable = factor(c("SR+", "SR+", "dCov", "dCov", "CV", 
                              "CV",   "SW",   "SW",   "SZ",   "SZ",
                              "RR",   "RR",   "SZ",   "SZ",   "RR",   "RR" ),
                            levels = c("SR+", "dCov", "CV", "SW", "SZ", "RR" ))


df_1979Q4 %<>% add_column(type = "sig")

df_all <- rbind(df_pro, df_1979Q4)

df_all %>% ggplot(aes(x = variable, y = value)) + 
  geom_hline(yintercept = 0, color = "red", alpha = 0.8) +
  facet_wrap(~period) +
  geom_bar(data=subset(df_all, lab == "old" & type == "agg"), aes(fill = variable),
           stat="identity", alpha = .5, width=0.5,  color="#999999", linetype = "solid", size=0.8) + 
  geom_bar(data=subset(df_all, lab == "new" & type == "agg"), aes(fill = variable),
           stat="identity", alpha = .5, width=0.5,  color="#999999", linetype = "dashed", size=0.8) + 
  geom_bar(data=subset(df_all, lab == "old" & type == "sig"), aes(fill = variable),
           stat="identity", alpha = 1, width=0.5,  color="#999999", linetype = "solid", size=0.8) + 
  geom_bar(data=subset(df_all, lab == "new" & type == "sig"), aes(fill = variable),
           stat="identity", alpha = 1, width=0.5,  color="#999999", linetype = "dashed", size=0.8) + 
  xlab("") + ylab(" ") + theme_bw() + theme(legend.position = "none") + scale_fill_ochre(palette="olsen_qual")


