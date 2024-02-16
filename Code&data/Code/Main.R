# Load functions ----------------------------------------------------------
rm(list = ls())
Functions <- devtools::as.package("Functions")
devtools::load_all("Functions")
library(ggfortify)
library(gridExtra)
library(magrittr)
library(svars)
library(readr)

#devtools::install_github("ropenscilabs/ochRe")
library(ochRe)

# Import data -------------------------------------------------------------
USA_Tri_raw   <- read.csv("../Data/USA_Tri.csv")
USA_Tri       <- ts(USA_Tri_raw[,-1], frequency = 4, start = c(floor(USA_Tri_raw$Year[1]), (USA_Tri_raw$Year[1] %% 1) / 0.25 + 1 ))
#autoplot(USA_Tri, ncol = 1) + theme_bw()

# VAR Spesicifation -------------------------------------------------------
## lag length
#VARselect(USA_Tri, lag.max = 8)
var4 <- VAR(USA_Tri, p = 4, type = "const")
#var4 %>% summary
#var4 %>% resid %>% t() %>% write.table(file = "NonFundamentalness/resid.txt", col.names = F, row.names = F)
#apply(USA_Tri, 1, function(x){x-colMeans(USA_Tri)}) %>% t() %>% write.table(file = "NonFundamentalness/dat.txt", col.names = F, row.names = F)
## LM test for serial correlation
#serial.test(var4, type = "BG")
#serial.test(var4, type = "ES")

## detecting structural break
#set.seed(1234)
#chow.test(var4, SB = c(1984, 1), frequency = 4, nboot = 2000) %>% summary 
#set.seed(1234)
#chow.test(var4, SB = c(1979, 3), frequency = 4, nboot = 2000) %>% summary 

## Gaussianity test
#ICtest::FOBIboot(var4 %>% resid, n.boot = 2000, 1) # test against the last 2 components are gaussian: reject
#ICtest::FOBIboot(var4 %>% resid, n.boot = 2000, 2) # test against the last 1 component is gaussian: reject

#tseries::jarque.bera.test(resid(var4)[,1])
#tseries::jarque.bera.test(resid(var4)[,2])
#tseries::jarque.bera.test(resid(var4)[,3])


# Identification ----------------------------------------------------------
## Sign restrictions (impact effect + taylor rule MP)
source('./aux_funs/sign_id1.R')
set.seed(1234)
SR = sign_id1(var4, MC = 2000, Epsname =  c("d", "s", "mp"))


## Identification by change in volatility
set.seed(1234)
CV          <- id.cv(var4, SB = c(1984, 1), frequency = 4)
CV          %>% summary

## Identification by minimizing distance covariance
set.seed(1234)
dCov        <- id.dc(var4)
dCov$B      <- dCov$B[,c(3,2,1)]
dCov$B[,3]  <- dCov$B[,3]*-1
dCov        %>% summary

## Proxy SVAR
# Smets-Wouters (2007)
SW.raw      <- read.csv("../Data/Instruments/SW.csv")
SW          <- ts(SW.raw[,-1], frequency = 4, start = c(floor(SW.raw$Year[1]), (SW.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.SW       <- get.id.iv(var4, instruments = SW, Synchro = T)
IV.SW$B
IV.SW$F_test
# write_csv(IV.SW$eps.ts, "../data/shock/SW_shock.csv")

# Guerkaynak, Sack and Swanson (2005)
GSS.raw     <- read.csv("../Data/Instruments/GSS.csv")
GSS         <- ts(GSS.raw[,-1], frequency = 4, start = c(floor(GSS.raw$Year[1]), (GSS.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.GSS      <- get.id.iv(var4, instruments = GSS, Synchro = T)
IV.GSS$B
IV.GSS$F_test

# Romer Romer (2004)
RR.raw      <- read.csv("../Data/Instruments/RR.csv")
RR          <- ts(RR.raw[,-1], frequency = 4, start = c(floor(RR.raw$Year[1]), (RR.raw$Year[1] %% 1) / 0.25 + 1 ))
RR          <- Remove_AC(RR, lag = 1) 
IV.RR       <- get.id.iv(var4, instruments = RR, Synchro = T)
IV.RR$B
IV.RR$F_test


# Sims Zha (2006)
SZ.raw      <- read.csv("../Data/Instruments/SZ.csv")
SZ          <- ts(SZ.raw[,-1], frequency = 4, start = c(floor(SZ.raw$Year[1]), (SZ.raw$Year[1] %% 1) / 0.25 + 1 ))
SZ          <- Remove_AC(SZ, 2)
IV.SZ       <- get.id.iv(var4, instruments = SZ, Synchro = T)
IV.SZ$B
IV.SZ$F_test


iv2ts <- function(x, varname){
  erg <- data.frame(x %>% time %>% as.numeric(),  x)
  colnames(erg) <- c("time", varname)
  return(erg)}
SW.ts <- iv2ts(SW, "SW")
GSS.ts <- iv2ts(GSS, "GSS")
RR.ts <- iv2ts(RR, "RR")
SZ.ts <- iv2ts(SZ, "SZ")

# Inference ---------------------------------------------------------------
set.seed(1234)
CV.boot     <- svars::wild.boot(CV, n.ahead = 20, nboot = 2000)
set.seed(1234)
dCov.boot   <- svars::wild.boot(dCov, n.ahead = 20, nboot = 2000)


set.seed(1234)
IV.SW.boot  <- get.MBB.fixed(IV.SW, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68)
set.seed(1234)
IV.RR.boot <- get.MBB.fixed(IV.RR, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68)
set.seed(1234)
IV.SZ.boot  <- get.MBB.fixed(IV.SZ, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68)


# Tab 
IV_id.tab <- function(ivid, ivboots){
  temp <- cbind(ivboots$Bmat$sd[,3], ivid$B)
  return(rbind(temp, rep(ivid$F_test,2))) 
}
tab_iv_id <- cbind(IV_id.tab(IV.SW, IV.SW.boot),
                   IV_id.tab(IV.RR, IV.RR.boot),
                   IV_id.tab(IV.SZ, IV.SZ.boot))

rownames(tab_iv_id) <- NULL
kableExtra::kable(tab_iv_id, format = "latex", booktabs = T, digits = 3, linesep = "")



# UMPS and Taylor rule ----------------------------------------------------
UMPS_freq = function(x, iv = FALSE){
  UMPS = Impact = Tayl = rep(0, 2000)
  
  if (iv){
    for (i in 1:2000) {
      Btemp = x$boot_B[,,i]
      if(sign(Btemp[2,3])*sign(Btemp[3,3]) == -1){ Impact[i] = 1}
      q3 = Btemp[,3, drop = F]
      q3 = q3 / sqrt(as.numeric(t(q3) %*% q3))
      B_inv = t(q3) %*% solve(x$boot_Sigma[,,i]) 
      Tayl[i] = all(c((B_inv[1]/B_inv[3]<0),(B_inv[2]/B_inv[3]<0)))
      if (Impact[i] == 1 && Tayl[i] == 1){UMPS[i] = 1}
      rm(list = c('Btemp','q3',  'B_inv'))
    }
  } else {
    for (i in 1:2000) {
      Btemp = x$boot_B[,,i]
      foo   = rep(0, 3)
      for (k in 1:3) {
        if(sign(Btemp[2,k])*sign(Btemp[3,k]) == -1){   foo[k] = 1}
      }
      if (sum(foo) == 1){Impact[i] = 1}
      B_inv = solve(Btemp)
      Tayl[i] = all(c((B_inv[3,1]/B_inv[3,3]<0),(B_inv[3,2]/B_inv[3,3]<0)))
      if (Impact[i] == 1 && Tayl[i] == 1){UMPS[i] = 1}
      rm(list = c('Btemp', 'foo', 'B_inv'))
    }
  }
  
  tibble('UMPS' = mean(UMPS), 'Impact' = mean(Impact),  'Tayl' = mean(Tayl))
}

UMPS_freq(dCov.boot)
UMPS_freq(CV.boot)
UMPS_freq(IV.SW.boot, T)
UMPS_freq(IV.SZ.boot, T)
UMPS_freq(IV.RR.boot, T)



# FAVAR -------------------------------------------------------------------
r <- 3
SDFM <- FALSE
source("Fa.R")
colnames(USA_FAVAR) <- c(colnames(USA_Tri), 1:r)
favar4 <- VAR(USA_FAVAR, p = 4, type = "const")

FaIV.SW      <- get.id.iv(favar4, instruments = SW, Synchro = T)
FaIV.RR       <- get.id.iv(favar4, instruments = RR, Synchro = T)
FaIV.SZ       <- get.id.iv(favar4, instruments = SZ, Synchro = T)
# write_csv(FaIV.SW$eps.ts, "../data/shock/SW_shock_FA.csv")

mp_iv_fa     <- left_join(FaIV.SW$eps.ts, FaIV.SZ$eps.ts, by = "time")
mp_iv_fa     <- left_join(mp_iv_fa, FaIV.RR$eps.ts, by = "time") %>% tidyr::drop_na()
colnames(mp_iv_fa) <- c("time", "SW", "SZ", "RR")

mp_iv_fa_tab <-  mp_iv_fa %>% filter(time == 1974.25 | time == 1979.75 | time == 1988.75 | time == 1994 | time == 1990.75 | time == 1998.75 | time == 2001.25| time == 2002.75)
mp_iv_fa_tab[c(4,5),] <- mp_iv_fa_tab[c(5,4),]
mp_iv_fa_tab %>% t() %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")

# Inference
set.seed(1234)
MP_fa.SW  <- get.MBB.fixed(FaIV.SW, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68) %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp", "X1", "X2", "X3"), Normalize = FaIV.SW$B[3]/0.25) + xlab("") + ylab("")

MP_fa.SW3 <- MP_fa.SW$data %>% filter(!stringr::str_detect(variable, "X")) %>% ggplot( aes(x = h, y = value)) + geom_line() + geom_hline(yintercept = 0, color = 'red') +
  facet_wrap(~variable, scales = "free_y", labeller = label_parsed, ncol = 1) +
  geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab(" ") + ylab(" ") +
  theme_bw() 


# Plot --------------------------------------------------------------------
# standardize IRFs plots
st_fct <- function(sboot, h){
  sboot$true$irf <- sboot$true$irf[1:(h+1),]
  sboot$true$irf[,1] <- 0:h
  for (i in 1:length(sboot[["bootstrap"]])) {
    sboot[["bootstrap"]][[i]]$irf <- sboot[["bootstrap"]][[i]]$irf[1:(h+1),]
    sboot[["bootstrap"]][[i]]$irf[,1] <- 0:h
  }
  return(sboot)}
CV.boot    %<>% st_fct(h = 15)
dCov.boot  %<>% st_fct(h = 15)

#change shock names
colnames(CV.boot$true$irf)[4] <- "epsilon[ mp ] %->% x"
colnames(CV.boot$true$irf)[7] <- "epsilon[ mp ] %->% pi"
colnames(CV.boot$true$irf)[10] <- "epsilon[ mp ] %->% GBR1"

colnames(dCov.boot$true$irf)[4] <- "epsilon[ mp ] %->% x"
colnames(dCov.boot$true$irf)[7] <- "epsilon[ mp ] %->% pi"
colnames(dCov.boot$true$irf)[10] <- "epsilon[ mp ] %->% GBR1"




SR.MP     <- sign_plot1(SR)
CV.MP     <- CV.boot %>% plot.sboot.mod(lowerq = 0.16, upperq = 0.84, Partial = 3, Normalize = CV$B[3,3]/0.25) + xlab("") + ylab("") 
dCov.MP   <- dCov.boot %>% plot.sboot.mod(lowerq = 0.16, upperq = 0.84, Partial = 3, Normalize = dCov$B[3,3]/0.25) + xlab("") + ylab("")
IRF_nonIV <- grid.arrange(SR.MP, CV.MP, dCov.MP, ncol = 3)

SW.MP     <- IV.SW.boot %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp"), Normalize = IV.SW$B[3]/0.25) + xlab("") + ylab("")
RR.MP     <- IV.RR.boot %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp"), Normalize = IV.RR$B[3]/0.25) + xlab("") + ylab("")
SZ.MP     <- IV.SZ.boot %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp"), Normalize = IV.SZ$B[3]/0.25) + xlab("") + ylab("")
IRF_IV <- grid.arrange(SW.MP, RR.MP, SZ.MP, ncol = 3)

# Correlation of MP shocks ------------------------------------------------
mp_sr     <- SR$epsilon.M %>% t() %>% as.data.frame() %>% select(3)
mp_cv     <- solve(CV$B)          %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(3)
mp_dcov   <- solve(dCov$B)        %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(3)
mp_noniv  <- cbind(mp_sr, mp_cv, mp_dcov); colnames(mp_noniv) <- c("SR", "CV", "dCov")
mp_iv     <- left_join(IV.SW$eps.ts, IV.SZ$eps.ts, by = "time")
mp_iv     <- left_join(mp_iv, IV.RR$eps.ts, by = "time") %>% tidyr::drop_na()
colnames(mp_iv) <- c("time", "SW", "SZ", "RR")
mp_noniv  <- cbind(time(var4$y)[-c(1:4)], mp_noniv); colnames(mp_noniv)[1] <- "time"
mp_all    <- left_join(mp_noniv, mp_iv, by = "time") %>% tidyr::drop_na() 
#mp_all <- mp_all %>% select(-1)

cor(mp_all) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 2, linesep = "")

mp_sds <- mp_all %>% select(-1) %>% apply( 2, sd)

mp_all_stand <- mp_all

for (j in 1:dim(mp_all)[1]) {
  for (i in 2:7) {
    mp_all_stand[j,i] <- mp_all[j,i] / mp_sds[i-1]
  } 
}

mp_all_stand %>% select(-1) %>% apply( 2, sd)


mp_all_plot3 <- reshape2::melt( mp_all_stand, id.var = "time")

library(ggplot2)
mp_all_plot3$variable <- factor(  mp_all_plot3$variable, levels = c(   "CV",   "dCov", "SW",   "SZ",   "RR","SR" ))
MP_shock_plot3 <- ggplot(mp_all_plot3, aes(x = time, y = value, group = variable)) +
  geom_line(aes(linetype = variable, color = variable), size = 0.7) + 
  scale_x_continuous(breaks = seq(mp_all_stand$time[1] %>% ceiling(),mp_all_stand$time[nrow(mp_all_stand)], 2)) + 
  xlab("") + ylab(" ")  + theme_bw() + theme(legend.title = element_blank()) + 
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + scale_color_ochre(palette="lorikeet") + ylim(-4.5, 4.5)
MP_shock_plot3

MP_shocks_plotSW_Fa_df <- IV.SW$eps.ts %>% left_join(FaIV.SW$eps.ts, by = "time") %>% filter(time >= mp_all_stand$time[1] & time <= mp_all_stand$time[nrow(mp_all_stand)])
MP_shocks_plotSW_Fa_df[,-1] <- scale(MP_shocks_plotSW_Fa_df[,-1], center = F)
colnames(MP_shocks_plotSW_Fa_df)[-1] <- c("VAR", "FAVAR")
MP_shocks_plotSW_Fa_df_long <-  reshape2::melt(MP_shocks_plotSW_Fa_df, id.var = "time")

MP_shocks_plotSW_Fa <- ggplot(MP_shocks_plotSW_Fa_df_long, aes(x = time, y = value, group = variable)) +
  geom_line(aes(linetype = variable)) + scale_linetype_manual(values = c("solid","longdash")) +
  scale_x_continuous(breaks = seq(mp_all_stand$time[1] %>% ceiling(),mp_all_stand$time[nrow(mp_all_stand)], 2)) +
  xlab("") + ylab(" ")  + theme_bw()  + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
  theme(legend.title = element_blank()) + ylim(-4.5, 4.5)

new_shock_plots <- grid.arrange(MP_shock_plot3, MP_shocks_plotSW_Fa, ncol = 1)


# Event analysis ----------------------------------------------------------
# April 1974, October 1979, December 1988 and February 1994 - were contractionary
# December 1990, October 1998, April 2001, and November 2002  -were expansionary 
mp_all_tab <- mp_all %>% filter(time == 1974.25 | time == 1979.75 | time == 1988.75 | time == 1994 | time == 1990.75 | time == 1998.75 | time == 2001.25| time == 2002.75)
#mp_all_stand_tab <- mp_all_stand_tab %>% select(-SR)
mp_all_tab <- mp_all_tab %>% add_column(SUM = rowSums(mp_all_tab[,-1]))
  
mp_all_tab[c(4,5),] <- mp_all_tab[c(5,4),]
mp_all_tab %>% t() %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")



# purging instruments against CV non-MP shocks ----------------------------
demand_cv   <- solve(CV$B)        %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(1)
demand_cv <- cbind(time(var4$y)[-c(1:4)], demand_cv); colnames(demand_cv) <- c("time", "demand_cv")
supply_cv     <- solve(CV$B)          %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(2)
supply_cv <- cbind(time(var4$y)[-c(1:4)], supply_cv); colnames(supply_cv) <- c("time", "supply_cv")

non_mp_cv <- left_join(demand_cv, supply_cv, by = "time")

# RR_pro
temp <- left_join(non_mp_cv, RR.ts, by = "time") %>% tidyr::drop_na()
lm1 <-  lm(temp$RR ~ temp[,2:3] %>% as.matrix()-1)

RR_pro <- data.frame("time" = temp$time, "RR_pro" =  lm1 %>% resid())
RR_pro <- ts(RR_pro$RR_pro, frequency = 4, start = c(floor(RR_pro$time[1]), (RR_pro$time[1] %% 1) / 0.25 + 1 ))

IV.RR_pro <- get.id.iv(var4, instruments = RR_pro, Synchro = T)
IV.RR_pro$B
IV.RR_pro$F_test

#SZ_pro
temp <- left_join(non_mp_cv, SZ.ts, by = "time") %>% tidyr::drop_na()
lm2 <- lm(temp$SZ ~ temp[,2:3] %>% as.matrix()-1)

SZ_pro <- data.frame("time" = temp$time, "SZ_pro" =  lm2 %>% resid())
SZ_pro <- ts(SZ_pro$SZ_pro, frequency = 4, start = c(floor(SZ_pro$time[1]), (SZ_pro$time[1] %% 1) / 0.25 + 1 ))

IV.SZ_pro <- get.id.iv(var4, instruments = SZ_pro, Synchro = T)
IV.SZ_pro$B
IV.SZ_pro$F_test

# Cormat
new_cormat <- mp_all %>% left_join(IV.RR_pro$eps.ts, by = "time") %>% left_join(IV.SZ_pro$eps.ts, by = "time") 
colnames(new_cormat)[c(8,9)] <- c("RR_pro", "SZ_pro")
new_cormat %>% select(-1) %>% cor %>% round(3) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")

# shock plot
mp_all_pro <- mp_all %>% left_join(IV.RR_pro$eps.ts, by = "time") %>% left_join(IV.SZ_pro$eps.ts, by = "time")
colnames(mp_all_pro)[c(8,9)] <- c("RR_pro", "SZ_pro")

# including cv demand and supply shocks
mp_all_pro %<>% left_join(non_mp_cv, by = "time") 

mp_sds_pro <- mp_all_pro %>% select(-1) %>% apply( 2, sd)
mp_all_pro_stand <- mp_all_pro
for (j in 1:dim(mp_all_pro)[1]) {
  for (i in 2:ncol(mp_all_pro)) {
    mp_all_pro_stand[j,i] <- mp_all_pro[j,i] / mp_sds_pro[i-1]
  } 
}


# Event analysis ----------------------------------------------------------
# April 1974, October 1979, December 1988 and February 1994 - were contractionary
# December 1990, October 1998, April 2001, and November 2002  -were expansionary 
mp_pros <- IV.RR_pro$eps.ts %>% left_join(IV.SZ_pro$eps.ts, by = "time")
colnames(mp_pros)[-1] <- c("RR", "SZ")
mp_pro_tab <- mp_pros %>% filter(time == 1974.25 | time == 1979.75 | time == 1988.75 | time == 1994 | time == 1990.75 | time == 1998.75 | time == 2001.25| time == 2002.75)
#mp_all_stand_tab <- mp_all_stand_tab %>% select(-SR)
#mp_all_pro_stand <- mp_all_pro_stand %>% add_column(SUM = rowSums(mp_all_pro_stand[,-1]))

mp_pro_tab[c(4,5),] <- mp_pro_tab[c(5,4),]
mp_pro_tab %>% t() %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")



mps_RR_pro <- mp_all_pro_stand %>% select(time, old=RR, new=RR_pro) %>% melt(id.var = "time")
mps_SZ_pro <- mp_all_pro_stand %>% select(time, old=SZ, new=SZ_pro) %>% melt(id.var = "time")


mps_plot.SW <- ggplot(mp_all_pro_stand, aes(x = time, y = SW)) + geom_line() + ylim(-4.5,  4.5)+
  scale_x_continuous(breaks = seq(mp_all_pro_stand$time[1] %>% ceiling(),mp_all_pro_stand$time[nrow(mp_all)], 2)) +
  xlab("") + ylab(" ")  + theme_bw()  + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
  annotate("rect", xmin = 1979.5, xmax = 1983, ymin = -Inf, ymax = Inf, alpha = 0.3)

mps_plot.RR <- ggplot(mps_RR_pro, aes(x = time, y = value, group = variable)) + ylim(-4.5,  4.5)+
  geom_line(aes(linetype = variable)) + scale_linetype_manual(values = c("longdash", "solid")) +
  scale_x_continuous(breaks = seq(mp_all_pro_stand$time[1] %>% ceiling(),mp_all_pro_stand$time[nrow(mp_all)], 2)) +
  xlab("") + ylab(" ")  + theme_bw()  + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
  theme(legend.position = "none") + annotate("rect", xmin = 1979.5, xmax = 1983, ymin = -Inf, ymax = Inf, alpha = 0.3)

mps_plot.SZ <- ggplot(mps_SZ_pro, aes(x = time, y = value, group = variable)) + ylim(-4.5,  4.5)+
  geom_line(aes(linetype = variable)) + scale_linetype_manual(values = c("longdash", "solid")) +
  scale_x_continuous(breaks = seq(mp_all_pro_stand$time[1] %>% ceiling(),mp_all_pro_stand$time[nrow(mp_all)], 2)) +
  xlab("") + ylab(" ")  + theme_bw()  + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
  theme(legend.position = "none") + annotate("rect", xmin = 1979.5, xmax = 1983, ymin = -Inf, ymax = Inf, alpha = 0.3)

mps_plot_pros <- grid.arrange(mps_plot.SW, mps_plot.RR, mps_plot.SZ)


# Other Instruments -------------------------------------------------------
# Caldara and Herbst 
CH.raw      <- read.csv("../Data/Instruments/CH.csv")
CH          <- ts(CH.raw[,-1], frequency = 4, start = c(floor(CH.raw$Year[1]), (CH.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.CH       <- get.id.iv(var4, instruments = CH, Synchro = T)
IV.CH$B
IV.CH$F_test
Box.test(CH, lag = 1, type = "Ljung-Box")

# RR83b
RR83b.raw      <- read.csv("../Data/Instruments/RR83b.csv")
RR83b          <- ts(RR83b.raw[,-1], frequency = 4, start = c(floor(RR83b.raw$Year[1]), (RR83b.raw$Year[1] %% 1) / 0.25 + 1 ))
RR83b          <- Remove_AC(RR83b, lag = 1) 
IV.RR83b       <- get.id.iv(var4, instruments = RR83b, Synchro = T)
IV.RR83b$B
IV.RR83b$F_test

# FF1vr
FF1vr.raw      <- read.csv("../Data/Instruments/FF1vr.csv")
FF1vr          <- ts(FF1vr.raw[,-1], frequency = 4, start = c(floor(FF1vr.raw$Year[1]), (FF1vr.raw$Year[1] %% 1) / 0.25 + 1 ))
FF1vr          <- Remove_AC(FF1vr, lag = 1) 
IV.FF1vr       <- get.id.iv(var4, instruments = FF1vr, Synchro = T)
IV.FF1vr$B
IV.FF1vr$F_test

# FF4vr
FF4vr.raw      <- read.csv("../Data/Instruments/FF4vr.csv")
FF4vr          <- ts(FF4vr.raw[,-1], frequency = 4, start = c(floor(FF4vr.raw$Year[1]), (FF4vr.raw$Year[1] %% 1) / 0.25 + 1 ))
FF4vr          <- Remove_AC(FF4vr, lag = 1) 
IV.FF4vr       <- get.id.iv(var4, instruments = FF4vr, Synchro = T)
IV.FF4vr$B
IV.FF4vr$F_test

# FF1gb
FF1gb.raw      <- read.csv("../Data/Instruments/FF1gb.csv")
FF1gb          <- ts(FF1gb.raw[,-1], frequency = 4, start = c(floor(FF1gb.raw$Year[1]), (FF1gb.raw$Year[1] %% 1) / 0.25 + 1 ))
FF1gb          <- Remove_AC(FF1gb, lag = 1) 
IV.FF1gb       <- get.id.iv(var4, instruments = FF1gb, Synchro = T)
IV.FF1gb$B
IV.FF1gb$F_test

# FF4gb
FF4gb.raw      <- read.csv("../Data/Instruments/FF4gb.csv")
FF4gb          <- ts(FF4gb.raw[,-1], frequency = 4, start = c(floor(FF4gb.raw$Year[1]), (FF4gb.raw$Year[1] %% 1) / 0.25 + 1 ))
FF4gb          <- Remove_AC(FF4gb, lag = 1) 
IV.FF4gb       <- get.id.iv(var4, instruments = FF4gb, Synchro = T)
IV.FF4gb$B
IV.FF4gb$F_test

# ED2vr
ED2vr.raw      <- read.csv("../Data/Instruments/ED2vr.csv")
ED2vr          <- ts(ED2vr.raw[,-1], frequency = 4, start = c(floor(ED2vr.raw$Year[1]), (ED2vr.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.ED2vr       <- get.id.iv(var4, instruments = ED2vr, Synchro = T)
IV.ED2vr$B
IV.ED2vr$F_test


# ED2gb
ED2gb.raw      <- read.csv("../Data/Instruments/ED2gb.csv")
ED2gb          <- ts(ED2gb.raw[,-1], frequency = 4, start = c(floor(ED2gb.raw$Year[1]), (ED2gb.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.ED2gb       <- get.id.iv(var4, instruments = ED2gb, Synchro = T)
IV.ED2gb$B
IV.ED2gb$F_test


# MA
MA.raw      <- read.csv("../Data/Instruments/MA.csv")
MA          <- ts(MA.raw[,-1], frequency = 4, start = c(floor(MA.raw$Year[1]), (MA.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.MA       <- get.id.iv(var4, instruments = MA, Synchro = T)
IV.MA$B
IV.MA$F_test


