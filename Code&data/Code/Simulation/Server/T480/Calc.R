rm(list = ls())
start_time = Sys.time()

# Server config -----------------------------------------------------------
server_mode = T
cl          = 32
cl_sub      = 16

Tob <- 480
AOA <- readRDS("../AOA.rds")

# Load functions ----------------------------------------------------------
Functions <- devtools::as.package("../../../Functions")
devtools::load_all("../../../Functions")
source("../../../aux_funs/sign_simu1.R")

# Distributional settings -------------------------------------------------
Dist      <- c("normdist", "tdist", "tdist", "mixed_b", "mixed_a")
Dist_lab  <- c("L1", "L2a", "L2b", "L3a", "L3b")

# Load configurations -----------------------------------------------------
for (l in 1:5) {

  path.temp   <- paste0(Dist_lab[l], "/Config.rds")
  Config.temp <- readRDS(path.temp)
  
  # Data generating ---------------------------------------------------------
  set.seed(1234)
  Y <- simu.datas(Config = Config.temp, DGP = 1)

  # Identification and evaluation -------------------------------------------
  ## Sign restrictions
  Signmat   <- matrix(c(1,1,1,-1,1,1,NA,-1,1), nrow = 3, ncol = 3) # agnostic sign restrictions
  TrueSign  <- matrix(c(1,1,1,-1,1,1,-1,-1,1), nrow = 3, ncol = 3) # true signs
  set.seed(1234)
  # SR        <- simulate.sign(Data = Y, AOA = AOA, Iter = 1000, Itermax = Inf, Signmat = Signmat, TrueSign = TrueSign, Bmat = Config.temp$Par$B, RHorizont = 0, MT = T, Ci = 0.68)
  SR0        <- simulate.sign1(Data = Y, AOA = AOA, Iter = 1000, Signmat = Signmat, TrueSign = TrueSign, Bmat = Config.temp$Par$B, Arias = F, cl = cl, cl_sub = cl_sub)
  SR1        <- simulate.sign1(Data = Y, AOA = AOA, Iter = 1000, Signmat = Signmat, TrueSign = TrueSign, Bmat = Config.temp$Par$B, Arias = T, cl = cl, cl_sub = cl_sub)
  
  saveRDS(SR0, paste0(Dist_lab[l], "/return/SR0.rds"))
  saveRDS(SR1, paste0(Dist_lab[l], "/return/SR1.rds"))
  
  
  ## Identification by change in volatility
  #set.seed(1234)
  #CV        <- simulate.cv(Data = Y, AOA = AOA, Break = Config.temp$Break, Bmat = Config.temp$Par$B, Guess = TRUE)
  #saveRDS(CV, paste0(Dist_lab[l], "/return/CV.rds"))
  
  ## Identification by minimizing distance covariance
  #set.seed(1234)
  #dCov      <- simulate.dcov(Data = Y, AOA = AOA, Bmat = Config.temp$Par$B)
  #saveRDS(dCov, paste0(Dist_lab[l], "/return/dCov.rds"))
  
  ## Identification by valid instruments
  #set.seed(1234)
  #IV        <- simulate.iv(Data = Y, AOA = AOA, Bmat = Config.temp$Par$B)
  #saveRDS(IV, paste0(Dist_lab[l], "/return/IV.rds"))
  
  ## Identification by weak instruments
  #set.seed(1234)
  #IV.W      <- simulate.iv(Data = Y, AOA = AOA, Weak = T, Bmat = Config.temp$Par$B)
  #saveRDS(IV.W, paste0(Dist_lab[l], "/return/IV_W.rds"))
  
  ## Identification by invalid instruments
  #set.seed(1234)
  #IV.P      <- simulate.iv(Data = Y, AOA = AOA, Combi = T, Bmat = Config.temp$Par$B)
  #saveRDS(IV.P, paste0(Dist_lab[l], "/return/IV_P.rds"))
  
  ## Identification by minimizing Cramer-von Mises distance 
  # set.seed(1234)
  # CvM       <- simulate.cvm(Data = Y, AOA = AOA, Iter = 750, Steptol = 200, Stage2 = 200, Bmat = Config.temp$Par$B)
  # saveRDS(CvM, paste0(Dist_lab[l], "/return/CvM.rds"))
}

end_time = Sys.time()
elapsed_time = difftime(end_time, start_time)
Journal = data.frame(start_time, end_time, elapsed_time)
write.csv(Journal, "journal.csv")
