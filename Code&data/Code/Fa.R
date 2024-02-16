
# Information variables ---------------------------------------------------
current <- read_csv("../Data/Factors/fred-database_code/current.csv", na = "NA", col_types = c(.default = col_double(), sasdate = col_character()))
current <- current[-c(1, 3:24),] # Beginn 
current <- current[1:178,] # End
#View(current)

miss_val <- which(apply(current, 2, function(x){all(!is.na(x))}) == FALSE)
current %<>% select(-miss_val)
current[1,1] <- NA
current$sasdate  <- current$sasdate %>% as.numeric()
current[-1,1] <- seq(1964.5, 2008.5, 0.25)


# Data Transformation -----------------------------------------------------
Transfct <- function(x){
  code <- x[1]
  data <- x[-1]
  Tob <- length(data)
  if(code == 1) {
    erg <- data[-(1:2)]
  } else if(code == 2){
    erg_temp <- diff(data)
    erg <- erg_temp[-1]
  } else if(code == 3){
    erg <- diff(diff(data))
  } else if(code == 4){
    erg_temp <- log(data)
    erg <- erg_temp[-(1:2)]
  } else if(code == 5){
    erg_temp <- diff(log(data))
    erg <- erg_temp[-1]
  } else if(code == 6){
    erg <- diff(diff(log(data)))
  } else if(code == 7){ 
    erg_temp <- rep(NA, Tob-1)
    for (i in 1:(Tob-1)) {
      erg_temp[i] <- x[i+1]/x[i] - 1
    }
    erg <- diff(erg_temp)
  } else {stop("Invalid transformation code!")}
  return(erg)
}


# Principal Component Estimation ------------------------------------------

Mccr <- cbind(current[-(1:3),1],apply(current[,-1], 2, Transfct))

pcfac <- function(data, r){
  
  ## PC factor extraction
  ## data must be (TxN)  
  nvars <- dim(data)[2]
  
  xx <- t(data)%*%data# X'X is proportional to covariance matrix
  eigen_x <- eigen(xx) 
  
  # applying the factor restricition
  lam <- sqrt(nvars)*eigen_x[[2]][,c(1:r)] # lam'lam/n = I. loadings are, where n = number of time series in x
  factors <- data%*%lam/nvars     
  
  return(list("Fact" = factors, "Lambda" = lam))
}

PCs <- pcfac(Mccr[,-1] %>% as.matrix(), r)

# Import data -------------------------------------------------------------
if (!SDFM){
  if (is_in("USA_Tri", ls())){
    Fact_orth <- PCs$Fact - USA_Tri %*% solve(t(USA_Tri) %*% USA_Tri) %*% t(USA_Tri) %*% PCs$Fact
    USA_FAVAR <- ts(cbind(USA_Tri, Fact_orth), frequency = 4, start = c(floor(USA_Tri_raw$Year[1]), (USA_Tri_raw$Year[1] %% 1) / 0.25 + 1 ))
    VARselect(USA_FAVAR, lag.max = 8) %>% print()
  } else if (is_in("USA_Qua", ls())){
    Fact_orth <- PCs$Fact - USA_Qua %*% solve(t(USA_Qua) %*% USA_Qua) %*% t(USA_Qua) %*% PCs$Fact
    USA_FAVAR <- ts(cbind(USA_Qua, Fact_orth), frequency = 4, start = c(floor(USA_Qua.raw$Year[1]), (USA_Qua.raw$Year[1] %% 1) / 0.25 + 1 ))
    VARselect(USA_FAVAR, lag.max = 8) %>% print()
  } else {stop("undefined VAR!") }

} else {
  Fact <- ts(PCs$Fact, frequency = 4, start = c(floor(USA_Tri_raw$Year[1]), (USA_Tri_raw$Year[1] %% 1) / 0.25 + 1 ))
  colnames(Fact) <- 1:r
  VARselect(Fact, lag.max = 8, type = "none") %>% print()
}
rm(list = c("current", "miss_val", "Mccr"))

