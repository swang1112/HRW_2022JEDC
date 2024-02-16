# Load functions ----------------------------------------------------------
rm(list = ls())
Functions <- devtools::as.package("../Functions")
devtools::load_all("../Functions")

IDM   <- c("SR0", "SR1", "CV", "dCov", "CvM", "IV", "IV_W", "IV_P")
Lab   <- c("SR",  "SR+", "CV", "dCov", "CvM", "IV", "IV-W", "IV-E")
Dist  <- c("L1", "L2a", "L2b", "L3a", "L3b")

# Graphics ----------------------------------------------------------------
Tobs <- 240
for (i in 1:5) {
  aoa_temp <- list()
  for (j in 1:length(IDM)) {
    Path_In <- paste0("Server/", "T", Tobs, "/", Dist[i], "/return/", IDM[j], ".rds")
    temp <- readRDS(Path_In)
    if (IDM[j] == "SR0" || IDM[j] == "SR1" ){
      aoa_temp[[j]] <- temp$M
    } else {
      aoa_temp[[j]] <- temp$aoa
    }
  }
  plot_temp <- plot.simu.multi(aoa_temp, Label = Lab, shape_manual = c(5,1,2, 15:17,3,4), level_manual = c("IV", "IV-W", "IV-E", "CV", "CvM", "dCov", "SR",  "SR+"))
  Path_Out <- paste0("../Plots/")
  ggsave(paste0(Path_Out, "T", Tobs, "_", Dist[i], ".pdf"), plot = plot_temp, width = 12, height = 8)
  rm(temp, aoa_temp, plot_temp)
}

# Tables ------------------------------------------------------------------
TT <- c(120, 240, 480)

# remove leading zero for frequencies
.nozero = function(x) substring(sprintf("%4.3f", x), 2)

## AoA tables
aoa_s <- matrix(NA, nrow = length(IDM)*5, ncol = 12)
aoa_l <- matrix(NA, nrow = length(IDM)*5, ncol = 12)
for (i in 1:5) {
  for (t in 1:3) {
    aoa_temp <- list()
    for (j in 1:length(IDM)) {
      Path_In <- paste0("Server/", "T", TT[t], "/", Dist[i], "/return/", IDM[j], ".rds")
      temp <- readRDS(Path_In)
      if (IDM[j] == "SR0" || IDM[j] == "SR1" ){
        aoa_temp[[j]] <- temp$M
      } else {
        aoa_temp[[j]] <- temp$aoa
      }
    }
    tab_temp <- tabulate.simu.aoa(aoa_temp, Label = Lab) 
    aoa_s[c((length(IDM)*(i-1)+1):(i*length(IDM))),c((4*(t-1)+1):(t*4))] <- tab_temp$short
    aoa_l[c((length(IDM)*(i-1)+1):(i*length(IDM))),c((4*(t-1)+1):(t*4))] <- tab_temp$long
  }
  rm(temp, aoa_temp, tab_temp)
}
rownames(aoa_s) <- rownames(aoa_l) <- rep(Lab, 5)  

# aoa_s %>% .nozero

# i = 1
# aoa_s[c((8*i-2):(8*i), (8*(i-1)+1):(8*(i-1)+5)),]

kableExtra::kable(aoa_s, format = "latex", col.names = 1:12, booktabs = T, digits = 3, linesep = "") %>% kableExtra::row_spec(c(1:5)*length(IDM), hline_after = T)
kableExtra::kable(aoa_l, format = "latex", col.names = 1:12, booktabs = T, digits = 3, linesep = "") %>% kableExtra::row_spec(c(1:5)*length(IDM), hline_after = T)

## UMP table
ump_tab <- matrix(NA, nrow = length(IDM), ncol = 15)
for (i in 1:5) {
  for (t in 1:3) {
    ump_temp <- list()
    for (j in 1:length(IDM)) {
      Path_In <- paste0("Server/", "T", TT[t], "/", Dist[i], "/return/", IDM[j], ".rds")
      ump_temp[[j]] <- readRDS(Path_In)
    }
    ump_tab[, ((i-1)*3+t)] <- tabulate.simu.ump(ump_temp, Lab)$ump
  }
  rm(ump_temp)
}
rownames(ump_tab) <- Lab
kableExtra::kable(ump_tab, format = "latex", col.names = 1:15, booktabs = T, digits = 3, linesep = "")

## RMSE table
rmse_tab <- matrix(NA, nrow = length(IDM)*5, ncol = 12)

for (i in 1:5) {
  for (t in 1:3) {
    rmse_temp <- list()
    for (j in 1:length(IDM)) {
      Path_In <- paste0("Server/", "T", TT[t], "/", Dist[i], "/return/", IDM[j], ".rds")
      rmse_temp[[j]] <- readRDS(Path_In)
    }
    rmse_tab[c((length(IDM)*(i-1)+1):(i*length(IDM))),c((4*(t-1)+1):(t*4))] <- tabulate.simu.ump(rmse_temp, Lab)$rmse
  }
  rm(rmse_temp)
}
rownames(rmse_tab) <- rep(Lab, 5)
kableExtra::kable(rmse_tab, format = "latex", col.names = 1:12, booktabs = T, digits = 3, linesep = "") %>% kableExtra::row_spec(c(1:5)*length(IDM), hline_after = T)

# UMPS Plot ---------------------------------------------------------------
ump_plot_int <- as.data.frame(ump_tab)
ump_plot_int["SR",] <- NA
ump_plot_int["SR+",] <- NA
ump_plot_int <- t(ump_plot_int) %>% as.data.frame()
ump_plot_int <- cbind("Dist" = rep(Dist, each = 3), ump_plot_int)  %>% reshape2::melt(id = "Dist")
ump_plot_int <- ump_plot_int %>% add_column(Tob = rep(TT, 5*8))
ump_plot_int <- ump_plot_int[,c(1,2,4,3)]
ump_plot_int <- ump_plot_int %>% add_column(measure = "Frequency")

mean_rmse <- rmse_tab[,c(1,5,9)] %>% as.data.frame()
rmse_plot_int <- mean_rmse[1:8,]
for (i in 2:5) {
  rmse_plot_int <- cbind(rmse_plot_int, mean_rmse[(8*(i-1)+1):(8*i),])
}
rownames(rmse_plot_int) <- Lab
#rmse_plot_int["IV-P",] <- NA
rmse_plot_int <- t(rmse_plot_int) %>% as.data.frame()
rmse_plot_int <- cbind("Dist" = rep(Dist, each = 3), rmse_plot_int)  %>% reshape2::melt(id = "Dist")
rmse_plot_int <- rmse_plot_int %>% add_column(Tob = rep(TT, 5*8))
#rmse_plot_int <- rmse_plot_int[,c(1,2,4,3)]
rmse_plot_int <- rmse_plot_int %>% add_column(measure = "RMSE")

UMPS_plot <- rbind(ump_plot_int, rmse_plot_int)


library(latex2exp)
levels(UMPS_plot$Dist) <- c(TeX("$\\textbf{L}_I$"), TeX("$\\textbf{L}^{(a)}_{II}$"), 
                            TeX("$\\textbf{L}^{(b)}_{II}$"), TeX("$\\textbf{L}^{(a)}_{III}$"),
                            TeX("$\\textbf{L}^{(b)}_{III}$"))

colnames(UMPS_plot)[2] <- "Methods"

UMPS_plot$measure<- UMPS_plot$measure %>% as.factor()
levels(UMPS_plot$measure)[2] <- TeX("$\\bar{RMSE}$")

UMPS_plot$Methods <- UMPS_plot$Methods %>% as.character()
UMPS_plot$Tob <- UMPS_plot$Tob %>% as.character()

UMPS_plot <- UMPS_plot %>% mutate(Methods = replace(Methods, Methods == "IV-P", "IV-E"))

UMPS_plot$Methods <- factor(  UMPS_plot$Methods, levels = c("IV", "IV-W", "IV-E", "CV", "CvM", "dCov", "SR", "SR+"))


ggplot(data = UMPS_plot, aes(x = Tob, y = value, group = Methods)) + geom_line(aes(linetype = Methods, color = Methods)) +
  geom_point(aes(shape = Methods, color = Methods), size=3) +
  facet_grid(measure~Dist, scales = "free_y", labeller = label_parsed) +
  xlab(" T ") + ylab(" ")  + theme_bw() + #scale_shape_manual(values = c(15:17, 5, 1:3)) +
  scale_shape_manual(values = c(5,1,2, 15:17,3,4))

ggsave("../Plots/UMPS_all_22.pdf", plot = last_plot(), width = 12, height = 7)



