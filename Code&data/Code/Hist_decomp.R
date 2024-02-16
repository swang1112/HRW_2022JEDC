
Epsname = c("d", "s", "mp")

# pi ----------------------------------------------------------------------
hd_pi_SR <- hd.mod(SR, series = 2, Partial = NULL, Epsname)
hd_pi_SR_mp <- hd.mod(SR, series = 2, Partial = 3, Epsname)
hd_pi_CV <- hd.mod(CV, series = 2, Partial = NULL, Epsname)
hd_pi_CV_mp <- hd.mod(CV, series = 2, Partial = 3, Epsname)
hd_pi_dCov <- hd.mod(dCov, series = 2, Partial = NULL, Epsname)
hd_pi_dCov_mp <- hd.mod(dCov, series = 2, Partial = 3, Epsname)
hd_pi_SW_mp <- hd.mod(IV.SW, series = 2, Partial = 3, Epsname)
hd_pi_RR_mp <- hd.mod(IV.RR, series = 2, Partial = 3, Epsname)
hd_pi_SZ_mp <- hd.mod(IV.SZ, series = 2, Partial = 3, Epsname)


hd_pi_RR_mp_pro <- hd.mod(IV.RR_pro, series = 2, Partial = 3, Epsname)
hd_pi_SZ_mp_pro <- hd.mod(IV.SZ_pro, series = 2, Partial = 3, Epsname)


## combined plot
hd_pi_mp <- hd_pi_SR_mp %>% left_join(hd_pi_dCov_mp, by = "t") %>% 
  left_join(hd_pi_CV_mp, by = "t") %>% left_join(hd_pi_SW_mp, by = "t") %>% 
  left_join(hd_pi_SZ_mp, by = "t") %>% left_join(hd_pi_RR_mp, by = "t") %>% 
  left_join(hd_pi_SZ_mp_pro, by = "t") %>% left_join(hd_pi_RR_mp_pro, by = "t") %>% tidyr::drop_na()

hd_pi_mp %<>% dplyr::filter(t >= 1973 & t <= 1983) 

colnames(hd_pi_mp)[-1] <- c("SR+", "dCov", "CV", "SW", "SZ", "RR", "SZ_pro", "RR_pro")

recessions_start <- c( 1973.75, 1980,1981.5)
recessions_end <- c( 1975, 1980.5,1982.75)

hd_pi_mp_df = hd_pi_mp %>% reshape2::melt(id = "t") 
hd_pi_mp_df %<>% mutate(Label = ifelse(variable %in% c("SZ_pro", "RR_pro"), "pro", "orig"))

hd_pi_mp_df$variable[which(hd_pi_mp_df$variable  == 'SZ_pro')] = "SZ"
hd_pi_mp_df$variable[which(hd_pi_mp_df$variable  == 'RR_pro')] = "RR"

hd_pi_mp_plot <- hd_pi_mp_df %>% ggplot(aes(x = t, y = value, linetype = Label)) +
  geom_line() + facet_wrap(~variable, ncol = 3) + xlab("") + ylab("") + 
  theme_bw() + geom_hline(yintercept = 0, color = "red", alpha = 0.5) +
  scale_linetype_manual(values = c('solid', 'dashed')) + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(hd_pi_mp$t[1] %>% ceiling(),hd_pi_mp$t[nrow(hd_pi_mp)], by = 2)) +  
  annotate("rect", xmin = recessions_start, xmax = recessions_end, ymin = -Inf, ymax = Inf, alpha = 0.3)
hd_pi_mp_plot

