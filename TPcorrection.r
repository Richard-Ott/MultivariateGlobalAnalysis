library(MASS)
library(ggplot2)
library(RColorBrewer)
setwd("E:/Richard/Global_Hypsometry/GIS/sinusoidal/T_P_correction")

# NORMAL CORRECTION ###########################################
Cvars <- c("P",'T','tetrapods')
Cdata <- read.table("BID_for_lm.txt", header = FALSE, sep = ",", col.names	= Cvars)  

# Tetrapods
ClinMod <- lm(Cdata$tetrapods ~ Cdata$T + Cdata$P, data= Cdata)
Clm_results <- summary(ClinMod)
# Empirical expected tetrapod richness at every T/P pair
CTETemp <- Clm_results$coefficients[[1,1]] + Cdata$T * Clm_results$coefficients[[2,1]] + Cdata$P * Clm_results$coefficients[[3,1]]
# correct tetrapods
CTETcorr <- Cdata$tetrapods - CTETemp

# write(CTETcorr, file = "BIDcorr_for_matlab.txt", sep = ",")
# write.table(CTETcorr, "BIDcorr_for_matlab.txt", append = FALSE, sep = ",", dec = ".",
#             row.names = F, col.names = F)

# THRESHOLD CORRECTION ########################################
###############################################################

Cvars <- c("P",'T','tetrapods')
Cdata <- read.table("MAM_for_lm.txt", header = FALSE, sep = ",", col.names	= Cvars)  

Cdata$P.thres = Cdata$P
Cdata$P.thres[Cdata$P > 3500] = 3500
# Tetrapods
ClinModT <- lm(Cdata$tetrapods ~ Cdata$T + Cdata$P.thres, data= Cdata)
Clm_resultsT <- summary(ClinModT)
# Empirical expected tetrapod richness at every T/P pair
CTETempT <- Clm_resultsT$coefficients[[1,1]] + Cdata$T * Clm_resultsT$coefficients[[2,1]] + Cdata$P.thres * Clm_resultsT$coefficients[[3,1]]
# correct tetrapods
CTETcorrT <- Cdata$tetrapods - CTETempT

# write(CTETcorrT, file = "BIDcorrThreshold_for_matlab.txt", sep = ",")
# write.table(CTETcorrT, "MAMcorrThreshold_for_matlab.txt", append = FALSE, sep = ",", dec = ".",
#             row.names = F, col.names = F)


# PLOT BI-RELATIONSSHIPS ######################################
# T- BID ######################################################
TF <-ggplot(data = Cdata, aes(x=Cdata$T, y=Cdata$tetrapods)) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
    geom_bin2d(bins = 100)+geom_contour()
    scale_fill_distiller(palette= "Spectral", direction=1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
    theme_bw() 
TF + labs(title = "T-correlation",x = "T",y= "Tetrapod richness")

# P- BID ######################################################
P <-ggplot(data = Cdata, aes(x=Cdata$P, y=Cdata$tetrapods)) +
  geom_density_2d()
  # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
  # geom_bin2d(bins = 200) + geom_contour(data = Cdata)+
  # stat_density_2d(aes(fill = ..density..), geom = "raster", contour = T) +
  # scale_fill_distiller(palette= "Spectral", direction=1) +
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0)) +
  theme_bw() 
P + labs(title = "P-correlation",x = "P",y= "Tetrapod richness")


