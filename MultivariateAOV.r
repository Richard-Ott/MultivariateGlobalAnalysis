# My very first R Script!!!!!!
# Richard Ott, 2019
setwd("E:/Richard/Global_Hypsometry/matlab/global correlations/200m_relief")
rm(list = ls())
library(DescTools)
library(nnet)
library(graphics)
library(car)
library(ggplot2)
library(plyr)
library(nlme)


# LOAD DATA ############################################################################
vars <- c("Elev","Locrel","slope","P",'T','Pse','NDVI','Lat','ksn','tetrapods'
          ,"tet_corr",'Amphibians',"Amph_corr",'Mammals',"Mam_corr",'KG','Geo','LC')
dat <- read.table("global_data_200mLC.txt", header = FALSE, sep = ",", col.names	= vars)
vars[19] <- "NDVIc"
dat$Amphibians[dat$Amphibians == -9999 | dat$Amphibians == -32768] = 0   # Amph have lots of nan values that should be set to zero for analysis
dat$Amph_corr[dat$Amph_corr == -9999 | dat$Amph_corr == -32768] = 0   # Amph have lots of nan values that should be set to zero for analysis
dat$tet_corr[dat$tet_corr == -9999 | dat$tet_corr == -32768] = 0 
dat$Mam_corr[dat$Mam_corr == -9999 | dat$Mam_corr == -32768] = 0 
corr = 1 # Do you want biological variables to be corrected for T and P influence?

# MODIFY DATA ########################################################################
depVar <- which(vars == 'Geo')                     # dependent variable, should be last variable
#         E Loc Sl P  T Pse  ND L ksn te teC am  amC ma maC KG Geo LC  NDVIc
inds <- c(T, T, T, T, F, F,  F, T, T, F , T ,F ,  T,  F, T,  F, F , T ,  T)

GeoIDf <- vector()
GeoIDf[dat[,depVar] == 1] <- 1                      # su
GeoIDf[is.element(dat[,depVar],c(3,7,8,9))] <- 2    # vc
GeoIDf[is.element(dat[,depVar],c(10,11,12))] <- 3   # pl
GeoIDf[dat[,depVar] == 2] <- 4                      # ss
GeoIDf[dat[,depVar] == 4] <- 5                      # sm
GeoIDf[dat[,depVar] == 5] <- 6                      # sc
GeoIDf[dat[,depVar] == 13] <- 7                     # mt
GeoIDf <- factor(GeoIDf)
g <- length(unique(GeoIDf))
g <- 7
lab = c("su","vc","pl","ss","sm","sc","mt")

########################################################################################
# TEMPERATURE AND PRECIPITATION CORRECTION OF BIOLOGY ##################################
if (corr == 1){
  thrs <- 2500
  # NDVI #######################
  linMod <- lm(dat$NDVI ~ dat$T + dat$P, data= dat)
  lm_results <- summary(linMod)
  # Empirical expected NDVI at every T/P pair with 3500mm threshold of rain
  Pthres = dat$P
  Pthres[Pthres > thrs] = thrs
  # NDVIemp <- lm_results$coefficients[[1,1]] + Pthres * lm_results$coefficients[[2,1]]
  NDVIemp <- lm_results$coefficients[[1,1]] + dat$T * lm_results$coefficients[[2,1]] + Pthres * lm_results$coefficients[[3,1]]
  # correct NDVI
  NDVIcorr <- dat$NDVI - NDVIemp
  dat$NDVIc <- NDVIcorr
  
  # AMPHiIBIANS #######################
  linMod <- lm(dat$Amphibians ~ dat$T + dat$P, data= dat)
  lm_results <- summary(linMod)
  Pthres = dat$P
  Pthres[Pthres > thrs] = thrs
  Amph.emp <- lm_results$coefficients[[1,1]] + dat$T * lm_results$coefficients[[2,1]] + Pthres * lm_results$coefficients[[3,1]]
  Amph.corr <- dat$Amphibians - Amph.emp
  dat$Amph_corr <- Amph.corr
  
  # TETRAPODS #######################
  linMod <- lm(dat$tetrapods ~ dat$T + dat$P, data= dat)
  lm_results <- summary(linMod)
  Pthres = dat$P
  Pthres[Pthres > thrs] = thrs
  Tet.emp <- lm_results$coefficients[[1,1]] + dat$T * lm_results$coefficients[[2,1]] + Pthres * lm_results$coefficients[[3,1]]
  Tet.corr <- dat$Tetrapods - Tet.emp
  dat$tet_corr <- Tet.corr
  
  rm(linMod,lm_results)
}

# group and assign factor for landcover data ###########################################
lcID <- vector()
lcID[is.element(dat$LC,c(1,2,14))] <-  1 # bFor, broadleaf forest
lcID[is.element(dat$LC,c(3,4))] <-  2       # nFor, needleleaf forest
#lcID[dat$LC == 5] <- 3                      # mFor, mixed forest
lcID[is.element(dat$LC,c(6,7,8,9))] <-  3   # Medium Vegetation
lcID[is.element(dat$LC,c(10,16,17))] <-  4  # Sparse/no Vegetation
lcID[is.element(dat$LC,c(11,12))] <-  5     # Cropland, 13 seems to be likely cropland but it seemed weird
lcID[dat$LC == 18] <- 6                     # Urban
#lcID[dat$LC == 19] <- 8                     # Snow/Ice, so small can be left out
#lcID[dat$LC == 20] <- 9                     # Water, so small can be left out
lcID <- factor(lcID)
# LClab <- c("bFor","nFor","mFor","medVeg","noVeg","Crop","Urban","Snow","Water")
LClab <- c("bFor","nFor","medVeg","noVeg","Crop","Urban")


# STATISTICS OF DATA-SET ###############################################################
stats <- summary(dat)  # entire data-set
stats.by.group <- array(NA,dim=c(length(vars),g*2))
for (i in 1:length(vars)){
  stats.by.group[i,1:g]<- tapply(dat[,i],GeoIDf,mean)
  stats.by.group[i,seq(g+1,g*2)]<- tapply(dat[,i],GeoIDf,median)
}
rownames(stats.by.group) <- vars
colnames(stats.by.group) <- c(lab,lab)
write.table(stats.by.group, "MeanStats_200m.csv", append = FALSE, sep = ",", dec = ".",
            row.names = T, col.names = T)
# Correlation coefficicient matrix
corrcoefs <- cor(dat,use=)

##################################################################
## PLOT DATA #####################################################


# BOXPLOT OF ALL DATA #################################################################
#######################################################################################
dev.new()
setEPS()
postscript("boxplot_all_corrected.eps")
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(3,4),cex=0.35,lwd=0.35)
llim <- c(0,    0,  0, 300, 0,  -150, 0,  -100, -15,  -100, -15)     # corrected limits
ulim <- c(2000,400,120,1500,100, 500, 250, 500, 3,   150, 10)
cols = c("#FFEBAF","#38A800","#727272","#FFAA00","#BED2FF","#0070FF","#734C00","#000000")
# dev.set(which == 2)
for (i in 1:length(vars[inds])){ # loop through variables
  # rearrange data matrix
  data.by.geo <-list() 
  for (j in 1:g){               # loop through lithologies to rearrange data
    data.by.geo[[j]] <- dat[,inds][GeoIDf == j,i]
  }
  
  # boxplot(data.by.geo, xaxt = "n", main = vars[inds][i], outline = F, boxlty = 0, whisklty = 0, staplelty = 0, ylim = c(llim[i],ulim[i]))
  # boxplot(data.by.geo, xaxt = "n", main = vars[inds][i])
  
  boxplot(data.by.geo, xaxt = "n", main = vars[inds][i], outline = F, staplelty = 0, whisklty = 1, col = cols)
  axis(1, at=1:length(lab), labels=lab)
}
dev.copy(which = 4)
dev.off(which = 4)
dev.off()

# # VIOLIN PLOT OF ALL DATA ############################################################
# ######################################################################################
# 
# ind = 1:length(vars[inds])
# uplims = c(4e3,0.6e3,2e2,5e3,250,255,600,500,600,60,100)
# lolims = c(0,0,0,0,0,0,-550,0,0,0,0)
# dev.new()
# for (i in 1:length(vars[inds])){
#   g1 <- data.frame(dat[GeoIDf == 1,ind[i]],variable = lab[1])
#   g2 <- data.frame(dat[GeoIDf == 2,ind[i]],variable = lab[2])
#   g3 <- data.frame(dat[GeoIDf == 3,ind[i]],variable = lab[3])
#   g4 <- data.frame(dat[GeoIDf == 4,ind[i]],variable = lab[4])
#   g5 <- data.frame(dat[GeoIDf == 5,ind[i]],variable = lab[5])
#   g6 <- data.frame(dat[GeoIDf == 6,ind[i]],variable = lab[6])
#   g7 <- data.frame(dat[GeoIDf == 7,ind[i]],variable = lab[7])
#   names(g2) <- names(g1) 
#   names(g3) <- names(g1)
#   names(g4) <- names(g1)
#   names(g5) <- names(g1)
#   names(g6) <- names(g1)
#   names(g7) <- names(g1)
#   dat.by.geo <- rbind(g1,g2,g3,g4,g5,g6,g7)
#   rm(g1,g2,g3,g4,g5,g6,g7)
#   names(dat.by.geo)[names(dat.by.geo)=="dat.GeoIDf....1..ind.i.."] <- "val"
#   
# 
#   setEPS()
#   postscript(paste("violinplot_", vars[inds][i] ,".eps",sep=""))
#   dev.set(which = 2)
#   if (i == 10){
#     P<-ggplot(dat.by.geo, aes(variable, val)) + geom_violin(aes(fill = val),adjust = 2,trim = T,scale = "width") + labs(main =vars[inds][i], y = vars[inds][i]) + ylim(lolims[i],uplims[i])
#   }else{
#     P<-ggplot(dat.by.geo, aes(variable, val)) + geom_violin(aes(fill = val),trim = T,scale= "width") + labs(main =vars[inds][i], y = vars[inds][i]) + ylim(lolims[i],uplims[i])
#   }
#   # rm(dat.by.geo)
#   P <- P + geom_boxplot(width=0.05, outlier.shape = NA)
#   P <- P+ scale_x_discrete(labels = lab) + theme(panel.background = element_blank())
#   print(P)
#   dev.copy(which = 4)
#   dev.off(which = 4)
# }


# boxplot of data per climate zone ###################################################
######################################################################################
# the next two lines change the NDVIcorrected to NDVI for climate zone plots, uncomment
# if not wished
inds[7] <- T
inds[length(inds)] <- F
# chnage biodiversity corrected to normal
inds[c(6,8,11,13,15)] <- F
inds[c(10,12,14)] <- T
# inds[c(6,8)] <- F
# dev.new()
KGlab = c("tropical", "dry", "temp_cont","polar")
for (k in 1:4){
  # setEPS()
  # postscript(paste("boxplotKG_", KGlab[k] ,".eps", sep  = "", collapse = NULL))
  # dev.set(which = 2)
  # par(oma=c(0,0,2,0),mfrow=c(3,4),cex=0.5) 
  
  if (k == 3){   # This makes climate zones temperate and continental merge together
    dat.kg <-dat[dat$KG == 3 | dat$KG == 4 ,inds]
    GeoKG <- GeoIDf[dat$KG == 3 | dat$KG == 4]
    stat.kg <- dat[dat$KG == 3 | dat$KG == 4 ,]
  }else{
    dat.kg <-dat[dat$KG == k,inds]
    GeoKG <- GeoIDf[dat$KG == k]
    stat.kg <- dat[dat$KG == k,]
  }
  
  # do statistics on this KG ########################
  stats.by.group <- array(NA,dim=c(length(vars),g*2))
  for (i in 1:length(vars)){
    stats.by.group[i,1:g]<- tapply(stat.kg[,i],GeoKG,mean)
    stats.by.group[i,seq(g+1,g*2)]<- tapply(stat.kg[,i],GeoKG,median)
  }
  rownames(stats.by.group) <- vars
  colnames(stats.by.group) <- c(lab,lab)
  write.table(stats.by.group, paste("MeanKG_", KGlab[k] ,".csv", sep  = "", collapse = NULL), append = FALSE, sep = ",", dec = ".",
              row.names = T, col.names = T)
  ####################################################
  
  # for (i in 1:length(vars[inds])){
  #   # rearrange data matrix
  #   areas <- vector()
  #   dat.by.geo <-list()
  #   for (j in 1:g){
  #     Ginds = GeoKG == j
  #     dat.by.geo[[j]] <- dat.kg[Ginds,i]
  #     areas[j] <- length(Ginds[Ginds == T])/nrow(dat.kg)
  #   }
  #   boxplot(dat.by.geo, xaxt = "n", main = vars[inds][i], outline = F, staplelty = 0, whisklty = 1)
  #   axis(1, at=1:g, labels=lab)
  # }
  # pie(areas, labels = lab)
  # title(main= KGlab[k] , outer=T)
  # # dev.copy(which = 4)
  # # dev.off(which = 4)
}
rm(Ginds,dat.kg,areas,GeoKG)
# violinplot of data per climate zone ###################################################
#########################################################################################
# 
# KGlab = c("tropical", "dry", "temp_cont","polar")
# uplims = c(4e3,0.6e3,2e2,5e3,250,255,80,500,600,50,200)
# lolims = c(0,0,0,0,0,0,-70,0,0,0,0)
# for (k in 1:4){
#   setEPS()
#   postscript(paste("violinKG_", KGlab[k] ,".eps", sep  = "", collapse = NULL))
#   dev.set(which = 2)
#   par(oma=c(0,0,2,0),mfrow=c(3,4),cex=0.5)
#   
#   if (k == 3){
#     dat.kg <-dat[dat$KG == 3 | dat$KG == 4 ,inds]
#     GeoKG <- GeoIDf[dat$KG == 3 | dat$KG == 4]
#   }else{
#     dat.kg <-dat[dat$KG == k,inds]
#     GeoKG <- GeoIDf[dat$KG == k]
#   }
#   for (i in 1:length(vars[inds])){
#     # rearrange data matrix
#     dat.by.geo <-list()
# 
#     for (j in 1:g){
#       dat.by.geo[[j]] <- dat.kg[GeoKG == j,i]
#     }
#     if (i == 10){
#       P<-ggplot(dat.by.geo, aes(variable, val)) + geom_violin(aes(fill = val),adjust = 2,trim = T,scale = "width") + labs(main =vars[inds][i], y = vars[inds][i]) + ylim(lolims[i],uplims[i])
#     }else{
#       P<-ggplot(dat.by.geo, aes(variable, val)) + geom_violin(aes(fill = val),trim = T,scale= "width") + labs(main =vars[inds][i], y = vars[inds][i]) + ylim(lolims[i],uplims[i])
#     }
#     P <- P + geom_boxplot(width=0.05, outlier.shape = NA)
#     P <- P+ scale_x_discrete(labels = lab) + theme(panel.background = element_blank())
#     print(P)
#   }
#   title(main= KGlab[k] , outer=T)
#   dev.copy(which = 4)
#   dev.off(which = 4)
# }

# MEDIAN PLOT PER CLIMATE ZONE ################################################################
###############################################################################################
library(ggsci)
polar = F                       # do you want polar values to be displayed?
pp = 3 # polar plot -pp
if (polar){
  pp = 4
}
KGlab = c("tropical", "dry", "temp_cont","polar")
#          E Loc Sl P  T Pse  ND L ksn te teC am  amC ma maC KG Geo  LC   NDVIc
pinds <- c(T, F, T, F, F, F,  F, F, T, F , T ,F ,  T,  F, F,  F, F , F ,   T)  #plot indices
dev.new()
setEPS()
postscript(paste("boxplotKG_Medians.eps", sep  = "", collapse = NULL))
dev.set(which = 2)
par(oma=c(0,0,0,0),mfrow=c(2,3),cex=0.8)
fKG <- as.factor(dat$KG)
levels(fKG) <- c("1","2","3","3","4")
cols = c("#FFEBAF","#38A800","#727272","#FFAA00","#BED2FF","#0070FF","#734C00","#000000")
plotchar <- c(seq(16,18),15,15,19,17)
# llim <- c(500,100,40,0,50,100,150,0,30)    # uncorrected limits
# ulim <- c(4000,650,160,2200,165,350,600,50,150)   
# llim <- c(500,100,40,0,50,-100,-15,-20,-30)     # corrected limits
# ulim <- c(4000,650,160,2200,300,100,10,20,30)
llim <- c(400, 40, 140  ,-210,-8, -40)
ulim <- c(2400,120,230,   700, 8  ,21)

for (i in 1:length(vars[pinds])){
  mdata <- dat[,which(pinds)[i]]   # data for this variable
  msplit <- split(mdata,GeoIDf)   # split data by Geo
  mKG <- split(fKG,GeoIDf)        # split KG data by Geo
  mean.Geo.KG = array(-9999,dim=c(g,length(KGlab)))
  # set up empty plot
  plot(c(0.8,pp+0.2), c( llim[i], ulim[i]) ,type="n",ylab = vars[pinds][i], xlab="")
  axis(1, at=1:pp, labels=KGlab[seq(1,pp)])
  for (k in 1:g){                 # loop through lithologies
    mk <- msplit[[k]]               # take a certain lithology
    mkKG <- mKG[[k]]                 # take its KG distribution
    mGeoKG.split <- split(mk,mkKG) # split by kG
    
    dummy <- lapply(mGeoKG.split, median, simplify = TRUE) # do means of this lithology in different KG's
    
    mean.Geo.KG[k,] <- simplify2array(dummy)
    lines(1:pp,mean.Geo.KG[k,seq(1,pp)], type = "p", col=cols[k], lwd=1.5, pch=plotchar[k], cex=2.5)
  }
  # plot total median
  tot.split <- split(mdata,fKG)
  tot.med <- lapply(tot.split, median, simplify = TRUE)
  lines(1:pp,tot.med[seq(1,pp)], type = "p", col=cols[g+1], lwd=1.5, pch=plotchar[k],cex=2.5)
  
  # add a legend
  if (i == length(vars[pinds])){
    legend(x = "topleft", y= 0.9 , lab, cex=0.8, col=cols, pch=plotchar)
  }
}
dev.copy(which = 4)
dev.off(which = 4)
rm(mk,mkKG,mKG,msplit,mdata,mean.Geo.KG,tot.split,tot.med)

############################################################################
# Plot only medians for global data-set
dev.new()
setEPS()
postscript(paste("Global_medians.eps", sep  = "", collapse = NULL))
dev.set(which = 2)
par(oma=c(0,0,0,0),mfrow=c(2,3),cex=0.8)
for (i in 1:length(vars[pinds])){
  mdata <- dat[,which(pinds)[i]]
  Geo.split <- split(mdata,GeoIDf)
  Geo.med <- lapply(Geo.split, median, simplify = TRUE)
  Geo.med <- simplify2array(Geo.med)
  plot(c(seq(1,g)),Geo.med, pch=plotchar, col = cols,ylab = vars[pinds][i], xlab="",
        cex = 2.5)
}
dev.copy(which = 4)
dev.off(which = 4)
rm(Geo.med,mdata)
# LANDCOVER DATA #########################################################
#########################################################################
# histograms, PLEASE GO BACK TO THE BEGINNING AND RERUN GEOIDF TO HAVE THE NATURAL ORDER

LC.split <- split(as.numeric(lcID),GeoIDf)
# par(oma=c(0,0,0,0),mar = c(3,3,3,2),mfrow=c(3,3),cex=0.8)
# lapply(LC.split, function(i) hist(i,breaks = seq(0.5,6.5),labels = LClab, xlab="",
#                                   ylab="", xaxt = "n"))
# calculate percentages
LC.counts <- lapply(LC.split, function (i) table(i))  # counts LC per Geo
LC.counts <- simplify2array(LC.counts)                # simplifiy to array
LC.counts <- rbind(LC.counts,apply(LC.counts,2,sum))  # calculate total area
colnames(LC.counts) <- lab
rownames(LC.counts) <- c(LClab,"total")
perc.LC.Geo <- LC.counts                       # copy to keep the names
for (i in 1:g){                                # calculate percentages for every Geo
  perc.LC.Geo[1:(length(LClab)+1),i] <- perc.LC.Geo[1:(length(LClab)+1),i]/perc.LC.Geo[
    (length(LClab))+1,i]*1e2
}
tot.count <- table(lcID)                       # total count of different LC's
tot.perc <- tot.count / sum(tot.count) *100    # convert to %
perc.LC.Geo <- cbind(perc.LC.Geo,c(tot.perc,100))
colnames(perc.LC.Geo)[8] <- "Total LC"

write.table(perc.LC.Geo, "LC_stats_200m.csv", append = FALSE, sep = ",", dec = ".",
            row.names = T, col.names = T)

# plot LC per Geo #####################################################
library(tidyverse)
cols = c("#FFEBAF","#38A800","#727272","#FFAA00","#BED2FF","#0070FF","#734C00","#000000")
plotchar <- c(seq(16,18),15,15,19,17,15)
ps <- 1.5
dev.new()
setEPS()
postscript(paste("LC_per_Geo.eps", sep  = "", collapse = NULL))
dev.set(which = 2)
perc.LC.Geo <- as.data.frame(perc.LC.Geo)
Glist <- colnames(perc.LC.Geo)
perc.LC.Geo <- perc.LC.Geo[-c(7),]        # remove total row since its only 100 anyway...
pp <- ggplot(perc.LC.Geo, aes(y = perc.LC.Geo$su, x = seq(1,6))) +
  geom_point(aes(color = cols[1], shape = plotchar[1], size = ps)) + 
  geom_point(aes(y = perc.LC.Geo$vc, color = cols[2], shape = plotchar[2], size = ps)) + 
  geom_point(aes(y = perc.LC.Geo$pl, color = cols[3], shape = plotchar[3], size = ps)) +
  geom_point(aes(y = perc.LC.Geo$ss, color = cols[4], shape = plotchar[4], size = ps)) +
  geom_point(aes(y = perc.LC.Geo$sm, color = cols[5], shape = plotchar[5], size = ps)) +
  geom_point(aes(y = perc.LC.Geo$sc, color = cols[6], shape = plotchar[6], size = ps)) +
  geom_point(aes(y = perc.LC.Geo$mt, color = cols[7], shape = plotchar[7], size = ps)) +
  geom_point(aes(y = perc.LC.Geo$`Total LC`, color = cols[8], shape = plotchar[8], size = ps)) +
  scale_shape_identity() + 
  scale_color_identity() + 
  labs(x = "LC type", y = "% frequency") + 
  scale_x_continuous(breaks = seq(1,6),labels = LClab) + 
  theme_bw()
pp
dev.copy(which = 4)
dev.off(which = 4)


#############################################################################################
#############################################################################################
# REGRESSION AND MORE STATISTICS #######################################################


# Do ANOVA #####################################################
# do ANOVA for all variables manually check output
library(mvtnorm)
library(zoo)
library(multcomp)
anova = list()
levTest = list()
anova.sum = list()
lev.sum = list()
glht.sum = list()
shap = list()
for (i in 1:length(vars[inds])){
  anova[[i]] = list()
  anova.sum[[i]] = list()
  levTest[[i]] = list()
  lev.sum[[i]] = list()
  glht.sum[[i]] =list()
  shap[[i]] = list()
  # Do ANOVA, ghlt, Levene test
  anova[[i]] <- aov(dat[,vars[inds][i]] ~ GeoIDf)
  
  # Levene Test checks the AOV assumption of homogenoeus variance, p should be > 0.05! 
  # to show that the variances are not significantly different between groups
  levTest[[i]] <- leveneTest(dat[,vars[inds][i]] ~ GeoIDf) 
  
  # Turkey HSD test shows between which groups the differences are significant, using a linear model
  glht.sum[[i]] <- summary(glht(anova[[i]], linfct = mcp(GeoIDf= "Tukey")))
  
  # Extract the residuals
  aov_residuals <- residuals(object = anova[[i]] )
  # Run Shapiro-Wilk test to check AOV normality assumption, p should be > 0.05 for normality to be probably not violated
  # shap[[i]] = shapiro.test(x = aov_residuals ) # only possible for < 5000 samples
  
  # save results
  anova.sum[[i]] <- summary(anova[[i]])
  lev.sum[[i]] <- summary(levTest[[i]])
}

# KOLMOGOROW-SMIRNOV-TEST & 2-SAMPLE KUIPER TEST ############################################
library("remotes")
library("nonpar")
library("kuiper.2samp")
ks.stats <- list()
cc.stats <- list()
kui.stats <- list()
for (i in 1:length(vars[inds])){  # Loop through variables
  Test.data <- dat[,which(inds)[i]]   # data for this variable
  data.split <- split(Test.data,GeoIDf)   # split data by Geo
  ks.stats[[i]] <- matrix(, nrow = g, ncol = length(data.split)-1)
  cc.stats[[i]] <- matrix(, nrow = g, ncol = length(data.split)-1)
  kui.stats[[i]] <- matrix(, nrow = g, ncol = length(data.split)-1)
  for (j in 1:g){                  # Loop through different rock types
    sample1 <- data.split[[1]]
    for (h in 2:length(data.split)){
      sample2 <- data.split[[h]]
      
      ks <- ks.test(sample1,sample2)
      ks.stats[[i]][j,h-1] <- ks[[2]]
      
      # cc <- cucconi.test(sample1,sample2)
      # cc.stats[[i]][j,h-1] <- cc[[2]]
      
      # kui <- kuiper.2samp(sample1,sample2)
      # kui.stats[[i]][j,h-1] <- kui[[2]]
    }
  }    
}


#############################################################################################
# MULTINOMINAL REGRESSION ###################################################################
#############################################################################################
GeoIDf<- relevel(GeoIDf, ref = 4)  # make ss reference level

z = 1
num <- 5 #??? number of predictor variables
if (z == 1){
  # regression in z-score
  Z <- scale(dat[inds], center = TRUE, scale = TRUE)
  mod <- multinom(GeoIDf ~ dat$Locrel + dat$ksn + dat$NDVI + dat$tet_corr + dat$Amph_corr)
  mod_data <-summary(mod)
  }else{
  # normal regression
  mod <- multinom(GeoIDf ~ dat$Locrel + dat$ksn + dat$NDVI + dat$tet_corr + dat$Amph_corr)
  mod_data <- summary(mod)
  # exp(coef(mod))  # exp because returns log coefficients
  psr = PseudoR2(mod,c("all"))
}


# PLOT REGRESSION RESULTS ##################################################################
dev.new()
setEPS()
postscript("MNR_coefficientsPV5.eps")
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(2,3),cex=0.8,lwd=0.8) 
for (i in 1:(g-1)){
  plot(mod_data$coefficients[i,2:ncol(mod_data$coefficients)],1:num, xlab = "coefficient",
       ylab = "variable" , pch = 19,yaxt = "n",  main = lab[as.numeric(row.names(mod_data$coefficients)[i])])
  arrows(mod_data$coefficients[i,2:ncol(mod_data$coefficients)]-mod_data$standard.errors[i,2:ncol(mod_data$standard.errors)]
         , 1:num,mod_data$coefficients[i,2:ncol(mod_data$coefficients)]+mod_data$standard.errors[i,2:ncol(mod_data$standard.errors)]
         , length=0.05, angle=180, code=3)
  abline(v = 0, col="black")
  axis(2, at=1:5, labels= c("Locrel", "ksn", "NDVI", "Tetrapods", "Amphibians"))
}
dev.copy(which = 4)
dev.off(which = 4)


# PV = 4
# MULTINOMINAL REGRESSION ###################################################################
z = 1
num <- 4 #??? number of predictor variables
if (z == 1){
  # regression in z-score
  z1 <- scale(dat$Locrel, center = TRUE, scale = TRUE)
  z2 <- scale(dat$NDVI, center = TRUE, scale = TRUE)
  z3 <- scale(dat$tet_corr, center = TRUE, scale = TRUE)
  z4 <- scale(dat$Amph_corr, center = TRUE, scale = TRUE)
  
  mod <- multinom(GeoIDf ~ z1 + z2+ z3 +z4)
  mod_data <-summary(mod)
}else{
  # normal regression
  mod <- multinom(GeoIDf ~ dat$Locrel + dat$NDVI + dat$tet_corr + dat$Amph_corr)
  mod_data <- summary(mod)
  # exp(coef(mod))  # exp because returns log coefficients
  psr = PseudoR2(mod,c("all"))
}

# vertical coefficient plot 
dev.new()
setEPS()
postscript("MNR_coefficientsPV4.eps")
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(2,3),cex=0.8,lwd=0.8) 
for (i in 1:(g-1)){
  plot(mod_data$coefficients[i,2:ncol(mod_data$coefficients)],1:num, xlab = "coefficient",
       ylab = "variable" , pch = 19,yaxt = "n",  main = lab[as.numeric(row.names(mod_data$coefficients)[i])])
  arrows(mod_data$coefficients[i,2:ncol(mod_data$coefficients)]-mod_data$standard.errors[i,2:ncol(mod_data$standard.errors)]
         , 1:num,mod_data$coefficients[i,2:ncol(mod_data$coefficients)]+mod_data$standard.errors[i,2:ncol(mod_data$standard.errors)]
         , length=0.05, angle=180, code=3)
  abline(v = 0, col="black")
  axis(2, at=1:4, labels= c("locrel", "NDVI", "Tetrapods", "Amphibians"))
}
dev.copy(which = 4)
dev.off(which = 4)


# PV = 3
# MULTINOMINAL REGRESSION ###################################################################
z = 1
num <- 3 #??? number of predictor variables
if (z == 1){
  # regression in z-score
  Z <- scale(dat[inds], center = TRUE, scale = TRUE)
  mod <- multinom(GeoIDf ~ dat$Locrel + dat$NDVI + dat$Amph_corr)
  mod_data <-summary(mod)
}else{
  # normal regression
  mod <- multinom(GeoIDf ~ dat$Locrel + dat$NDVI + dat$Amph_corr)
  mod_data <- summary(mod)
  # exp(coef(mod))  # exp because returns log coefficients
  psr = PseudoR2(mod,c("all"))
}

# vertical coefficient plot NORMAL SCORE
dev.new()
setEPS()
postscript("MNR_coefficientsPV3.eps")
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(2,3),cex=0.8,lwd=0.8) 
for (i in 1:(g-1)){
  plot(mod_data$coefficients[i,2:ncol(mod_data$coefficients)],1:num, xlab = "coefficient",
       ylab = "variable" , pch = 19,yaxt = "n",  main = lab[as.numeric(row.names(mod_data$coefficients)[i])])
  arrows(mod_data$coefficients[i,2:ncol(mod_data$coefficients)]-mod_data$standard.errors[i,2:ncol(mod_data$standard.errors)]
         , 1:num,mod_data$coefficients[i,2:ncol(mod_data$coefficients)]+mod_data$standard.errors[i,2:ncol(mod_data$standard.errors)]
         , length=0.05, angle=180, code=3)
  abline(v = 0, col="black")
  axis(2, at=1:3, labels= c("locrel", "NDVI", "Amphibians"))
}
dev.copy(which = 4)
dev.off(which = 4)
