# My very first R Script!!!!!!
# Richard Ott, 2019
rm(list = ls())
mthres <- 50
# setwd(paste('C:/Richard/PhD_ETH/matlab/Multivariate Analysis/', mthres, 'm/for_stat_consulting',sep=''))
setwd(paste('E:/Richard/Global_Hypsometry/matlab/global correlations/', mthres, 'm_relief/for_stat_consulting',sep=''))
library(DescTools)
library(nnet)
library(graphics)
library(car)
library(ggplot2)
library(plyr)
library(nlme)


# LOAD DATA ############################################################################
vars <- c("Elev","Locrel","slope","P",'T','NDVI','Lat','ksn','tetrapods'
          ,'Amphibians','Mammals','KG','Geo','LC','SR')
dat <- read.table(paste('global_data_',mthres,'mLC2.txt',sep=''), header = FALSE, sep = ",", col.names	= vars)
dat$Amphibians[dat$Amphibians == -9999 | dat$Amphibians == -32768] = 0   # Amph have lots of nan values that should be set to zero for analysis
corr = 1 # Do you want biological variables to be corrected for T and P influence?
setwd(paste('E:/Richard/Global_Hypsometry/matlab/global correlations/', mthres, 'm_relief/pixel',sep=''))


# GEO FACTORS ########################################################################
depVar <- which(vars == 'Geo')                     # dependent variable, should be last variable
#         E Loc Sl P  T  ND L ksn te am   ma  KG Geo LC   SR  NDVIc Amphc Tetc
inds <- c(T, T, T, T, F, F, F, T, F  ,F ,  F,  F, F , F ,  T,   T,    T,    T)

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
  rm(lm_results,NDVIemp,NDVIcorr)
  
  # AMPHiIBIANS #######################
  linMod <- lm(dat$Amphibians ~ dat$T + dat$P, data= dat)
  lm_results <- summary(linMod)
  Pthres = dat$P
  Pthres[Pthres > thrs] = thrs
  Amph.emp <- lm_results$coefficients[[1,1]] + dat$T * lm_results$coefficients[[2,1]] + Pthres * lm_results$coefficients[[3,1]]
  Amph.corr <- dat$Amphibians - Amph.emp
  dat$Amph_corr <- Amph.corr
  rm(lm_results,Amph.corr,Amph.emp)
  
  # TETRAPODS #######################
  linMod <- lm(dat$tetrapods ~ dat$T + dat$P, data= dat)
  lm_results <- summary(linMod)
  Pthres = dat$P
  Pthres[Pthres > thrs] = thrs
  Tet.emp <- lm_results$coefficients[[1,1]] + dat$T * lm_results$coefficients[[2,1]] + Pthres * lm_results$coefficients[[3,1]]
  Tet.corr <- dat$tetrapods - Tet.emp
  dat$tet_corr <- Tet.corr
  
  rm(linMod,lm_results,Tet.emp,Tet.corr,Pthres)
  vars[16:18] <- c('NDVIc','Amphc','Tetc')
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
postscript("boxplot_all_corrected.eps", horizontal = F, pointsize = 7)
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(3,3),cex=0.5,lwd=0.35)
llim <- c(0,    0,  0, 300, 0,  -150, 0,  -100, -15,  -100, -15)     # corrected limits
ulim <- c(2000,400,120,1500,100, 500, 250, 500, 3,   150, 10)
cols = c("#FFEBAF","#38A800","#727272","#FFAA00","#BED2FF","#0070FF","#734C00","#000000")
# dev.set(which == 2)
for (i in 1:length(vars[inds])){ # loop through variables
  data.by.geo <- split(dat[,which(inds)[i]], GeoIDf)     # split by geology
  
  # boxplot(data.by.geo, xaxt = "n", main = vars[inds][i], outline = F, boxlty = 0, whisklty = 0, staplelty = 0, ylim = c(llim[i],ulim[i]))
  # boxplot(data.by.geo, xaxt = "n", main = vars[inds][i])
  
  boxplot(data.by.geo, xaxt = "n", main = vars[inds][i], outline = F, staplelty = 0, whisklty = 1, col = cols)
  axis(1, at=1:length(lab), labels=lab)
}
dev.copy(which = 4)
dev.off(which = 4)
dev.off()


# TOPO-BOXPLOTS TECTONICALLY INACTIVE VS ACTIVE #######################################
#######################################################################################
Tect <- dat$SR == 0              # inactive = True, active = False
Geo.by.tect <- split(GeoIDf,Tect)
#              E Loc Sl P  T  ND L ksn te am   ma  KG Geo LC   SR  NDVIc Amphc Tetc
topo.inds <- c(T, T, T, F, F, F, F, T, F  ,F ,  F,  F, F , F ,  F,   F,    F,    F)

dev.new()
setEPS()
postscript("boxplot_tectonic.eps", horizontal = F, pointsize = 7)
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(2,4),cex=0.35,lwd=0.35)
cols = c("#FFEBAF","#38A800","#727272","#FFAA00","#BED2FF","#0070FF","#734C00","#000000")
# dev.set(which == 2)

for (i in 1:length(which(topo.inds))){ # loop through variables
  # rearrange data matrix
  data.by.tect <- split(dat[,which(topo.inds)[i]], Tect)     # split by tectonics
  
  # first do active areas
  data.by.geo <- split(data.by.tect[[1]],Geo.by.tect[[1]])        # split by GEO
  boxplot(data.by.geo,names = lab ,main = paste(vars[topo.inds][i], '_act', sep='')
          , outline = F, boxlty = 1, whisklty = 1, staplelty = 0, col = cols, cex.axis = 2)
  
  # now do inactive areas
  rm(data.by.geo)
  data.by.geo <- split(data.by.tect[[2]],Geo.by.tect[[2]])        # split by GEO
  lims <- par("usr")[3:4]
  boxplot(data.by.geo,names = lab ,main = paste(vars[topo.inds][i], '_in', sep='')
          , outline = F, boxlty = 1, whisklty = 1, staplelty = 0, col = cols,cex.axis = 2, ylims = lims)
}
dev.copy(which = 4)
dev.off(which = 4)
dev.off()
rm(Geo.by.tect,topo.inds,data.by.geo,lims,data.by.tect)



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


#############################################################################################
## P-T-RASTER FOR BIO-VARIABLES #############################################################
#############################################################################################
library(Hmisc)
library(mapplots)
library(fields)
library(spatstat)
#             E Loc Sl P  T  ND L ksn te am   ma  KG Geo LC   SR  NDVIc Amphc Tetc
bio.inds <- c(F, F, F, F, F, T, F, F, T  ,T ,  F,  F, F , F ,  F,   F,    F,    F)
T.quart <- quantile(dat$T)   # determine quartile boundaries
P.quart <- quantile(dat$P)

T.intervals <- findInterval(dat$T,T.quart)             # find P & T intervals for every value
T.intervals[which(T.intervals == 5)] <- 4              # intervals macht größer-gleich, deshalb ist höchster Wert in neuer Kategorie und muss manuell angepasst werden
P.intervals <- findInterval(dat$P,P.quart)
P.intervals[which(P.intervals == 5)] <- 4

TP.val <- rep(NA,length(GeoIDf))
n <- 0
for (i in 1:4){                                # put all values in 16 grid raster defined by quartiles
  P.inds <- P.intervals == i                   # T=0 & P=0 == 1; then rows increase in T drection, new column with higher precip
  TP.val[P.inds] <- n + T.intervals[P.inds]
  n <- n + 4
}
TP.val <- factor(TP.val)
TP.frequencies <- list()

# IMPORTANT: I only plot data were there are more samples than a certain threshold
plot.thres <- 1e3                        # below this number of samples a certain T-P- region is not viewed as sampled enough

# Calculate the means for all bio variables in the TP-grids
all.means <- list()
diff.means <- list()
zlims <- list('1' = c(-25,25), '2'= c(-80,80), '3' = c(-12,12))
cols <- colorRampPalette(c("firebrick","yellow","green4"))(100)
dev.new()
for (i in 1:length(which(bio.inds))){     # loop through all bio variables
  setEPS()
  postscript(paste('TP_', vars[which(bio.inds)[i]] ,'.eps', sep=''), horizontal = F, pointsize = 7,
             width = 185/25.4, height = 100/25.4, paper = "special", onefile = FALSE)
  dev.set(which = 2)
  
  all.means[[i]] <- list()
  diff.means[[i]] <- list()
  # Calculate total mean
  var <- dat[,which(bio.inds)[i]]
  split.var <- split(var,TP.val)
  TP.frequencies <- sapply(split.var, function(y) length(y))
  split.means <- sapply(1:length(split.var), function(i) mean(split.var[[i]]))     # can change this to median
  # split.means <- sapply(1:length(split.weight), function(i) median(split.var[[i]]))     # can change this to median
  all.means[[i]][[g+1]] <- matrix(split.means,4,4, byrow = T)
  TP.frequencies <- matrix(TP.frequencies,4,4, byrow = T)
  
  par(oma=c(0,0,2,0),mfrow=c(2,4),cex=0.5)
  for (j in 1:g){                                  # loop through GEO's
    var <- dat[GeoIDf == j, which(bio.inds)[i]]    # get values of this lith and variable
    split.var <- split(var, TP.val[GeoIDf == j])   # split by TP-index
    id <- sapply(split.var, function(i) length(i))
    id <- id < plot.thres
    split.means <- sapply(1:length(split.var), function(i) mean(split.var[[i]]))     # can change this to median
    # split.means <- sapply(1:length(split.weight), function(i) median(split.var[[i]]))     # can change this to median
    split.means[id] <- NA                              # This removes regions that dont have enough data from the plot
    split.means <- matrix(split.means,4,4, byrow = T)
    
    all.means[[i]][[j]] <- split.means
    diff.means[[i]][[j]] <- split.means - all.means[[i]][[g+1]]
    # values below or above the z-limit should be shown by th max color. By default R will leave those values out so I will correct them for the image process
    diff.means[[i]][[j]][which(diff.means[[i]][[j]] > zlims[[i]][2])] <- zlims[[i]][2] # replace above maximum zlim values with max zlim
    diff.means[[i]][[j]][which(diff.means[[i]][[j]] < zlims[[i]][1])] <- zlims[[i]][1]
    
    image(x =seq(1,4), y =seq(1,4), z =diff.means[[i]][[j]], main = lab[j],xaxt = 'n', yaxt = 'n', zlim = zlims[[i]],
          xlab = 'T', ylab = 'P', col = cols)
    axis(1,at=seq(0.5,4.5,1), labels= as.character(T.quart/10))
    axis(2,at=seq(0.5,4.5,1), labels= as.character(P.quart))
  }
  # draw frequencies of all data
  image(x =seq(1,4), y =seq(1,4), z = TP.frequencies, main = 'sample frequencies',xaxt = 'n', yaxt = 'n', xlab = 'T', ylab = 'P',
        col = cols)
  axis(1,at=seq(0.5,4.5,1), labels= as.character(T.quart/10))
  axis(2,at=seq(0.5,4.5,1), labels= as.character(P.quart))

  dev.copy(which = 4)
  dev.off(which = 4)
}
rm(split.var,diff.means,all.means,zlims,var,T.intervals,P.intervals)
dev.off()

# BOOTSTRAP SIGNIFICANCE BETWEEN T-P-INTERVALS ###########################################
library(multcomp)
library(openxlsx)
library(broom)
m <- 1                           # Method for bootstrap sampling: 0 - same from every T-P-bin, 1- identical to original T-P-bin distribution
b <- 5e3                         # number of bootstrap samples you want to take
# SIZE OF BOOTSTRAP SAMPLE #
CI.sig <- 95                       # confidence interval difference level
# Based on POWER ANALYSIS
library(pwr)
sig <- 0.05                      # significance level, Type I error
alph <- 0.95                     # Power, (1-Type II error), chance of finding an effect
effect.size <- 0.1               # levels by Cohen 1992, but can vary fro field to field, 0.1 -small, 0.3 medium, 0.5 large
samp.size <- pwr.anova.test(k = g, f = effect.size, sig.level = sig, power = alph)
bs.size <- ceiling(samp.size$n)


# Method 1 ####################################################################
mean.median.diff <- lapply(1:length(which(bio.inds)), function(i) matrix(NA, nrow = g, ncol = g))      # mean difference in the medians between 2 groups
low.CI.diff <- lapply(1:length(which(bio.inds)), function(i) matrix(NA, nrow = g, ncol = g))           # lower bound oconfidence interval between 2 groups
up.CI.diff <- lapply(1:length(which(bio.inds)), function(i) matrix(NA, nrow = g, ncol = g))            # upper bound of CI 
sig.matrix <- lapply(1:length(which(bio.inds)), function(i) matrix(NA, nrow = g, ncol = g))            # upper bound of CI 

dev.new()
setEPS()
postscript(paste('Bootstrap_M1_TP_ES' ,effect.size, '.eps',sep=''), horizontal = F, pointsize = 7)
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(2,2),cex=1)
for (i in 1:length(which(bio.inds))){
  TP.freq <- TP.frequencies/sum(TP.frequencies)
  bs.data <- lapply(1:g, function(i) vector('numeric',b))
  for (j in 1:g){
    vals <- dat[which(GeoIDf == j), which(bio.inds)[i]]          # select the data froom the desired variable and lithology 
    TP.bin <- TP.val[which(GeoIDf == j)]
    vals.TP <- split(vals,TP.bin)
    if(m == 0){
      samp.per.bin <- ceiling(bs.size/16)   # divide by number of bins to know how many samples from every bin should be drawn
    }else{
      samp.per.bin <- round(TP.freq*bs.size)   # how many samples per bin
      samp.per.bin <- as.vector(t(samp.per.bin))         # transform to vector
    }
    for (k in 1:b){
      if (m == 0){
        bs.sample <- unlist(sapply(1:16,function(z) sample(vals.TP[[z]],size= samp.per.bin, replace = T)))         # bootstrap sample, this is slow...
      }else{
        bs.sample <- unlist(sapply(1:16,function(z) sample(vals.TP[[z]],size= samp.per.bin[z], replace = T)))      # bootstrap sample
      }
      bs.data[[j]][k] <- median(bs.sample)   # THE STATISTIC!!!
    }
  }
  for (j in 1:(g-1)){    # loop through columns of the final matrices
    for (k in (j+1):g){  #  loop through rows of lower triangle
      get.data.1 <- bs.data[[j]]    # data from column
      get.data.2 <- bs.data[[k]]    # data from row
      delta <- get.data.2 - get.data.1
      mean.median.diff[[i]][k,j] <- mean(delta)
      low.CI.diff[[i]][k,j] <- quantile(delta,probs = (1e2-CI.sig)/1e2)
      up.CI.diff[[i]][k,j] <- quantile(delta,probs = (CI.sig)/1e2)
      sig.matrix[[i]][k,j] <- findInterval(0,c(low.CI.diff[[i]][k,j],up.CI.diff[[i]][k,j]))   
      # This means 0 - row is higher than column, 1 - insignificant, 2- column higher than row
    }
  }
  # plot the signifcance matrices ############################################
  rotate <- function(x) t(apply(x, 2, rev))
  
  image(x =seq(1,g), y =seq(1,g), z =rotate(sig.matrix[[i]]), main = vars[which(bio.inds)[i]],xaxt = 'n', yaxt = 'n', xlab = '', ylab = ''
        , col = colorRampPalette(c("blue", "white", "red"))(3), zlim = c(0,2))
  axis(1,at=seq(1,g,1), labels= lab[1:g])
  axis(2,at=seq(1,g,1), labels= rev(lab[1:g]))
  grid(nx=g, ny=g, lty=1)
}
write.xlsx(sig.matrix, file= paste("TP_LogPmat_M1_ES", effect.size, ".xlsx",sep=''),asTable = F, sheetName= vars[bio.inds], append=TRUE, row.names=T, col.names = T)
dev.copy(which = 4)
dev.off(which = 4)
dev.off()

rm(TP.val,get.data.1,get.data.2,delta,mean.median.diff,low.CI.diff,up.CI.diff,samp.per.bin,bs.data,bs.sample,vals,m,b,bs.size,effect.size,alph,CI.sig)

# METHOD 2 ################################################################
wilcox.p.perc <- lapply(1:length(which(bio.inds)), function (i) i)
dev.new()
setEPS()
postscript(paste('Bootstrap_M1_TP_ES', effect.size, '.eps'), horizontal = F, pointsize = 7)
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(3,3),cex=0.5)
for (i in 1:length(which(bio.inds))){
  split.by.geo <- split(dat[,which(bio.inds)[i]], GeoIDf)   # make groups
  TP.per.geo <- split(TP.vals, GeoIDf)
  TP.geo.split <- lapply(1:g, function (l) split(split.by.geo[[l]],TP.per.geo[[l]]))
  wilcox.p <- lapply(1:b, function(l) matrix(NA, nrow = g-1, ncol = g-1))
  for (j in 1:b){
    bs.sample <- list()
    bs.sample <- lapply(1:g, function (z) sample(split.by.geo[[z]], size= bs.size, replace = T))   # take a bootstrap sample for every GEO
    bs.sample <- unlist(bs.sample)
    bs.group <- rep(seq(1,g), each = bs.size)
    wilcox <- pairwise.wilcox.test(bs.sample, bs.group, p.adjust.method = 'BH')  # Do pairwise wilcox test to test for differences in the distributions
    wilcox.p[[j]] <- wilcox$p.value                                              # save only p-values
  }
  wilcox.binary <- lapply(wilcox.p, function(z) z < sig)                         # make binary matrices with significant or not
  wilcox.sum <- Reduce('+',wilcox.binary)                                        # sum them all up
  wilcox.p.perc[[i]] <- wilcox.sum/b                                             # Calculate percentage of significant outcome
  
  # plot the signifcance matrices ############################################
  rotate <- function(x) t(apply(x, 2, rev))
  
  image(x =seq(1,6), y =seq(1,6), z =rotate(wilcox.p.perc[[i]]), main = vars[which(bio.inds)[i]],xaxt = 'n', yaxt = 'n', xlab = '', ylab = '',
        col = colorRampPalette(c("blue","white","red"))(100), zlim = c(0,1))
  axis(1,at=seq(1,6,1), labels= lab[1:(g-1)])
  axis(2,at=seq(1,6,1), labels= rev(lab[2:g]))
  grid(nx=(g-1), ny=(g-1), lty=1)
}
dev.copy(which = 4)
dev.off(which = 4)
write.xlsx(wilcox.p.perc, file= paste("Pmat_M2_TP_ES", effect.size, ".xlsx",sep=''),asTable = F, sheetName= vars[bio.inds], append=TRUE, row.names=T, col.names = T)

rm(split.by.geo,wilcox.p,wilcox.binary,wilcox.p.perc,wilcox,bio.inds)
dev.off()


##########################################################################
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

# BOOTSTRAP STATISTICs #####################################################################
###########################################################################################
# FOR ALL DATA ###########################################################################
library(multcomp)
library(openxlsx)
library(broom)
#              E Loc Sl P  T  ND L ksn te am   ma  KG Geo LC   SR  NDVIc Amphc Tetc
boot.inds <- c(T, T, T, T, F, F, F, T, F  ,F ,  F,  F, F , F ,  F,   T,    T,    T)
b <- 5e3                         # number of bootstrap samples you want to take

# SIZE OF BOOTSTRAP SAMPLE #
# METHOD 1: An arbirtary small subset of the total data
bs.size <- nrow(dat)/g/1e3         # size of bootstrap sample based on % of pixel number
CI.sig <- 95                       # confidence interval difference level
# METHOD 2: Based on POWER ANALYSIS
library(pwr)
sig <- 0.05                      # significance level, Type I error
alph <- 0.95                     # Power, (1-Type II error), chance of finding an effect
effect.size <- 0.1               # levels by Cohen 1992, but can vary fro field to field, 0.1 -small, 0.3 medium, 0.5 large
samp.size <- pwr.anova.test(k = g, f = effect.size, sig.level = sig, power = alph)
bs.size <- ceiling(samp.size$n)

# METHOD 1, only keep statistic ########################################################
mean.median.diff <- lapply(1:length(which(boot.inds)), function(i) matrix(NA, nrow = g, ncol = g))      # mean difference in the medians between 2 groups
low.CI.diff <- lapply(1:length(which(boot.inds)), function(i) matrix(NA, nrow = g, ncol = g))           # lower bound oconfidence interval between 2 groups
up.CI.diff <- lapply(1:length(which(boot.inds)), function(i) matrix(NA, nrow = g, ncol = g))            # upper bound of CI 
sig.matrix <- lapply(1:length(which(boot.inds)), function(i) matrix(NA, nrow = g, ncol = g))            # upper bound of CI 

dev.new()
setEPS()
postscript(paste('Bootstrap_M1_all_ES' ,effect.size, '.eps',sep=''), horizontal = F, pointsize = 7)
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(3,3),cex=1)
for (i in 1:length(which(boot.inds))){
  split.by.geo <- split(dat[,which(boot.inds)[i]], GeoIDf)   # make groups
  bs.data <- lapply(1:g, function(l) vector('numeric',b))
  for (j in 1:g){                                                                # loop through list of split.by.geo
    for (k in 1:b){
      bs.sample <- sample(split.by.geo[[j]], size= bs.size, replace = T)         # bootstrap sample
      bs.data[[j]][k] <- median(bs.sample)                                       # THE STATISTIC YOU WISH FOR
    }
  }

  for (j in 1:(g-1)){    # loop through columns of the final matrices
    for (k in (j+1):g){  #  loop through rows of lower triangle
      get.data.1 <- bs.data[[j]]    # data from column
      get.data.2 <- bs.data[[k]]    # data from row
      delta <- get.data.2 - get.data.1
      mean.median.diff[[i]][k,j] <- mean(delta)
      low.CI.diff[[i]][k,j] <- quantile(delta,probs = (1e2-CI.sig)/1e2)
      up.CI.diff[[i]][k,j] <- quantile(delta,probs = (CI.sig)/1e2)
      sig.matrix[[i]][k,j] <- findInterval(0,c(low.CI.diff[[i]][k,j],up.CI.diff[[i]][k,j]))   
      # This means 0 - row is higher than column, 1 - insignificant, 2- column higher than row
    }
  }
  
  # plot the signifcance matrices ############################################
  rotate <- function(x) t(apply(x, 2, rev))
  
  image(x =seq(1,g), y =seq(1,g), z =rotate(sig.matrix[[i]]), main = vars[which(boot.inds)[i]],xaxt = 'n', yaxt = 'n', xlab = '', ylab = '',
        col = colorRampPalette(c("blue", "white", "red"))(3), zlim = c(0,2))
  axis(1,at=seq(1,g,1), labels= lab[1:g])
  axis(2,at=seq(1,g,1), labels= rev(lab[1:g]))
  grid(nx=g, ny=g, lty=1)
}
dev.copy(which = 4)
dev.off(which = 4)
write.xlsx(sig.matrix, file= paste("LogPmat_M1_ES", effect.size, ".xlsx",sep=''),asTable = F, sheetName= vars[boot.inds], append=TRUE, row.names=T, col.names = T)
dev.off()

rm(mean.median.diff,low.CI.diff,up.CI.diff,sig.matrix,split.by.geo,bs.data,get.data.1,get.data.2,delta,sig.plot)


# METHOD 2, calculate statistic directly on whole distribution #########################
wilcox.p.perc <- lapply(1:length(which(boot.inds)), function (i) i)
dev.new()
setEPS()
postscript(paste('Bootstrap_M2_all_ES', effect.size, '.eps'), horizontal = F, pointsize = 7)
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(3,3),cex=0.5)
for (i in 1:length(which(boot.inds))){
  split.by.geo <- split(dat[,which(boot.inds)[i]], GeoIDf)   # make groups
  wilcox.p <- lapply(1:b, function(l) matrix(NA, nrow = g-1, ncol = g-1))
  for (j in 1:b){
    bs.sample <- list()
    bs.sample <- lapply(1:g, function (z) sample(split.by.geo[[z]], size= bs.size, replace = T))   # take a bootstrap sample for every GEO
    bs.sample <- unlist(bs.sample)
    bs.group <- rep(seq(1,g), each = bs.size)
    wilcox <- pairwise.wilcox.test(bs.sample, bs.group, p.adjust.method = 'BH')  # Do pairwise wilcox test to test for differences in the distributions
    wilcox.p[[j]] <- wilcox$p.value                                              # save only p-values
  }
  wilcox.binary <- lapply(wilcox.p, function(z) z < sig)                         # make binary matrices with significant or not
  wilcox.sum <- Reduce('+',wilcox.binary)                                        # sum them all up
  wilcox.p.perc[[i]] <- wilcox.sum/b                                             # Calculate percentage of significant outcome
  
  # plot the signifcance matrices ############################################
  rotate <- function(x) t(apply(x, 2, rev))
  
  image(x =seq(1,6), y =seq(1,6), z =rotate(wilcox.p.perc[[i]]), main = vars[which(boot.inds)[i]],xaxt = 'n', yaxt = 'n', xlab = '', ylab = '',
        col = colorRampPalette(c("blue","white","red"))(100), zlim = c(0,1))
  axis(1,at=seq(1,6,1), labels= lab[1:(g-1)])
  axis(2,at=seq(1,6,1), labels= rev(lab[2:g]))
  grid(nx=(g-1), ny=(g-1), lty=1)
}
dev.copy(which = 4)
dev.off(which = 4)
write.xlsx(wilcox.p.perc, file= paste("Pmat_M2_ES", effect.size, ".xlsx",sep=''),asTable = F, sheetName= vars[boot.inds], append=TRUE, row.names=T, col.names = T)

rm(split.by.geo,wilcox.p,wilcox.binary,wilcox.p.perc,wilcox)
dev.off()


# BOOTSTRAP BY TECTONIC ZONES ################################################
##############################################################################
library(multcomp)
library(openxlsx)
library(broom)
#              E Loc Sl P  T  ND L ksn te am   ma  KG Geo LC   SR  NDVIc Amphc Tetc
topo.inds <- c(T, T, T, F, F, F, F, T, F  ,F ,  F,  F, F , F ,  F,   F,    F,    F)
b <- 5e3                         # number of bootstrap samples you want to take

# SIZE OF BOOTSTRAP SAMPLE #
# METHOD 1: An arbirtary small subset of the total data
bs.size <- nrow(dat)/g/1e3         # size of bootstrap sample based on % of pixel number
CI.sig <- 95                       # confidence interval difference level
# METHOD 2: Based on POWER ANALYSIS
library(pwr)
sig <- 0.05                      # significance level, Type I error
alph <- 0.95                     # Power, (1-Type II error), chance of finding an effect
effect.size <- 0.1               # levels by Cohen 1992, but can vary fro field to field, 0.1 -small, 0.3 medium, 0.5 large
samp.size <- pwr.anova.test(k = g, f = effect.size, sig.level = sig, power = alph)
bs.size <- ceiling(samp.size$n)

# METHOD 1, only keep statistic ########################################################
mean.median.diff <- lapply(1:(length(which(topo.inds))*2), function(i) matrix(NA, nrow = g, ncol = g))      # mean difference in the medians between 2 groups
low.CI.diff <- lapply(1:(length(which(topo.inds))*2), function(i) matrix(NA, nrow = g, ncol = g))           # lower bound oconfidence interval between 2 groups
up.CI.diff <- lapply(1:(length(which(topo.inds))*2), function(i) matrix(NA, nrow = g, ncol = g))            # upper bound of CI 
sig.matrix <- lapply(1:(length(which(topo.inds))*2), function(i) matrix(NA, nrow = g, ncol = g))            # upper bound of CI 

dev.new()
setEPS()
postscript(paste('Bootstrap_M1_Tect_ES' ,effect.size, '.eps',sep=''), horizontal = F, pointsize = 7)
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(2,4),cex=0.5)

split.tect <- split(dat,dat$SR == 0)
split.geo.idf <- split(GeoIDf, dat$SR == 0)
for(t in 1:2){
  for (i in 1:length(which(topo.inds))){
      split.by.geo <- split(split.tect[[t]][,which(topo.inds)[i]], split.geo.idf[[t]])   # make groups
      bs.data <- lapply(1:g, function(l) vector('numeric',b))
      for (j in 1:g){                                                        # loop through list of split.by.geo
        for (k in 1:b){
          bs.sample <- sample(split.by.geo[[j]], size= bs.size, replace = T)       # bootstrap sample
          bs.data[[j]][k] <- median(bs.sample)                               # THE STATISTIC YOU WISH FOR
        }
      }
      
      for (j in 1:(g-1)){    # loop through columns of the final matrices
        for (k in (j+1):g){  #  loop through rows of lower triangle
          get.data.1 <- bs.data[[j]]    # data from column
          get.data.2 <- bs.data[[k]]    # data from row
          delta <- get.data.2 - get.data.1
          mean.median.diff[[i+length(which(topo.inds))*(t-1)]][k,j] <- mean(delta)
          low.CI.diff[[i+length(which(topo.inds))*(t-1)]][k,j] <- quantile(delta,probs = (1e2-CI.sig)/1e2)
          up.CI.diff[[i+length(which(topo.inds))*(t-1)]][k,j] <- quantile(delta,probs = (CI.sig)/1e2)
          sig.matrix[[i+length(which(topo.inds))*(t-1)]][k,j] <- findInterval(0,c(low.CI.diff[[i]][k,j],up.CI.diff[[i]][k,j]))   
          # This means 0 - row is higher than column, 1 - insignificant, 2- column higher than row
        }
      }
      
      # plot the signifcance matrices ############################################
      rotate <- function(x) t(apply(x, 2, rev))
      
      image(x =seq(1,g), y =seq(1,g), z =rotate(sig.matrix[[i+length(which(topo.inds))*(t-1)]]), main = paste(vars[which(topo.inds)[i]], tecLab[t], sep='_')
            ,xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', col = colorRampPalette(c("blue", "white", "red"))(3), zlim = c(0,2))
      axis(1,at=seq(1,g,1), labels= lab[1:g])
      axis(2,at=seq(1,g,1), labels= rev(lab[1:g]))
      grid(nx=g, ny=g, lty=1)
    }
}

dev.copy(which = 4)
dev.off(which = 4)
write.xlsx(sig.matrix, file= paste("LogPmat_M1_tect_ES", effect.size, ".xlsx",sep=''),asTable = F, sheetName= vars[topo.inds], append=TRUE, row.names=T, col.names = T)
dev.off()

# Plot Median differences and CI's in respect to ref level ##################
#############################################################################
library(ggplot2)
library(gridExtra)
ref <- 6        # choose which lithology should be the reference level
# assemble the data
means <- lapply(mean.median.diff, function(i) i[,ref])         # get the values from column
up.CI <- lapply(up.CI.diff, function(i) i[,ref])
low.CI <- lapply(low.CI.diff, function(i) i[,ref])
for (i in 1:length(means)){
  means[[i]][is.na(means[[i]])] <- mean.median.diff[[i]][ref,is.na(means[[i]])]*(-1) # get the values that werent in column from rows, but for those the signs needs to be switched
  up.CI[[i]][is.na(up.CI[[i]])] <- low.CI.diff[[i]][ref,is.na(low.CI[[i]])]*(-1) # get the values that werent in column from rows, but for those the signs needs to be switched
  low.CI[[i]][is.na(low.CI[[i]])] <- up.CI.diff[[i]][ref,is.na(low.CI[[i]])]*(-1) # get the values that werent in column from rows, but for those the signs needs to be switched
}
means <- lapply(means, function(i) i[!is.na(i)])
up.CI <- lapply(up.CI, function(i) i[!is.na(i)])
low.CI <- lapply(low.CI, function(i) i[!is.na(i)])
cols = c("#FFEBAF","#38A800","#727272","#FFAA00","#BED2FF","#0070FF","#734C00","#000000")
plotchar <- c(seq(16,18),15,15,19,17,15)
ps <- 1.5

dev.new()
setEPS()
postscript(paste('Bootstrap_M1_Tect_refplot_ES' ,effect.size, '.eps',sep=''), horizontal = F, pointsize = 7,
           width = 185/25.4, height = 90/25.4, paper = "special", onefile = FALSE)
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(2,4),cex=0.5)
plot.titles <- c(sapply(vars[which(topo.inds)], function(i) paste(i, '_act',sep='')),
                 sapply(vars[which(topo.inds)], function(i) paste(i, '_inact',sep='')))
for( t in 1:length(means)){
  plim <- max(c(abs(min(low.CI[[t]])), abs(max(up.CI[[t]]))))
  plot(seq(1:6), means[[t]],pch = plotchar[-ref], xlab='', ylab= paste(lab[ref], '_ref_level',sep=''),main = plot.titles[t], 
       col = cols[-ref], ylim = c(plim*(-1),plim), cex = 5)
  axis(1,at=seq(1,(g-1),1), labels= lab[-ref])
  arrows(seq(1,6), low.CI[[t]], seq(1,6),up.CI[[t]], length= 0, angle = 90, code = 3)
  abline(a = 0, b = 0, lty = 2)
}
dev.copy(which = 4)
dev.off(which = 4)
rm(mean.median.diff,low.CI.diff,up.CI.diff,sig.matrix,split.by.geo,bs.data,get.data.1,get.data.2,delta,sig.plot)
dev.off()

# METHOD 2, calculate statistic directly on whole distribution #########################
wilcox.p.perc <- lapply(1:(length(which(boot.inds))*2), function (i) i)
dev.new()
setEPS()
postscript(paste('Bootstrap_M2_tect_ES', effect.size, '.eps'), horizontal = F, pointsize = 7)
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(2,4),cex=0.5)

split.tect <- split(dat,dat$SR == 0)                                         # split into active vs inactive
split.geo.idf <- split(GeoIDf,dat$SR == 0)
tecLab <- c('act','inact')
for(t in 1:2){
  for (i in 1:length(which(topo.inds))){
    split.by.geo <- split(split.tect[[t]][,which(topo.inds)[i]], split.geo.idf[[t]])   # make groups
    wilcox.p <- lapply(1:b, function(l) matrix(NA, nrow = g-1, ncol = g-1))
    for (j in 1:b){
      bs.sample <- list()
      bs.sample <- lapply(1:g, function (z) sample(split.by.geo[[z]], size= bs.size, replace = T))   # take a bootstrap sample for every GEO
      bs.sample <- unlist(bs.sample)
      bs.group <- rep(seq(1,g), each = bs.size)
      wilcox <- pairwise.wilcox.test(bs.sample, bs.group, p.adjust.method = 'BH')  # Do pairwise wilcox test to test for differences in the distributions
      wilcox.p[[j]] <- wilcox$p.value                                              # save only p-values
    }
    wilcox.binary <- lapply(wilcox.p, function(z) z < sig)                         # make binary matrices with significant or not
    wilcox.sum <- Reduce('+',wilcox.binary)                                        # sum them all up
    wilcox.p.perc[[i]] <- wilcox.sum/b                                             # Calculate percentage of significant outcome
    
    # plot the signifcance matrices ############################################
    rotate <- function(x) t(apply(x, 2, rev))
    
    image(x =seq(1,6), y =seq(1,6), z =rotate(wilcox.p.perc[[i]]), main = paste(vars[which(topo.inds)[i]], tecLab[t], sep='_'),
          xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', col = colorRampPalette(c("blue","white","red"))(100), zlim = c(0,1))
    axis(1,at=seq(1,6,1), labels= lab[1:(g-1)])
    axis(2,at=seq(1,6,1), labels= rev(lab[2:g]))
    grid(nx=(g-1), ny=(g-1), lty=1)
  }
}
dev.copy(which = 4)
dev.off(which = 4)
write.xlsx(wilcox.p.perc, file= paste("Pmat_M2_tect_ES", effect.size, ".xlsx",sep=''),asTable = F, sheetName= vars[boot.inds], append=TRUE, row.names=T, col.names = T)

rm(split.by.geo,bs.data,wilcox.p,wilcox.binary,wilcox.p.perc,wilcox,w.plot)
dev.off()

################################################################
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
