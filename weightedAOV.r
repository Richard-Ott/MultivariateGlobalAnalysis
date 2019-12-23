# My very first R Script!!!!!!
# Richard Ott, 2019
rm(list = ls())
mthres <- 50
setwd(paste('C:/Richard/PhD_ETH/matlab/Multivariate Analysis/', mthres, 'm/for_stat_consulting',sep=''))
library(DescTools)
library(nnet)
library(graphics)
library(car)
library(ggplot2)
library(plyr)
library(nlme)


# LOAD DATA ############################################################################
vars <- c("Elev","Locrel","slope","P",'T','NDVI','ksn','tetrapods'
          ,'Amphibians','Mammals', 'KG','Geo','LC', 'SR', 'weight')
dat <- read.table(paste('global_poly_',mthres,'mLC.txt',sep=''), header = FALSE, sep = ",", col.names	= vars)
dat$Amphibians[dat$Amphibians == -9999 | dat$Amphibians == -32768] = 0   # Amph have lots of nan values that should be set to zero for analysis
corr = 1 # Do you want biological variables to be corrected for T and P influence?

# MODIFY DATA ########################################################################
depVar <- which(vars == 'Geo')                     # dependent variable, should be last variable
#         E Loc Sl P  T  ND ksn te  am  ma KG Geo LC SR  area NDVIc Amphc Tetc
inds <- c(T, T, T, T, F, F,  T, F  ,F , F, F, F , F , T,  F,   T     ,T,   T)

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
  Tet.corr <- dat$tetrapods - Tet.emp
  dat$tet_corr <- Tet.corr
  
  rm(linMod,lm_results,Amph.corr,Amph.emp,NDVIcorr,NDVIemp,Tet.emp,Tet.corr,Pthres)
  vars[16:18] <- c('NDVIc','Amphc','Tetc')
}

# group and assign factor for landcover data ###########################################
lcID <- vector()
lcID[is.element(dat$LC,c(1,2,14))] <-  1    # bFor, broadleaf forest
lcID[is.element(dat$LC,c(3,4))] <-  2       # nFor, needleleaf forest
#lcID[dat$LC == 5] <- 3                     # mFor, mixed forest
lcID[is.element(dat$LC,c(6,7,8,9))] <-  3   # Medium Vegetation
lcID[is.element(dat$LC,c(10,16,17))] <-  4  # Sparse/no Vegetation
lcID[is.element(dat$LC,c(11,12))] <-  5     # Cropland, 13 seems to be likely cropland but it seemed weird
lcID[dat$LC == 18] <- 6                     # Urban
#lcID[dat$LC == 19] <- 8                    # Snow/Ice, so small can be left out
#lcID[dat$LC == 20] <- 9                    # Water, so small can be left out
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
write.table(stats.by.group, paste("MeanStats_",mthres,"poly_m.csv"), append = FALSE, sep = ",", dec = ".",
            row.names = T, col.names = T)

##################################################################
## PLOT DATA #####################################################


# BOXPLOT OF ALL DATA #################################################################
#######################################################################################
library(epade)
dev.new()
setEPS()
postscript("boxplot_all_corrected_poly.eps")
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(3,3),cex=0.5,lwd=0.35)
cols = c("#FFEBAF","#38A800","#727272","#FFAA00","#BED2FF","#0070FF","#734C00","#000000")
# dev.set(which == 2)
weight.by.geo <- split(dat$weight,GeoIDf)
for (i in 1:length(vars[inds])){ # loop through variables
  # rearrange data matrix
  data.by.geo <- split(dat[,which(inds)[i]], GeoIDf)     # split by geology
  weighted.by.geo <-list()
  for (j in 1:g){
    weighted.by.geo[[j]] <- rep(data.by.geo[[j]],weight.by.geo[[j]])     # repeat every value as many times as weight, to make boxplot correct
  }
  # box.plot.wtd(dat[,which(inds)[i]],GeoIDf, w = dat$weight, vnames= lab, col = cols ,main = vars[inds][i])
  boxplot(weighted.by.geo,names = lab ,main = vars[inds][i], outline = F, boxlty = 1, whisklty = 1, staplelty = 0, col = cols)
}
dev.copy(which = 4)
dev.off(which = 4)
dev.off()

# TOPO-BOXPLOTS TECTONICALLY INACTIVE VS ACTIVE #######################################
#######################################################################################
Tect <- dat$SR == 0              # inactive = True, active = False
weight.by.tect <- split(dat$weight,Tect) # split weights by tectonic activity
Geo.by.tect <- split(GeoIDf,Tect)
act.geo.weights <- split(weight.by.tect[[1]],Geo.by.tect[[1]])    # weights for active areas, split by geo 
in.geo.weights <- split(weight.by.tect[[2]],Geo.by.tect[[2]])

dev.new()
setEPS()
postscript("boxplot_tectonic_poly.eps")
dev.set(which = 2)
par(oma=c(0,0,2,0),mfrow=c(2,4),cex=0.35,lwd=0.35)
cols = c("#FFEBAF","#38A800","#727272","#FFAA00","#BED2FF","#0070FF","#734C00","#000000")
# dev.set(which == 2)
#              E Loc Sl P  T  ND ksn te  am  ma KG Geo LC SR  area NDVIc Amphc Tetc
topo.inds <- c(T, T, T, F, F, F,  T, F  ,F , F, F, F , F , F,  F,   F     ,F,   F)
for (i in 1:length(which(topo.inds))){ # loop through variables
  # rearrange data matrix
  data.by.tect <- split(dat[,which(topo.inds)[i]], Tect)     # split by tectonics
  
  # first do active areas
  data.by.geo <- split(data.by.tect[[1]],Geo.by.tect[[1]])        # split by GEO
  weighted.by.tect <-list()
  for (j in 1:g){
    weighted.by.tect[[j]] <- rep(data.by.geo[[j]],act.geo.weights[[j]])     # repeat every value as many times as weight, to make boxplot correct
  }
  boxplot(weighted.by.tect,names = lab ,main = paste(vars[topo.inds][i], '_act', sep='')
          , outline = F, boxlty = 1, whisklty = 1, staplelty = 0, col = cols, cex.axis = 2)
  
  # now do inactive areas
  rm(data.by.geo)
  data.by.geo <- split(data.by.tect[[2]],Geo.by.tect[[2]])        # split by GEO
  weighted.by.tect <-list()
  for (j in 1:g){
    weighted.by.tect[[j]] <- rep(data.by.geo[[j]],in.geo.weights[[j]])     # repeat every value as many times as weight, to make boxplot correct
  }
  lims <- par("usr")[3:4]
  boxplot(weighted.by.tect,names = lab ,main = paste(vars[topo.inds][i], '_in', sep='')
          , outline = F, boxlty = 1, whisklty = 1, staplelty = 0, col = cols,cex.axis = 2, ylims = lims)
}
dev.copy(which = 4)
dev.off(which = 4)
dev.off()
rm(weight.by.tect,Geo.by.tect,act.geo.weights,in.geo.weights,topo.inds,data.by.geo,lims,data.by.tect,weighted.by.tect)

# BOXPLOT PER CLIMATE ZONE ###########################################################
######################################################################################

#             E Loc Sl P  T  ND ksn te  am  ma KG Geo LC SR  area NDVIc Amphc Tetc
bio.inds <- c(F, F, F, F, F, T,  F,  T  ,T , F, F, F , F , F,  F,   F     ,F,   F)

dev.new()
KGlab = c("tropical", "dry", "temp_cont","polar")
for (k in 1:4){
  # setEPS()
  # postscript(paste("boxplotKG_", KGlab[k] ,".eps", sep  = "", collapse = NULL))
  # dev.set(which = 2)
  par(oma=c(0,0,2,0),mfrow=c(3,4),cex=0.5)
  
  if (k == 3){   # This makes climate zones temperate and continental merge together
    dat.kg <-dat[dat$KG == 3 | dat$KG == 4 ,bio.inds]
    GeoKG <- GeoIDf[dat$KG == 3 | dat$KG == 4]
    stat.kg <- dat[dat$KG == 3 | dat$KG == 4 ,]
  }else{
    dat.kg <-dat[dat$KG == k,bio.inds]
    GeoKG <- GeoIDf[dat$KG == k]
    stat.kg <- dat[dat$KG == k,]
  }
  
  # do statistics on this KG ########################
  # stats.by.group <- array(NA,dim=c(length(vars),g*2))
  # for (i in 1:length(vars)){
  #   stats.by.group[i,1:g]<- tapply(stat.kg[,i],GeoKG,mean)
  #   stats.by.group[i,seq(g+1,g*2)]<- tapply(stat.kg[,i],GeoKG,median)
  # }
  # rownames(stats.by.group) <- vars
  # colnames(stats.by.group) <- c(lab,lab)
  # write.table(stats.by.group, paste("MeanKG_", KGlab[k] ,".csv", sep  = "", collapse = NULL), append = FALSE, sep = ",", dec = ".",
  #             row.names = T, col.names = T)
  ####################################################
  
  for (i in 1:length(which(bio.inds))){
    # rearrange data matrix
    areas <- vector()
    dat.by.geo <-list()
    for (j in 1:g){
      Ginds = GeoKG == j
      dat.by.geo[[j]] <- dat.kg[Ginds,i]
      areas[j] <- length(Ginds[Ginds == T])/nrow(dat.kg)
    }
    boxplot(dat.by.geo, xaxt = "n", main = vars[bio.inds][i], outline = F, boxlty = 1, whisklty = 1, staplelty = 0, col = cols)
    axis(1, at=1:g, labels=lab)
  }
  # pie(areas, labels = lab)
  # title(main= KGlab[k] , outer=T)
  # dev.copy(which = 4)
  # dev.off(which = 4)
}
rm(Ginds,dat.kg,areas,GeoKG)

## P-T-RASTER FOR BIO-VARIABLES #############################################################
#############################################################################################
library(Hmisc)
library(mapplots)
library(fields)
library(spatstat)
T.quart <- wtd.quantile(dat$T, weights = dat$weight)   # determine quartile boundaries
P.quart <- wtd.quantile(dat$P, weights = dat$weight)

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

# IMPORTANT: I only plot data were there are more samples than a certain threshold
plot.thres <- 1e3                        # below this number of samples a certain T-P- region is not viewed as sampled enough

# Calculate the means for all bio variables in the TP-grids
all.means <- list()
diff.means <- list()
zlims <- list('1' = c(-25,25), '2'= c(-80,80), '3' = c(-12,12))
dev.new()
for (i in 1:length(which(bio.inds))){     # loop through all bio variables
  setEPS()
  postscript(paste('TP_', vars[which(bio.inds)[i]] ,'.eps', sep=''))
  dev.set(which = 2)
  
  all.means[[i]] <- list()
  diff.means[[i]] <- list()
  # Calculate total mean
  var <- dat[,which(bio.inds)[i]]
  split.var <- split(var,TP.val)
  split.weight <- split(dat$weight,TP.val)
  TP.frequencies <- sapply(split.weight, function(i) sum(i))
  split.means <- sapply(1:length(split.weight), function(i)
    wtd.mean(split.var[[i]], split.weight[[i]]))     # can change this to median
  # split.means <- sapply(1:length(split.weight), function(i) 
  #   weighted.median(split.var[[i]], split.weight[[i]]))     # can change this to median
  all.means[[i]][[g+1]] <- matrix(split.means,4,4, byrow = T)
  TP.frequencies <- matrix(TP.frequencies,4,4, byrow = T)
  
  par(oma=c(0,0,2,0),mfrow=c(3,3),cex=0.5)
  for (j in 1:g){                                  # loop through GEO's
    var <- dat[GeoIDf == j, which(bio.inds)[i]]    # get values of this lith and variable
    split.var <- split(var, TP.val[GeoIDf == j])   # split by TP-index
    id <- sapply(split.var, function(i) length(i))
    id <- id < plot.thres
    split.weight <- split(dat$weight[GeoIDf == j], TP.val[GeoIDf == j])
    split.means <- sapply(1:length(split.weight), function(i)
      wtd.mean(split.var[[i]], split.weight[[i]]))     # can change this to median
    # split.means <- sapply(1:length(split.weight), function(i) 
    #   weighted.median(split.var[[i]], split.weight[[i]]))     # can change this to median
    split.means[id] <- NA                              # This removes regions that dont have enough data from the plot
    split.means <- matrix(split.means,4,4, byrow = T)
    
    all.means[[i]][[j]] <- split.means
    diff.means[[i]][[j]] <- split.means - all.means[[i]][[g+1]]
    # values below or above the z-limit should be shown by th max color. By default R will leave those values out so I will correct them for the image process
    diff.means[[i]][[j]][which(diff.means[[i]][[j]] > zlims[[i]][2])] <- zlims[[i]][2] # replace above maximum zlim values with max zlim
    diff.means[[i]][[j]][which(diff.means[[i]][[j]] < zlims[[i]][1])] <- zlims[[i]][1]
    
    image(x =seq(1,4), y =seq(1,4), z =diff.means[[i]][[j]], main = lab[j],xaxt = 'n', yaxt = 'n', zlim = zlims[[i]],
          xlab = 'T', ylab = 'P', col = tim.colors())
    axis(1,at=seq(0.5,4.5,1), labels= as.character(T.quart/10))
    axis(2,at=seq(0.5,4.5,1), labels= as.character(P.quart))
  }
  # draw frequencies of all data
  image(x =seq(1,4), y =seq(1,4), z = TP.frequencies, main = 'sample frequencies',xaxt = 'n', yaxt = 'n', xlab = 'T', ylab = 'P', col = tim.colors())
  axis(1,at=seq(0.5,4.5,1), labels= as.character(T.quart/10))
  axis(2,at=seq(0.5,4.5,1), labels= as.character(P.quart))
  
  dev.copy(which = 4)
  dev.off(which = 4)
}
rm(split.var,split.weight,diff.means,all.means,zlims,T.quart,P.quart,var,TP.val,T.intervals,P.intervals)

#############################################################################################
#############################################################################################
# REGRESSION AND MORE STATISTICS #######################################################


# Do weighted ANOVA #####################################################
# do ANOVA for all variables manually check output
dir.create("./AOV")
setwd("./AOV")
library(mvtnorm)
library(zoo)
library(multcomp)
library(openxlsx)
library(broom)
anov = list()
levTest = list()
anova.sum = list()
lev.sum = list()
glht.sum = list()
shap = list()
for (i in 1:length(vars[inds])){
  anov[[i]] = list()
  anova.sum[[i]] = list()
  levTest[[i]] = list()
  lev.sum[[i]] = list()
  glht.sum[[i]] =list()
  shap[[i]] = list()
  # Do ANOVA, ghlt, Levene test
  mod <- lm(dat[,which(inds)[i]] ~ GeoIDf, weights = dat$weight)
  anov[[i]] <- anova(mod)
  # anov[[i]] <- aov(dat[,vars[inds][i]] ~ GeoIDf)
  
  # Levene Test checks the AOV assumption of homogenoeus variance, p should be > 0.05! 
  # to show that the variances are not significantly different between groups
  levTest[[i]] <- leveneTest(dat[,which(inds)[i]] ~ GeoIDf) 
  
  # Turkey HSD test shows between which groups the differences are significant, using a linear model
  glht.sum[[i]] <- tidy(summary(glht(model = mod, linfct = mcp(GeoIDf= "Tukey"))))
  
  # Extract the residuals
  aov_residuals <- residuals(object = anov[[i]] )
  # Run Shapiro-Wilk test to check AOV normality assumption, p should be > 0.05 for normality to be probably not violated
  # shap[[i]] = shapiro.test(x = aov_residuals ) # only possible for < 5000 samples
  
  # save results
  anova.sum[[i]] <- as.data.frame(summary(anov[[i]]))
  lev.sum[[i]] <- as.data.frame(summary(levTest[[i]]))
  
  # save results as excel
  # write.xlsx(as.data.frame(anova.sum[[i]]), file="aov.xlsx", asTable = F,sheetName= vars[inds][i], append=TRUE, row.names=FALSE)
  # write.xlsx(as.data.frame(lev.sum[[i]]), file="lev.xlsx",asTable = F, sheetName= vars[inds][i], append=TRUE, row.names=FALSE)
  # write.xlsx(tidy(glht.sum[[i]]), file="tukey.xlsx",asTable = F, sheetName= vars[inds][i], append=TRUE, row.names=FALSE)
}
write.xlsx(anova.sum, file="aov.xlsx", asTable = F,sheetName= vars[inds], append=TRUE, row.names=FALSE)
write.xlsx(lev.sum, file="lev.xlsx",asTable = F, sheetName= vars[inds], append=TRUE, row.names=FALSE)
write.xlsx(glht.sum, file="tukey.xlsx",asTable = F, sheetName= vars[inds], append=TRUE, row.names=FALSE)

# %SR per GEO #########################################################################################
split.by.geo <- split(dat$SR,GeoIDf)
stat.by.geo <- lapply(split.by.geo, function(i) length(which(i > 0))/length(i))
stat.by.geo <- unlist(stat.by.geo)
names(stat.by.geo) <- lab
write.table(stat.by.geo, paste("SRperc",mthres,"m.csv"), append = FALSE, sep = ",", dec = ".",
            row.names = T, col.names = F)


# KOLMOGOROW-SMIRNOV-TEST & 2-SAMPLE KUIPER TEST ############################################
library("remotes")
library("nonpar")
library("kuiper.2samp")
ks.stats <- list()
cc.stats <- list()
kui.stats <- list()
for (i in 1:length(vars[inds])){    # Loop through variables
  Test.data <- dat[,which(inds)[i]] # data for this variable
  data.split <- split(Test.data,GeoIDf)   # split data by Geo
  ks.stats[[i]] <- matrix(, nrow = g, ncol = length(data.split)-1)
  cc.stats[[i]] <- matrix(, nrow = g, ncol = length(data.split)-1)
  kui.stats[[i]] <- matrix(, nrow = g, ncol = length(data.split)-1)
  for (j in 1:g){                   # Loop through different rock types
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