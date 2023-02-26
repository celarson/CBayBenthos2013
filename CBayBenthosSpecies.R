#species level analysis
#use output of CBayAmbiTaxa code

#Data analysis for species level taxonomy on the CBay 2013 benthos dataset

#Load packages
library(fossil)
library(tidyverse)
library(vegan)
library(ecodist)
library(indicspecies)
library(ggcorrplot)
library(psych)
library(polycor)
library(caret)
library(ape)
library(randomForest)
library(party)
library(terra)
library(gdm)
library(raster)

#functions
directionalSAC <- function(community, gradient) {
  require(vegan)
  if (!inherits(community, "data.frame")){
    if(!inherits(community, "matrix"))
      stop("Object community must be of class data.frame or matrix")
  }
  if (!inherits(gradient, "dist")){
    if(!is.vector(gradient)){
      if(!is.matrix(gradient)){
        if(!is.data.frame(gradient))
          stop("Object gradient must be of class vector, matrix, data.frame or
dist")
      }}}
  if (any(community < 0))
    stop("Negative value in community")
  if (any(rowSums(community) < 1e-16))
    stop("Empty plots in community")
  if(is.null(rownames(community))) stop("Object community must have row names")
  if(is.vector(gradient)){
    if(length(gradient)!=nrow(community)) stop("Incorrect definition of object
gradient")
    if(!is.null(names(gradient))){
      if(any(!rownames(community)%in%names(gradient))) stop("Names in gradient
must be the same as row names in community")
      gradient <- gradient[rownames(community)]
    }
    Rgradient <- rownames(community)[order(gradient)]
    agg <- community[Rgradient, ]
    3
    richness <- specnumber(agg)
    average_alfa <- cumsum(richness)/(1:length(richness))
    average <- specaccum(agg, method="collector")$richness
  }
  else {
    if(!is.matrix(gradient)) gradient <- as.matrix(gradient)
    if(ncol(gradient)==1) stop("A gradient with a single quantitative variable
must be coded as a vector rather than as a matrix")
    if(nrow(gradient)!=nrow(community)) stop("Incorrect definition of object
gradient")
    if(!is.null(rownames(gradient))){
      if(any(!rownames(community)%in%rownames(gradient))) stop("Row names in
gradient must be the same as row names in community")
      gradient <- gradient[rownames(community), ]
    }
    res <- array(NA, c(ncol(gradient), nrow(gradient)))
    for(i in 1:ncol(gradient)){
      nami <- rownames(community)
      res[i, ] <- nami[order(gradient[, i])]
    }
    spatial_order <- t(res)
    f <- nrow(spatial_order)
    n <- ncol(spatial_order)
    result <- array(dim = c(f, n))
    alfa_average <- array(dim = c(f, n))
    for(i in 1:n) {
      agg <- community[spatial_order[, i], ]
      richness <- specnumber(agg)
      alfa_s <- cumsum(richness)/(1:length(richness))
      c <- specaccum(agg, method="collector")
      result[, i] <- c$richness
      alfa_average[, i] <- alfa_s
    }
    average <- rowMeans(result)
    average_alfa <- rowMeans(alfa_average)
  }
  beta_s <- average/average_alfa
  beta_S <- (beta_s-1)/((1:length(beta_s))-1)
  exact <- specaccum(community, method = "exact")
  beta_exact <- exact$richness/exact$richness[1]
  beta_N <- (beta_exact-1)/((1:length(beta_exact))-1)
  beta_norm_autocor <- (beta_N- beta_S)/(beta_N+ beta_S)
  SCR <- data.frame(as.matrix(average), as.matrix(exact$richness),
                    as.matrix(average_alfa), as.matrix(beta_s), as.matrix(beta_S),
                    as.matrix(beta_exact), as.matrix(beta_N), as.matrix(beta_norm_autocor))
  names(SCR) <- c("N_SCR", "N_Exact", "Alpha_dir", "Beta_M_dir",
                  "Beta_N_dir","Beta_M", "Beta_N", "Beta_Autocor")
  return(SCR)
}

# Calculate a pairwise association between all variables in a data-frame. In particular nominal vs nominal with Chi-square, numeric vs numeric with Pearson correlation, and nominal vs numeric with ANOVA.
# Adopted from https://stackoverflow.com/a/52557631/590437
mixed_assoc = function(df, cor_method="spearman", adjust_cramersv_bias=TRUE){
  df_comb = expand.grid(names(df), names(df),  stringsAsFactors = F) %>% set_names("X1", "X2")
  
  is_nominal = function(x) class(x) %in% c("factor", "character")
  # https://community.rstudio.com/t/why-is-purr-is-numeric-deprecated/3559
  # https://github.com/r-lib/rlang/issues/781
  is_numeric <- function(x) { is.integer(x) || is_double(x)}
  
  f = function(xName,yName) {
    x =  pull(df, xName)
    y =  pull(df, yName)
    
    result = if(is_nominal(x) && is_nominal(y)){
      # use bias corrected cramersV as described in https://rdrr.io/cran/rcompanion/man/cramerV.html
      cv = cramerV(as.character(x), as.character(y), bias.correct = adjust_cramersv_bias)
      data.frame(xName, yName, assoc=cv, type="cramersV")
      
    }else if(is_numeric(x) && is_numeric(y)){
      correlation = cor(x, y, method=cor_method, use="complete.obs")
      data.frame(xName, yName, assoc=correlation, type="correlation")
      
    }else if(is_numeric(x) && is_nominal(y)){
      # from https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
      r_squared = summary(lm(x ~ y))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else if(is_nominal(x) && is_numeric(y)){
      r_squared = summary(lm(y ~x))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else {
      warning(paste("unmatched column type combination: ", class(x), class(y)))
    }
    
    # finally add complete obs number and ratio to table
    result %>% mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), complete_obs_ratio=complete_obs_pairs/length(x)) %>% rename(x=xName, y=yName)
  }
  
  # apply function to each variable combination
  map2_df(df_comb$X1, df_comb$X2, f)
}

#Color vectors
Zone_col_vec<-c("#e41a1c","#377eb8","#4daf4a","#984ea3")
veg_col_vec<-c("#99d8c9","#41ae76","#006d2c","#00441b","black")
Depth_col_vec<-c("#d0d1e6","#3690c0","#023858")
Marina4_col_vec<-c("#fcbba1","#cb181d","#67000d")

#Upload and refine datasets
#CBayAPISbenthos Chao1
CBayspChao1<-read.csv("output/all/ezD7782405_APTC_SG1_singdoub_Lymn.csv",header=T)
names(CBayspChao1)
#change dash to underscore for site names
CBayspChao1$InvEventFK<- gsub("_", "-", CBayspChao1$InvEventFK)
CBayspChao1$InvEventFK<- gsub("CBay", "CB", CBayspChao1$InvEventFK)
#CBayAPISbenthos Chao2
CBayspChao2<-read.csv("output/all/ezD7782400_APTC_SG1_unidup_Lymn.csv",header=T)
names(CBayspChao2)
#change dash to underscore for site names
CBayspChao2$InvEventFK<- gsub("_", "-", CBayspChao2$InvEventFK)
CBayspChao2$InvEventFK<- gsub("CBay", "CB", CBayspChao2$InvEventFK)
#Environment
CBayAPISEnv<-read.csv("ezD7782411_CBay_Env_Inf.csv",header=T)
names(CBayAPISEnv)
summary(CBayAPISEnv)
#rename site code key to match short name in benthos dataset
CBayAPISEnv<-rename(CBayAPISEnv, "InvEventFK"="SiteCodeKey")
#make factors
CBayAPISEnv$ZoneCode<-as.factor(CBayAPISEnv$ZoneCode)
#change order of sediment categories
CBayAPISEnv$SedComb<-as.factor(CBayAPISEnv$SedComb)
CBayAPISEnv$SedComb<-factor(CBayAPISEnv$SedComb,levels = c("organic_plus", "clay_silt", "sand_only","sand_plus","hard"))
#make row name site code for new env data frame only CBay and create zone code variable
CBayEnvrn<-subset(CBayAPISEnv, ZoneCode!="random" & ZoneCode!="target" & SiteComments!="Matched to sled")
CBayEnvrn$ZoneCode<-droplevels(CBayEnvrn$ZoneCode)
levels(CBayEnvrn$ZoneCode)
#combine zones so no band
CBayEnvrn$Zone<-CBayEnvrn$ZoneCode
levels(CBayEnvrn$Zone) <- list("Zone1"=c("Z1_BandNO", "Z1_BandYES"), 
                               "Zone2"=c("Z2_BandNO", "Z2_BandYES"),
                               "Zone3"=c("Z3_BandNO", "Z3_BandYES"),
                               "Zone4"=c("Z4_BandNO", "Z4_BandYES"))
levels(CBayEnvrn$Zone)
rownames(CBayEnvrn)<-CBayEnvrn$InvEventFK
#upload distance from marina
CBayMarina<-read.csv("ezD9001767_MarinaDistancesCBay.csv",header=T)
#merge with envirnmental
CBayEnvrnm<-merge(CBayEnvrn,CBayMarina,by="InvEventFK")

#sediment no organic switched
CBaySed<-read.csv("ezD7782420_CBaySedWideInf.csv",header=T)
names(CBaySed)
CBaySed<-rename(CBaySed, "InvEventFK"="SiteCodeKey")
#merge with environmental variabes
CBaySedEnv<-merge(CBaySed,CBayEnvrnm,by="InvEventFK")
rownames(CBaySedEnv)<-CBaySedEnv$InvEventFK
#sediment organic switched
CBaySedorg<-read.csv("ezD7782421_CBaySedWideInf_orgsplit.csv",header=T)
CBaySedorg<-rename(CBaySedorg, "InvEventFK"="SiteCodeKey")
#merge with environmental variabes
CBaySedEnvorg<-merge(CBaySedorg,CBayEnvrnm,by="InvEventFK")
rownames(CBaySedEnvorg)<-CBaySedEnvorg$InvEventFK
#upload distances between sites
CBayDist<-read.csv("ezD7562779_Distance_calulations_4_27_22.csv",header=T,row.names =1)
names(CBayDist)<- gsub("\\.", "-", names(CBayDist))
#Upload projected coordinates
CBayCoord<-read.csv("ezD7971880_CBay_Coords.csv",header=T)

#Create site by taxa matrix Chao1
names(CBayspChao1)
CBayPonarMatrixCh1<-data.frame(t(create.matrix(
  CBayspChao1,
  tax.name="newtaxa",
  locality="InvEventFK",
  time.col=NULL,
  time=NULL, 
  abund=TRUE,
  abund.col="InvCount.sum"
)))
#273 events and 253 species captured
#write csv to folder
write.csv(CBayPonarMatrixCh1,"CBayPonarMatrixChao1.csv", row.names=TRUE)

#Combine matrix and environmental parameters
CBayPonarMatrixCh1Env<-merge(CBayPonarMatrixCh1,CBaySedEnv,by=0)
names(CBayPonarMatrixCh1Env)

#reduce variable space in environmental parameters
#Water chemistry is incomplete, so just use depth, veg, and sediment
#check variances - don't use if variance is close to 0
var(CBayPonarMatrixCh1Env[,c(256:261,272,277,295)])
#rock has variance of 0.003, so don't use. Rest are >0.2
#Doesn't need to be normal, but assumes linear relationship. Plot to see if linear
plot(CBayPonarMatrixCh1Env[,c(256:260,272,277,295)])
#linear assumption looks okay
envcor<-data.frame(abs(cor(CBayPonarMatrixCh1Env[,c(256:260,272,277,295)])))
#use correlation cut off of |0.5|
#Sand correlated with clay (0.85)
#depth correlated with sand (0.57)
#depth correlated with clay (0.69)
diag(envcor) <- NA
mean(envcor$clay)
#0.33
mean(envcor$sand)
#0.37
mean(envcor$HabDepth)
#0.30
#delete sand and clay, leave depth
write.csv(CBayPonarMatrixCh1Env[,c(1:254,257:259,265:266,272,277,295)],"CBayPonarMatrixChao1Env.csv", row.names=F)

#run principal componenet
pcaenv<- prcomp(CBayPonarMatrixCh1Env[,c(257:259,272,277,295)], scale. = TRUE)
write.csv(pcaenv$x,"CBayEnvPCA.csv", row.names=TRUE)
pca_1_2 <- data.frame(pcaenv$x[, 1:2])
plot(pcaenv$x[,1], pcaenv$x[,2])
pca_var <- pcaenv$sdev^2
pca_var_perc <- round(pca_var/sum(pca_var) * 100, 1)
barplot(pca_var_perc, main = "Variation Plot", xlab = "PCs", ylab = "Percentage Variance", ylim = c(0, 100))
biplot(pcaenv)

#try again with organic split
#Combine matrix and environmental parameters
CBayPonarMatrixCh1Envorg<-merge(CBayPonarMatrixCh1,CBaySedEnvorg,by=0)
names(CBayPonarMatrixCh1Envorg)

#reduce variable space in environmental parameters
#Water chemistry is incomplete, so just use depth, veg, and sediment
#check variances - don't use if variance is close to 0
var(CBayPonarMatrixCh1Envorg[,c(256:262,273,278,296)])
#rock and wood have very small variance, so don't use. Rest are >0.2
envcororg<-data.frame(abs(cor(CBayPonarMatrixCh1Envorg[,c(256:257,259:261,273,278,296)])))
#Sand and depth correlated with clay
#delete sand and clay
#run principal componenet
pcaenvorg<- prcomp(CBayPonarMatrixCh1Envorg[,c(257,259:260,273,278,296)], scale. = TRUE)
pca_1_2org <- data.frame(pcaenvorg$x[, 1:2])
plot(pcaenvorg$x[,1], pcaenvorg$x[,2])
pca_varorg <- pcaenvorg$sdev^2
pca_var_percorg <- round(pca_varorg/sum(pca_varorg) * 100, 1)
barplot(pca_var_percorg, main = "Variation Plot", xlab = "PCs", ylab = "Percentage Variance", ylim = c(0, 100))
biplot(pcaenvorg)

#didn't make much difference, so just use the datset where organic is how originally coded
#merge PC coordinates with matrix
CBayPonarMatrixCh1Env$Row.names<-NULL
CBayPonarMatrixCh1EnvPC<-merge(CBayPonarMatrixCh1Env,pcaenv$x,by=0)

#subset into zones
levels(CBayPonarMatrixCh1Env$Zone)
CBayPonMatCh1EnvZon<-split(CBayPonarMatrixCh1Env,f=CBayPonarMatrixCh1Env$Zone)
#49 in zone 1, 76 in zone 2, 61 in zone 3, 87 in zone 4
#Create species matrices for each zone separately
CBayPonMatCh1EnvZon1<-CBayPonMatCh1EnvZon[[1]]
row.names(CBayPonMatCh1EnvZon1)<-CBayPonMatCh1EnvZon1$Row.names
CBayPonMatCh1Zon1<-CBayPonMatCh1EnvZon1[,c(0,1:253)]
CBayPonMatCh1Zon1<-CBayPonMatCh1Zon1[, which(colSums(CBayPonMatCh1Zon1) != 0)]
#168 species found in zone 1
CBayPonMatCh1EnvZon2<-CBayPonMatCh1EnvZon[[2]]
row.names(CBayPonMatCh1EnvZon2)<-CBayPonMatCh1EnvZon2$Row.names
CBayPonMatCh1Zon2<-CBayPonMatCh1EnvZon2[,c(0,1:253)]
CBayPonMatCh1Zon2<-CBayPonMatCh1Zon2[, which(colSums(CBayPonMatCh1Zon2) != 0)]
#170 species found in zone 2
CBayPonMatCh1EnvZon3<-CBayPonMatCh1EnvZon[[3]]
row.names(CBayPonMatCh1EnvZon3)<-CBayPonMatCh1EnvZon3$Row.names
CBayPonMatCh1Zon3<-CBayPonMatCh1EnvZon3[,c(0,1:253)]
CBayPonMatCh1Zon3<-CBayPonMatCh1Zon3[, which(colSums(CBayPonMatCh1Zon3) != 0)]
#153 species found in zone 3
CBayPonMatCh1EnvZon4<-CBayPonMatCh1EnvZon[[4]]
row.names(CBayPonMatCh1EnvZon4)<-CBayPonMatCh1EnvZon4$Row.names
CBayPonMatCh1Zon4<-CBayPonMatCh1EnvZon4[,c(0,1:253)]
CBayPonMatCh1Zon4<-CBayPonMatCh1Zon4[, which(colSums(CBayPonMatCh1Zon4) != 0)]
#153 species found in zone 4

#Species accumulation curve
#for all sites
CBaySpecAcCh1<-specaccum(CBayPonarMatrixCh1, method = "exact", conditioned =TRUE, 
                      gamma = "jack1",  w = NULL)
plot(CBaySpecAcCh1)
CBaypoolCh1<-poolaccum(CBayPonarMatrixCh1)
plot(CBaypoolCh1)
#max chao1 311.299, so S50=155.65 (reached by 14 samples) and S95=295.73 (reached by 160 samples)

#now repeat with subsetted datasets
#do zone 4 first because it has the most sites
CBaySpecAc4Ch1<-specaccum(CBayPonMatCh1Zon4, method = "exact", conditioned =TRUE, 
                       gamma = "jack1",  w = NULL)
plot(CBaySpecAc4Ch1,col="#984ea3",ylim=c(0,180))
CBaySpecAc1Ch1<-specaccum(CBayPonMatCh1Zon1, method = "exact", conditioned =TRUE, 
                       gamma = "jack1",  w = NULL)
plot(CBaySpecAc1Ch1,add=T,col="#e41a1c")
CBaySpecAc2Ch1<-specaccum(CBayPonMatCh1Zon2, method = "exact", conditioned =TRUE, 
                       gamma = "jack1",  w = NULL)
plot(CBaySpecAc2Ch1,add=T,col="#377eb8")
CBaySpecAc3Ch1<-specaccum(CBayPonMatCh1Zon3, method = "exact", conditioned =TRUE, 
                       gamma = "jack1",  w = NULL)
plot(CBaySpecAc3Ch1,add=T,col="#4daf4a")

#Try again using individuals rather than sites
sum(CBayPonarMatrixCh1)
#94223 individuals
CBaySpecAciCh1<-specaccum(CBayPonarMatrixCh1, method = "rarefaction", conditioned =TRUE,
                       gamma = "jack1",  w = NULL)
plot(CBaySpecAciCh1, xvar="individuals")

#now repeat with subsetted datasets
#figure out which one to do first by highest total individuals
sum(CBayPonMatCh1Zon1) #16508
sum(CBayPonMatCh1Zon2) #34564
sum(CBayPonMatCh1Zon3) #18216
sum(CBayPonMatCh1Zon4) #24935
#2 first
CBaySpecAc2iCh1<-specaccum(CBayPonMatCh1Zon2, method = "rarefaction", conditioned =TRUE, 
                        gamma = "jack1",  w = NULL)
plot(CBaySpecAc2iCh1,col="#377eb8",ylim=c(0,180),xvar="individuals")
CBaySpecAc1iCh1<-specaccum(CBayPonMatCh1Zon1, method = "rarefaction", conditioned =TRUE, 
                        gamma = "jack1",  w = NULL)
plot(CBaySpecAc1iCh1,add=T,col="#e41a1c",xvar="individuals")
CBaySpecAc4iCh1<-specaccum(CBayPonMatCh1Zon4, method = "rarefaction", conditioned =TRUE, 
                        gamma = "jack1",  w = NULL)
plot(CBaySpecAc4iCh1,add=T,col="#984ea3",xvar="individuals")
CBaySpecAc3iCh1<-specaccum(CBayPonMatCh1Zon3, method = "rarefaction", conditioned =TRUE, 
                        gamma = "jack1",  w = NULL)
plot(CBaySpecAc3iCh1,add=T,col="#4daf4a",xvar="individuals")

#Estimate richness
chao1(CBayPonarMatrixCh1,taxa.row=F)
#303
chao2(CBayPonarMatrixCh1,taxa.row=F)
#312

chao1(CBayPonMatCh1Zon1,taxa.row=F)
#236
chao2(CBayPonMatCh1Zon1,taxa.row=F)
#237

chao1(CBayPonMatCh1Zon2,taxa.row=F)
#192
chao2(CBayPonMatCh1Zon2,taxa.row=F)
#222

chao1(CBayPonMatCh1Zon3,taxa.row=F)
#173
chao2(CBayPonMatCh1Zon3,taxa.row=F)
#201

chao1(CBayPonMatCh1Zon4,taxa.row=F)
#181.4
chao2(CBayPonMatCh1Zon4,taxa.row=F)
#233.1818

#Now do all again with Chao2

#Create site by taxa matrix Chao2
names(CBayspChao2)
CBayPonarMatrixCh2<-data.frame(t(create.matrix(
  CBayspChao2,
  tax.name="newtaxa",
  locality="InvEventFK",
  time.col=NULL,
  time=NULL, 
  abund=TRUE,
  abund.col="InvCount"
)))
#273 events and 253 species captured
#write csv to folder
write.csv(CBayPonarMatrixCh2,"CBayPonarMatrixChao2.csv", row.names=TRUE)

#subset into zones
CBayPonarMatrixCh2Env<-merge(CBayPonarMatrixCh2,CBayEnvrn,by=0)
names(CBayPonarMatrixCh2Env)
CBayPonarMatrixCh2Env$SedComb<-as.factor(CBayPonarMatrixCh2Env$SedComb)
levels(CBayPonarMatrixCh2Env$Zone)
CBayPonMatCh2EnvZon<-split(CBayPonarMatrixCh2Env,f=CBayPonarMatrixCh2Env$Zone)
#49 in zone 1, 76 in zone 2, 61 in zone 3, 87 in zone 4
#Create species matrices for each zone separately
CBayPonMatCh2EnvZon1<-CBayPonMatCh2EnvZon[[1]]
row.names(CBayPonMatCh2EnvZon1)<-CBayPonMatCh2EnvZon1$Row.names
CBayPonMatCh2Zon1<-CBayPonMatCh2EnvZon1[,c(0,2:254)]
CBayPonMatCh2Zon1<-CBayPonMatCh2Zon1[, which(colSums(CBayPonMatCh2Zon1) != 0)]
#167 species found in zone 1
CBayPonMatCh2EnvZon2<-CBayPonMatCh2EnvZon[[2]]
row.names(CBayPonMatCh2EnvZon2)<-CBayPonMatCh2EnvZon2$Row.names
CBayPonMatCh2Zon2<-CBayPonMatCh2EnvZon2[,c(0,2:254)]
CBayPonMatCh2Zon2<-CBayPonMatCh2Zon2[, which(colSums(CBayPonMatCh2Zon2) != 0)]
#170 species found in zone 2
CBayPonMatCh2EnvZon3<-CBayPonMatCh2EnvZon[[3]]
row.names(CBayPonMatCh2EnvZon3)<-CBayPonMatCh2EnvZon3$Row.names
CBayPonMatCh2Zon3<-CBayPonMatCh2EnvZon3[,c(0,2:254)]
CBayPonMatCh2Zon3<-CBayPonMatCh2Zon3[, which(colSums(CBayPonMatCh2Zon3) != 0)]
#150 species found in zone 3
CBayPonMatCh2EnvZon4<-CBayPonMatCh2EnvZon[[4]]
row.names(CBayPonMatCh2EnvZon4)<-CBayPonMatCh2EnvZon4$Row.names
CBayPonMatCh2Zon4<-CBayPonMatCh2EnvZon4[,c(0,2:254)]
CBayPonMatCh2Zon4<-CBayPonMatCh2Zon4[, which(colSums(CBayPonMatCh2Zon4) != 0)]
#152 species found in zone 4

#Species accumulation curve
#for all sites
CBaySpecAcCh2<-specaccum(CBayPonarMatrixCh2, method = "exact", conditioned =TRUE, 
                         gamma = "jack1",  w = NULL)
plot(CBaySpecAcCh2)
CBaypoolCh2<-poolaccum(CBayPonarMatrixCh2)
plot(CBaypoolCh2)

#now repeat with subsetted datasets
#do zone 4 first because it has the most sites
 
#Try again using individuals rather than sites
sum(CBayPonarMatrixCh2)
#92889 individuals
CBaySpecAciCh2<-specaccum(CBayPonarMatrixCh2, method = "rarefaction", conditioned =TRUE,
                          gamma = "jack1",  w = NULL)
plot(CBaySpecAciCh2, xvar="individuals")

#now repeat with subsetted datasets
#figure out which one to do first by highest total individuals
sum(CBayPonMatCh2Zon1) #16430
sum(CBayPonMatCh2Zon2) #33397
sum(CBayPonMatCh2Zon3) #18212
sum(CBayPonMatCh2Zon4) #24850
#2 first
CBaySpecAc2iCh2<-specaccum(CBayPonMatCh2Zon2, method = "rarefaction", conditioned =TRUE, 
                           gamma = "jack1",  w = NULL)
plot(CBaySpecAc2iCh2,col="#377eb8",ylim=c(0,180),xvar="individuals")
CBaySpecAc1iCh2<-specaccum(CBayPonMatCh2Zon1, method = "rarefaction", conditioned =TRUE, 
                           gamma = "jack1",  w = NULL)
plot(CBaySpecAc1iCh2,add=T,col="#e41a1c",xvar="individuals")
CBaySpecAc4iCh2<-specaccum(CBayPonMatCh2Zon4, method = "rarefaction", conditioned =TRUE, 
                           gamma = "jack1",  w = NULL)
plot(CBaySpecAc4iCh2,add=T,col="#984ea3",xvar="individuals")
CBaySpecAc3iCh2<-specaccum(CBayPonMatCh2Zon3, method = "rarefaction", conditioned =TRUE, 
                           gamma = "jack1",  w = NULL)
plot(CBaySpecAc3iCh2,add=T,col="#4daf4a",xvar="individuals")

#Estimate richness
chao1(CBayPonarMatrixCh2,taxa.row=F)
#303
chao2(CBayPonarMatrixCh2,taxa.row=F)
#315

chao1(CBayPonMatCh2Zon1,taxa.row=F)
#235
chao2(CBayPonMatCh2Zon1,taxa.row=F)
#229

chao1(CBayPonMatCh2Zon2,taxa.row=F)
#191
chao2(CBayPonMatCh2Zon2,taxa.row=F)
#240

chao1(CBayPonMatCh2Zon3,taxa.row=F)
#166
chao2(CBayPonMatCh2Zon3,taxa.row=F)
#196

chao1(CBayPonMatCh2Zon4,taxa.row=F)
#193
chao2(CBayPonMatCh2Zon4,taxa.row=F)
#245

#examine spatial differences in richness

#Triangle plots
#Threshold for depth
median(CBayPonarMatrixCh1Env$HabDepth)
#2.6 - looks similar to spline on generalised dissimilarity
#1.05 and 8.25 on random forest

#threshold for distance to marina
median(CBayPonarMatrixCh1Env$DistMarina4)
#5691

#For geographic distance - use zone as a bins

#Create triangle plots with % in same zone on x axis and % deep on y axis
# deep is greater than 2.6
#33 shallow 16 deep in zone 1
#42 shallow 34 deep in zone 2
#34 shallow 27 deep in zone 3
#28 shallow 59 deep in zone 4
#iterate on n=14, greatest efficiencies in sampling occur in the early portion of the species-effort curve (Hoffman et al. 2011)
#10 random draws

#do this in a simplified way first.
#don't a priori specify depth, just do zone
#iterate the same code 4 times, once for each zone

#Add distances to the matrix/environmental dataset
CBayPonarMatrixCh1EnvDist<-merge(CBayPonarMatrixCh1Env,CBayDist,By="InvEventFK")

##################################
#Only Run if need randomization analysis
#############################

#subset to deep sites only
CBayPonMatCh1EnvZonDist234<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone1")
CBayPonMatCh1EnvZonDist134<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone2")
CBayPonMatCh1EnvZonDist124<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone3")
CBayPonMatCh1EnvZonDist123<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone4")
CBayPonMatCh1EnvZonDist234d<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone1" & HabDepth>2.6)
CBayPonMatCh1EnvZonDist134d<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone2" & HabDepth>2.6)
CBayPonMatCh1EnvZonDist124d<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone3" & HabDepth>2.6)
CBayPonMatCh1EnvZonDist123d<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone4" & HabDepth>2.6)
CBayPonMatCh1EnvZonDist234s<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone1" & HabDepth<=2.6)
CBayPonMatCh1EnvZonDist134s<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone2" & HabDepth<=2.6)
CBayPonMatCh1EnvZonDist124s<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone3" & HabDepth<=2.6)
CBayPonMatCh1EnvZonDist123s<-subset(CBayPonarMatrixCh1EnvDist,Zone!="Zone4" & HabDepth<=2.6)
CBayPonMatCh1EnvZonDist1d<-subset(CBayPonarMatrixCh1EnvDist,Zone=="Zone1" & HabDepth>2.6)
CBayPonMatCh1EnvZonDist2d<-subset(CBayPonarMatrixCh1EnvDist,Zone=="Zone2" & HabDepth>2.6)
CBayPonMatCh1EnvZonDist3d<-subset(CBayPonarMatrixCh1EnvDist,Zone=="Zone3" & HabDepth>2.6)
CBayPonMatCh1EnvZonDist4d<-subset(CBayPonarMatrixCh1EnvDist,Zone=="Zone4" & HabDepth>2.6)
CBayPonMatCh1EnvZonDist1s<-subset(CBayPonarMatrixCh1EnvDist,Zone=="Zone1" & HabDepth<=2.6)
CBayPonMatCh1EnvZonDist2s<-subset(CBayPonarMatrixCh1EnvDist,Zone=="Zone2" & HabDepth<=2.6)
CBayPonMatCh1EnvZonDist3s<-subset(CBayPonarMatrixCh1EnvDist,Zone=="Zone3" & HabDepth<=2.6)
CBayPonMatCh1EnvZonDist4s<-subset(CBayPonarMatrixCh1EnvDist,Zone=="Zone4" & HabDepth<=2.6)


dfrasubsamples<-expand.grid(n.deepz=seq(0,14,1), n.shallowz=seq(0,14,1),
                            ndeepoz=seq(0,14,1), nshallowoz=seq(0,14,1))
dfrasubsamples$sum<-rowSums(dfrasubsamples)
dfrasubsamples$zonesum<-dfrasubsamples$n.deepz+dfrasubsamples$n.shallowz

#zone 1
dfrasubsamples1<-subset(dfrasubsamples, sum==14 & zonesum>=3)
dfrasubsamples1$sum<-NULL
dfrasubsamples1$zonesum<-NULL
dfrasubsamples1$result<-0
dfrasubsamples1$sitenames<-0
vec1<-c(1:12)

#RUN LOOPS 
#This code takes a while to run 
for (i in 1:nrow(dfrasubsamples1))
{
  for (j in 1:length(vec1))
  {
    CBayPonMatCh1EnvZonDist1d.sample<-
      CBayPonMatCh1EnvZonDist1d[sample(1:nrow(CBayPonMatCh1EnvZonDist1d), 
                                                               dfrasubsamples1[i,1], 
                                   replace=FALSE),]
    CBayPonMatCh1EnvZonDist1s.sample<-
      CBayPonMatCh1EnvZonDist1s[sample(1:nrow(CBayPonMatCh1EnvZonDist1s), dfrasubsamples1[i,2],
                                   replace=FALSE),]
    CBayPonMatCh1EnvZonDist234d.sample<-
      CBayPonMatCh1EnvZonDist234d[sample(1:nrow(CBayPonMatCh1EnvZonDist234d), 
                                     dfrasubsamples1[i,3], replace=FALSE),]
    
    CBayPonMatCh1EnvZonDist234s.sample<-
      CBayPonMatCh1EnvZonDist234s[sample(1:nrow(CBayPonMatCh1EnvZonDist234s), 
                                     dfrasubsamples1[i,4], replace=FALSE),]
    sample1<-rbind(CBayPonMatCh1EnvZonDist1d.sample, CBayPonMatCh1EnvZonDist1s.sample, 
                  CBayPonMatCh1EnvZonDist234d.sample, CBayPonMatCh1EnvZonDist234s.sample)
    inv1<-chao1(sample1[,1:253],taxa.row=F)
    vec1[j]<-sum(inv1)
  }
  dfrasubsamples1[i,5]<-mean(vec1)

}
#zone 2
dfrasubsamples2<-subset(dfrasubsamples, sum==14 & zonesum>=3)
dfrasubsamples2$sum<-NULL
dfrasubsamples2$zonesum<-NULL
dfrasubsamples2$result<-0
vec2<-c(1:19)

#RUN LOOPS 
#This code takes a while to run 
for (i in 1:nrow(dfrasubsamples2))
{
  for (j in 1:length(vec2))
  {
    CBayPonMatCh1EnvZon2d.sample<-
      CBayPonMatCh1EnvZon2d[sample(1:nrow(CBayPonMatCh1EnvZon2d), 
                                   dfrasubsamples2[i,1], 
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon2s.sample<-
      CBayPonMatCh1EnvZon2s[sample(1:nrow(CBayPonMatCh1EnvZon2s), dfrasubsamples2[i,2],
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon134d.sample<-
      CBayPonMatCh1EnvZon134d[sample(1:nrow(CBayPonMatCh1EnvZon134d), 
                                     dfrasubsamples2[i,3], replace=FALSE),]
    
    CBayPonMatCh1EnvZon134s.sample<-
      CBayPonMatCh1EnvZon134s[sample(1:nrow(CBayPonMatCh1EnvZon134s), 
                                     dfrasubsamples2[i,4], replace=FALSE),]
    sample2<-rbind(CBayPonMatCh1EnvZon2d.sample, CBayPonMatCh1EnvZon2s.sample, 
                  CBayPonMatCh1EnvZon134d.sample, CBayPonMatCh1EnvZon134s.sample)
    inv2<-chao1(sample2[,1:253],taxa.row=F)
    vec2[j]<-sum(inv2)
  }
  dfrasubsamples2[i,5]<-mean(vec2)
}
#zone 3
dfrasubsamples3<-subset(dfrasubsamples, sum==14 & zonesum>=3)
dfrasubsamples3$sum<-NULL
dfrasubsamples3$zonesum<-NULL
dfrasubsamples3$result<-0
vec3<-c(1:15)
#RUN LOOPS 
#This code takes a while to run 
for (i in 1:nrow(dfrasubsamples3))
{
  for (j in 1:length(vec3))
  {
    CBayPonMatCh1EnvZon3d.sample<-
      CBayPonMatCh1EnvZon3d[sample(1:nrow(CBayPonMatCh1EnvZon3d), 
                                   dfrasubsamples3[i,1], 
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon3s.sample<-
      CBayPonMatCh1EnvZon3s[sample(1:nrow(CBayPonMatCh1EnvZon3s), dfrasubsamples3[i,2],
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon124d.sample<-
      CBayPonMatCh1EnvZon124d[sample(1:nrow(CBayPonMatCh1EnvZon124d), 
                                     dfrasubsamples3[i,3], replace=FALSE),]
    CBayPonMatCh1EnvZon124s.sample<-
      CBayPonMatCh1EnvZon124s[sample(1:nrow(CBayPonMatCh1EnvZon124s), 
                                     dfrasubsamples3[i,4], replace=FALSE),]
    sample3<-rbind(CBayPonMatCh1EnvZon3d.sample, CBayPonMatCh1EnvZon3s.sample, 
                  CBayPonMatCh1EnvZon124d.sample, CBayPonMatCh1EnvZon124s.sample)
    inv3<-chao1(sample3[,1:253],taxa.row=F)
    vec3[j]<-sum(inv3)
  }
  dfrasubsamples3[i,5]<-mean(vec3)
}
#zone 4
dfrasubsamples4<-subset(dfrasubsamples, sum==14 & zonesum>=3)
dfrasubsamples4$sum<-NULL
dfrasubsamples4$zonesum<-NULL
dfrasubsamples4$result<-0
vec4<-c(1:22)
#RUN LOOPS 
#This code takes a while to run 
for (i in 1:nrow(dfrasubsamples4))
{
  for (j in 1:length(vec4))
  {
    CBayPonMatCh1EnvZon4d.sample<-
      CBayPonMatCh1EnvZon4d[sample(1:nrow(CBayPonMatCh1EnvZon4d), 
                                   dfrasubsamples4[i,1], 
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon4s.sample<-
      CBayPonMatCh1EnvZon4s[sample(1:nrow(CBayPonMatCh1EnvZon4s), dfrasubsamples4[i,2],
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon123d.sample<-
      CBayPonMatCh1EnvZon123d[sample(1:nrow(CBayPonMatCh1EnvZon123d), 
                                     dfrasubsamples4[i,3], replace=FALSE),]
    CBayPonMatCh1EnvZon123s.sample<-
      CBayPonMatCh1EnvZon123s[sample(1:nrow(CBayPonMatCh1EnvZon123s), 
                                     dfrasubsamples4[i,4], replace=FALSE),]
    sample4<-rbind(CBayPonMatCh1EnvZon4d.sample, CBayPonMatCh1EnvZon4s.sample, 
                   CBayPonMatCh1EnvZon123d.sample, CBayPonMatCh1EnvZon123s.sample)
    inv4<-chao1(sample4[,1:253],taxa.row=F)
    vec4[j]<-sum(inv4)
  }
  dfrasubsamples4[i,5]<-mean(vec4)
}
#combine 4 zones
dfrasubsamples<-rbind(dfrasubsamples1,dfrasubsamples2,dfrasubsamples3,dfrasubsamples4)
dfrasubsamplesagg<-aggregate(dfrasubsamples,
                             by=list(n.deepz=dfrasubsamples$n.deepz,
                            n.shallowz=dfrasubsamples$n.shallowz,
                            ndeepoz=dfrasubsamples$ndeepoz,
                            nshallowoz=dfrasubsamples$nshallowoz),FUN="mean")
#in the aggregated each result represents mean of 68 random iterations
dfrasubsamplesagg$percshallow<-((dfrasubsamplesagg$n.shallowz+dfrasubsamplesagg$nshallowoz)/14)
dfrasubsamplesagg$perczone<-((dfrasubsamplesagg$n.shallowz+dfrasubsamplesagg$n.deepz)/14)
dfrasubsamplesaggred<-dfrasubsamplesagg
dfrasubsamplesaggred$n.shallowz<-NULL
dfrasubsamplesaggred$n.deepz<-NULL
dfrasubsamplesaggred$nshallowoz<-NULL
dfrasubsamplesaggred$ndeepoz<-NULL

ggplot(dfrasubsamplesaggred, aes(percshallow, perczone))+
  geom_tile(aes(fill = result))+
  scale_fill_gradient(low = "#ffffd9", high = "#081d58", name="Chao1 Richness")+
  labs(x="% Shallow Stations", y="% Stations within same zone")+
  theme_bw()


#Try again without zone 4
#subset to deep sites only
CBayPonMatCh1EnvZon23<-subset(CBayPonarMatrixCh1Env,Zone!="Zone1"&Zone!="Zone4")
CBayPonMatCh1EnvZon13<-subset(CBayPonarMatrixCh1Env,Zone!="Zone2"&Zone!="Zone4")
CBayPonMatCh1EnvZon12<-subset(CBayPonarMatrixCh1Env,Zone!="Zone3"&Zone!="Zone4")
CBayPonMatCh1EnvZon23d<-subset(CBayPonarMatrixCh1Env,Zone!="Zone1" &Zone!="Zone4"& HabDepth>2.6)
CBayPonMatCh1EnvZon13d<-subset(CBayPonarMatrixCh1Env,Zone!="Zone2"  &Zone!="Zone4"& HabDepth>2.6)
CBayPonMatCh1EnvZon12d<-subset(CBayPonarMatrixCh1Env,Zone!="Zone3"  &Zone!="Zone4"& HabDepth>2.6)
CBayPonMatCh1EnvZon23s<-subset(CBayPonarMatrixCh1Env,Zone!="Zone1" &Zone!="Zone4" & HabDepth<=2.6)
CBayPonMatCh1EnvZon13s<-subset(CBayPonarMatrixCh1Env,Zone!="Zone2" &Zone!="Zone4" & HabDepth<=2.6)
CBayPonMatCh1EnvZon12s<-subset(CBayPonarMatrixCh1Env,Zone!="Zone3"  &Zone!="Zone4"& HabDepth<=2.6)

#zone 1
dfrasubsamples1no4<-subset(dfrasubsamples, sum==14 & zonesum>=3)
dfrasubsamples1no4$sum<-NULL
dfrasubsamples1no4$zonesum<-NULL
dfrasubsamples1no4$result<-0
vec1no4<-c(1:12)

#RUN LOOPS 
#This code takes a while to run 
for (i in 1:nrow(dfrasubsamples1no4))
{
  for (j in 1:length(vec1no4))
  {
    CBayPonMatCh1EnvZon1d.sampleno4<-
      CBayPonMatCh1EnvZon1d[sample(1:nrow(CBayPonMatCh1EnvZon1d), 
                                   dfrasubsamples1no4[i,1], 
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon1s.sampleno4<-
      CBayPonMatCh1EnvZon1s[sample(1:nrow(CBayPonMatCh1EnvZon1s), dfrasubsamples1no4[i,2],
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon23d.sampleno4<-
      CBayPonMatCh1EnvZon23d[sample(1:nrow(CBayPonMatCh1EnvZon23d), 
                                     dfrasubsamples1no4[i,3], replace=FALSE),]
    
    CBayPonMatCh1EnvZon23s.sampleno4<-
      CBayPonMatCh1EnvZon23s[sample(1:nrow(CBayPonMatCh1EnvZon23s), 
                                     dfrasubsamples1no4[i,4], replace=FALSE),]
    sample1no4<-rbind(CBayPonMatCh1EnvZon1d.sampleno4, CBayPonMatCh1EnvZon1s.sampleno4, 
                   CBayPonMatCh1EnvZon23d.sampleno4, CBayPonMatCh1EnvZon23s.sampleno4)
    inv1no4<-chao1(sample1no4[,1:253],taxa.row=F)
    vec1no4[j]<-sum(inv1no4)
  }
  dfrasubsamples1no4[i,5]<-mean(vec1no4)
}
#zone 2
dfrasubsamples2no4<-subset(dfrasubsamples, sum==14 & zonesum>=3)
dfrasubsamples2no4$sum<-NULL
dfrasubsamples2no4$zonesum<-NULL
dfrasubsamples2no4$result<-0
vec2no4<-c(1:19)

#RUN LOOPS 
#This code takes a while to run 
for (i in 1:nrow(dfrasubsamples2no4))
{
  for (j in 1:length(vec2no4))
  {
    CBayPonMatCh1EnvZon2d.sampleno4<-
      CBayPonMatCh1EnvZon2d[sample(1:nrow(CBayPonMatCh1EnvZon2d), 
                                   dfrasubsamples2no4[i,1], 
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon2s.sampleno4<-
      CBayPonMatCh1EnvZon2s[sample(1:nrow(CBayPonMatCh1EnvZon2s), dfrasubsamples2no4[i,2],
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon13d.sampleno4<-
      CBayPonMatCh1EnvZon13d[sample(1:nrow(CBayPonMatCh1EnvZon13d), 
                                     dfrasubsamples2no4[i,3], replace=FALSE),]
    
    CBayPonMatCh1EnvZon13s.sampleno4<-
      CBayPonMatCh1EnvZon13s[sample(1:nrow(CBayPonMatCh1EnvZon13s), 
                                     dfrasubsamples2no4[i,4], replace=FALSE),]
    sample2no4<-rbind(CBayPonMatCh1EnvZon2d.sampleno4, CBayPonMatCh1EnvZon2s.sampleno4, 
                   CBayPonMatCh1EnvZon13d.sampleno4, CBayPonMatCh1EnvZon13s.sampleno4)
    inv2no4<-chao1(sample2no4[,1:253],taxa.row=F)
    vec2no4[j]<-sum(inv2no4)
  }
  dfrasubsamples2no4[i,5]<-mean(vec2no4)
}
#zone 3
dfrasubsamples3no4<-subset(dfrasubsamples, sum==14 & zonesum>=3)
dfrasubsamples3no4$sum<-NULL
dfrasubsamples3no4$zonesum<-NULL
dfrasubsamples3no4$result<-0
vec3no4<-c(1:15)
#RUN LOOPS 
#This code takes a while to run 
for (i in 1:nrow(dfrasubsamples3no4))
{
  for (j in 1:length(vec3no4))
  {
    CBayPonMatCh1EnvZon3d.sampleno4<-
      CBayPonMatCh1EnvZon3d[sample(1:nrow(CBayPonMatCh1EnvZon3d), 
                                   dfrasubsamples3no4[i,1], 
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon3s.sampleno4<-
      CBayPonMatCh1EnvZon3s[sample(1:nrow(CBayPonMatCh1EnvZon3s), dfrasubsamples3no4[i,2],
                                   replace=FALSE),]
    CBayPonMatCh1EnvZon12d.sampleno4<-
      CBayPonMatCh1EnvZon12d[sample(1:nrow(CBayPonMatCh1EnvZon12d), 
                                     dfrasubsamples3no4[i,3], replace=FALSE),]
    CBayPonMatCh1EnvZon12s.sampleno4<-
      CBayPonMatCh1EnvZon12s[sample(1:nrow(CBayPonMatCh1EnvZon12s), 
                                     dfrasubsamples3no4[i,4], replace=FALSE),]
    sample3no4<-rbind(CBayPonMatCh1EnvZon3d.sampleno4, CBayPonMatCh1EnvZon3s.sampleno4, 
                   CBayPonMatCh1EnvZon12d.sampleno4, CBayPonMatCh1EnvZon12s.sampleno4)
    inv3no4<-chao1(sample3no4[,1:253],taxa.row=F)
    vec3no4[j]<-sum(inv3no4)
  }
  dfrasubsamples3no4[i,5]<-mean(vec3no4)
}
#combine 4 zones
dfrasubsamplesno4<-rbind(dfrasubsamples1no4,dfrasubsamples2no4,dfrasubsamples3no4)
dfrasubsamplesaggno4<-aggregate(dfrasubsamplesno4,
                             by=list(n.deepz=dfrasubsamplesno4$n.deepz,
                                     n.shallowz=dfrasubsamplesno4$n.shallowz,
                                     ndeepoz=dfrasubsamplesno4$ndeepoz,
                                     nshallowoz=dfrasubsamplesno4$nshallowoz),FUN="mean")
#in the aggregated each result represents mean of 68 random iterations
dfrasubsamplesaggno4$percshallow<-((dfrasubsamplesaggno4$n.shallowz+
                                      dfrasubsamplesaggno4$nshallowoz)/14)
dfrasubsamplesaggno4$perczone<-((dfrasubsamplesaggno4$n.shallowz+dfrasubsamplesaggno4$n.deepz)/14)
dfrasubsamplesaggredno4<-dfrasubsamplesaggno4
dfrasubsamplesaggredno4$n.shallowz<-NULL
dfrasubsamplesaggredno4$n.deepz<-NULL
dfrasubsamplesaggredno4$nshallowoz<-NULL
dfrasubsamplesaggredno4$ndeepoz<-NULL

ggplot(dfrasubsamplesaggredno4, aes(percshallow, perczone))+
  geom_tile(aes(fill = result))+
  scale_fill_gradient(low = "#ffffd9", high = "#081d58", name="Chao1 Richness")+
  labs(x="% Shallow Stations", y="% Stations within same zone")+
  theme_bw()

##################################
#####Start here for directional species accumulation
######################################

#Singletons (Ch1) first
#subset CBayDist so same sites as CBayPonar Matrix
uneeded<-which(!rownames(CBayDist) %in% rownames(CBayPonarMatrixCh1))    
CBayDist<-CBayDist[-uneeded,-uneeded]
CBaybetas <- directionalSAC(CBayPonarMatrixCh1, CBayDist)

##Environmental gradient
betaPC1<- directionalSAC(CBayPonarMatrixCh1, CBayPonarMatrixCh1EnvPC$PC1)

#Do again with PC2
betaPC2<- directionalSAC(CBayPonarMatrixCh1, CBayPonarMatrixCh1EnvPC$PC2)

#plot distance and PC1 in the same plot
#create merged dataframe
CBaybetas$gradient<-rep("Distance",273)
CBaybetas$NumberofSites<-as.numeric(row.names(CBaybetas))
betaPC1$gradient<-rep("PC1",273)
betaPC1$NumberofSites<-as.numeric(row.names(CBaybetas))
betaPC2$gradient<-rep("PC2",273)
betaPC2$NumberofSites<-as.numeric(row.names(CBaybetas))

#Create null dataframe
CBaybetanull<-CBaybetas
CBaybetanull$N_SCR<-CBaybetanull$N_Exact
CBaybetanull$gradient<-rep("Random",273)
CBaybetanull$Beta_Autocor<-rep("NA",273)
CBaybetanull$N_Exact<-NULL
CBaybetas$N_Exact<-NULL
betaPC1$N_Exact<-NULL
betaPC2$N_Exact<-NULL

betadistPC12<-rbind(CBaybetas,betaPC1,betaPC2,CBaybetanull)
betadistPC12$Beta_Autocor<-as.numeric(betadistPC12$Beta_Autocor)
ggplot(betadistPC12, aes(x=NumberofSites, y=N_SCR,colour=gradient)) +
  geom_point(size=2) +
  xlab("Number of sites")+
  ylab("Species richness")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.text=element_text(size=12),legend.title=element_blank())
#PC1 accumulates faster than random, but PC2 is essentially the same as random

ggplot(betadistPC12, aes(x=NumberofSites, y=Beta_Autocor,colour=gradient)) +
  geom_point(size=2) +
  xlab("Number of sites")+
  ylab("Autocorrelation")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.text=element_text(size=12),legend.title=element_blank())

#####################################
#start here for random forest for richness
#########################################

#create actual richness column
CBayPonarMatrixCh1EnvPC$Richness<-specnumber(CBayPonarMatrixCh1EnvPC[,2:254])
Ch1rf <-randomForest(Richness~Organic+Gravel+Silt+HabDepth+
                       VegCovScore,data=CBayPonarMatrixCh1EnvPC, importance=TRUE,ntree=1000)
Ch1imp = importance(Ch1rf, type=1)
varImpPlot(Ch1rf, type=1)
reprtree:::plot.getTree(Ch1rf)
#veg cover and depth most important
#since organic is important, try again with organic split dataset
CBayPonarMatrixCh1Envorg$Richness<-specnumber(CBayPonarMatrixCh1Envorg[,2:254])
Ch1rforg<-randomForest(Richness~Organic+Gravel+Silt+HabDepth+
                         VegCovScore,data=CBayPonarMatrixCh1Envorg, importance=TRUE,ntree=1000)
Ch1imporg= importance(Ch1rforg, type=1)
varImpPlot(Ch1rforg, type=1)
reprtree:::plot.getTree(Ch1rforg)

Ch1rftree<-randomForest::getTree(Ch1rf)

############################
#community differences
########################################
adonis2(CBayPonarMatrixCh1 ~ Organic+Gravel+Silt+HabDepth+
                                VegCovScore+DistMarina4, 
        data=CBayPonarMatrixCh1Env,
        method="bray", permutations = 999)

CBayPonarMatrixCh1Env$VegCovScorecat<-as.factor(CBayPonarMatrixCh1Env$VegCovScore)
CBayPonarMatrixCh1Env$HabDepthcat<- cut(CBayPonarMatrixCh1Env$HabDepth,
                       breaks=c(0, 5, 10,40),
                       labels=c('shallow', 'medium', 'deep'))
summary(CBayPonarMatrixCh1Env$DistMarina4)
CBayPonarMatrixCh1Env$DistMarina4cat<- cut(CBayPonarMatrixCh1Env$DistMarina4,
                                        breaks=c(0, 3000, 9900,16000),
                                        labels=c('Close', 'Medium', 'Far'))

#nmds all samples
CBay_NMDSCh1<-metaMDS(CBayPonarMatrixCh1, distance = "bray", k = 3, 
                      maxit=1000, trymax = 300, wascores = FALSE,
                      autotransform = FALSE, trace = 2, noshare = FALSE)
#stress = 0.18
write.csv(CBay_NMDSCh1[["points"]],"NMDSCBayPonar.csv")

#Stressplot 
stressplot(CBay_NMDSCh1)

#NMDS plot for zone
ordiplot(CBay_NMDSCh1, type="n")
with(CBay_NMDSCh1, points(CBay_NMDSCh1, display="sites", col=Zone_col_vec[CBayPonarMatrixCh1Env$Zone], pch=19))
with(CBay_NMDSCh1, legend("topleft", legend=levels(CBayPonarMatrixCh1Env$Zone), bty="n", col=Zone_col_vec, pch=19, pt.bg=Zone_col_vec))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$Zone, kind="se", conf=0.95, lwd=2, col="#e41a1c", show.groups = "Zone1"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$Zone, kind="se", conf=0.95, lwd=2, col="#377eb8", show.groups = "Zone2"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$Zone, kind="se", conf=0.95, lwd=2, col="#4daf4a", show.groups = "Zone3"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$Zone, kind="se", conf=0.95, lwd=2, col="#984ea3", show.groups = "Zone4"))

#nmds veg
CBayPonarMatrixCh1Env$VegCovScorecat<-as.factor(CBayPonarMatrixCh1Env$VegCovScore)
ordiplot(CBay_NMDSCh1, type="n")
with(CBay_NMDSCh1, points(CBay_NMDSCh1, display="sites", col=veg_col_vec[CBayPonarMatrixCh1Env$VegCovScorecat], pch=19))
with(CBay_NMDSCh1, legend("topleft", legend=levels(CBayPonarMatrixCh1Env$VegCovScorecat), bty="n", col=veg_col_vec, pch=19, pt.bg=veg_col_vec))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$VegCovScorecat, kind="se", conf=0.95, lwd=2, col="#99d8c9", show.groups = "0"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$VegCovScorecat, kind="se", conf=0.95, lwd=2, col="#41ae76", show.groups = "1"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$VegCovScorecat, kind="se", conf=0.95, lwd=2, col="#006d2c", show.groups = "2"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$VegCovScorecat, kind="se", conf=0.95, lwd=2, col="#00441b", show.groups = "3"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$VegCovScorecat, kind="se", conf=0.95, lwd=2, col="black", show.groups = "4"))

#nmds plot for depth
ordiplot(CBay_NMDSCh1, type="n")
with(CBay_NMDSCh1, points(CBay_NMDSCh1, display="sites", col=Depth_col_vec[CBayPonarMatrixCh1Env$HabDepthcat], pch=19))
with(CBay_NMDSCh1, legend("topleft", legend=levels(CBayPonarMatrixCh1Env$HabDepthcat), bty="n", col=Depth_col_vec, pch=19, pt.bg=Depth_col_vec))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$HabDepthcat, kind="se", conf=0.95, lwd=2, col="#d0d1e6", show.groups = "shallow"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$HabDepthcat, kind="se", conf=0.95, lwd=2, col="#3690c0", show.groups = "medium"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$HabDepthcat, kind="se", conf=0.95, lwd=2, col="#023858", show.groups = "deep"))

#nmds plot for distance to marina
ordiplot(CBay_NMDSCh1, type="n")
with(CBay_NMDSCh1, points(CBay_NMDSCh1, display="sites", col=Marina4_col_vec[CBayPonarMatrixCh1Env$DistMarina4cat], pch=19))
with(CBay_NMDSCh1, legend("topleft", legend=levels(CBayPonarMatrixCh1Env$DistMarina4cat), bty="n", col=Marina4_col_vec, pch=19, pt.bg=Marina4_col_vec))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$DistMarina4cat, kind="se", conf=0.95, lwd=2, col="#fcbba1", show.groups = "Close"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$DistMarina4cat, kind="se", conf=0.95, lwd=2, col="#cb181d", show.groups = "Medium"))
with(CBay_NMDSCh1, ordiellipse(CBay_NMDSCh1, CBayPonarMatrixCh1Env$DistMarina4cat, kind="se", conf=0.95, lwd=2, col="#67000d", show.groups = "Far"))

#Create bray curtis distance matrix 
bcCBayPonarCh1<-as.matrix(vegdist(CBayPonarMatrixCh1, "bray"))
write.csv(bcCBayPonar,"braycurtisCBayPonarCh1.csv")

#indicator species

#indicator species analysis for zone
CBay_Com_z_indicch1<-data.frame(signassoc(CBayPonarMatrixCh1, cluster=CBayPonarMatrixCh1Env$ZoneCode,  mode=0, alternative = "two.sided",control = how(nperm=999)))
CBay_Com_z_indicch1_sig<-subset(CBay_Com_z_indicch1, psidak<=0.05)
#67 indicator species for vegetation cover

######################################
#Start here for Generalised dissimilarity modelling
###########################################

#format dataset for gdm
#Convert distances to geographic coordinates
mdsCBayDist<-data.frame(cmdscale(CBayDist, k = 2))
mdsCBayDist$InvEventFK<-rownames(mdsCBayDist)

#Create environmental only dataset
CBayPonarMatrixCh1EnvCoord<-merge(CBayPonarMatrixCh1Env,mdsCBayDist, by="InvEventFK")
names(CBayPonarMatrixCh1EnvCoord)
CBayPonarEnvCh1<-CBayPonarMatrixCh1EnvCoord[,c(1,256:258,271,276,294,296:297)]
write.csv(CBayPonarEnvCh1,"GDMCBay.csv", row.names=TRUE)
CBayPonarCh1SpSite<-CBayPonarMatrixCh1EnvCoord[,c(1:254)]
CBayPonargdmformat<-formatsitepair(CBayPonarCh1SpSite, bioFormat=1, dist="bray", 
                                   abundance=TRUE, siteColumn="InvEventFK", sppColumn=NULL,
                                   abundColumn=NULL, sppFilter=0, predData=CBayPonarEnvCh1,
                                   XColumn="X1",YColumn="X2",distPreds = NULL,
                                   weightType="equal", custWeights=NULL,
                                   sampleSites=1, verbose=FALSE)
#default of 3 spline per predictor
CBayCh1gdm.1 <- gdm(CBayPonargdmformat, geo=T)
#validate model
CBayCh1gdm.1cv<-gdm.crossvalidation(CBayPonargdmformat, train.proportion=0.9, n.crossvalid.tests=1,
                                    geo=T, splines=NULL, knots=NULL)

summary(CBayCh1gdm.1)
#[1] NULL Deviance:  4918.469
#[1] GDM Deviance:  3993.602
#[1] Percent Deviance Explained:  18.804
#pretty low deviance explained - should expect between 20-50, but also have large number of points
#Intercept:  1.1

plot(CBayCh1gdm.1)

##This line of code takes a long time to run
##########
CBayCh1modTest <- gdm.varImp(CBayPonargdmformat, geo=T, nPerm=50, parallel=T, cores=10)
########
#perc deviance explained by the full model 18.804
#predictors used: geographic, silt, habdepth, vegcovscore, distmarina4
#p values:             All predictors
#Geographic            0.00
#Silt                  0.48
#HabDepth              0.00
#VegCovScore           0.00
#DistMarina4           0.02

#Predictor importance
#All predictors
#Geographic           4.643
#Silt                 0.119
#HabDepth            69.325
#VegCovScore          8.011
#DistMarina4          3.820

#Try again with only geographic, habdepth, vegcovscore, and marina as predictors
names(CBayPonarEnvCh1)
CBayPonarEnvCh1a<-CBayPonarEnvCh1[,c(1,5:9)]
#convert to raster
x <- raster(xmn=-15762.63, xmx=30424.78, ymn=-7019.836, ymx=12412.76, res=100, 
            crs="+proj=utm +zone=15 +datum=NAD84")

us_fire <- rasterize(CBayPonarEnvCh1a[, c('X1', 'X2')], x, CBayPonarEnvCh1a[, c('HabDepth','VegCovScore','DistMarina4')], fun=mean)
envRast <- stack(us_fire)
#Back to original code
CBayPonargdmformata<-formatsitepair(CBayPonarCh1SpSite, bioFormat=1, dist="bray", 
                                   abundance=TRUE, siteColumn="InvEventFK", sppColumn=NULL,
                                   abundColumn=NULL, sppFilter=0, predData=CBayPonarEnvCh1a,
                                   XColumn="X1",YColumn="X2",distPreds = NULL,
                                   weightType="equal", custWeights=NULL,
                                   sampleSites=1, verbose=FALSE)
#default of 3 spline per predictor
CBayCh1agdm.1 <- gdm(CBayPonargdmformata, geo=T)
#validate model
CBayCh1agdm.1cv<-gdm.crossvalidation(CBayPonargdmformata, train.proportion=0.9, n.crossvalid.tests=1,
                                     geo=T, splines=NULL, knots=NULL)

summary(CBayCh1agdm.1)
#[1] NULL Deviance:  4918.469
#[1] GDM Deviance:  3994.706
#[1] Percent Deviance Explained:  18.782
#pretty low deviance explained - should expect between 20-50, but also have large number of points
#Intercept:  1.111

plot(CBayCh1agdm.1)
#######
##This line of code takes a long time to run
CBayCh1amodTest <- gdm.varImp(CBayPonargdmformata, geo=T, nPerm=50)
#all p values <0.00
#predictor importance All predictors
#Geographic           4.655
#HabDepth            69.637
#VegCovScore          8.175
#DistMarina4          3.964

#make predictions
CBayCh1agdm.1.pred <- predict(object=CBayCh1agdm.1,
                      data=CBayPonargdmformata)
plot(CBayPonargdmformata$distance,
     CBayCh1agdm.1.pred,
     xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity",
     xlim=c(0,1),
     ylim=c(0,1),
     pch=20,
     col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))
transRasts<- gdm.transform(model=CBayCh1agdm.1,
                            data=envRast)
plot(transRasts$HabDepth)
plot(envRast)
#Do again with Lat Long
#Upload interpolated dataset
CBayInt<-read.csv("ezD9214116_CBO122722_TableToExcel.csv",header=T)
names(CBayInt)
CBayPonarEnvCh1aLLint<-CBayInt[,c(272,273,277,278,280)]
#convert to raster
x <- raster(xmn=657278.8, xmx=676378.8, ymn=5161002, ymx=5202402, res=300, 
            crs="+proj=utm +zone=15 +datum=NAD84")

us_fire <- rasterize(CBayPonarEnvCh1aLLint[, c('X1', 'X2')], x, CBayPonarEnvCh1aLLint[, c('HabDepth','VegCovScore',"DistMarina4")], fun=mean)
envRast <- stack(us_fire)
plot(envRast)
transRastsint<- gdm.transform(model=CBayCh1agdm.1,
                           data=envRast)
plot(transRastsint)
writeRaster(transRastsint, filename = "CBaygdmint.grd")
##########################
#Other code here
##########################

#GRAPH OF RANDOMIZATION RESULTS
df12$prop.trawl<-df12$n.trawl/20
df12$prop.fyke<-df12$n.fyke/20
df12$prop.efish<-df12$n.efish/20
Fig.6.A.Location<-ggplot(df12, aes(prop.trawl, prop.efish))+
  geom_tile(aes(fill = result), color='white')+
  scale_fill_gradient2(name="", low ="#0072B2", mid='white', high ="red4", space = 'rgb', guide = "colourbar", midpoint = (max(df12$result)+min(df12$result))/2, breaks=c(min(df12$result),(max(df12$result)+min(df12$result))/2, max(df12$result)), labels=c(round(min(df12$result)), round((max(df12$result)+min(df12$result))/2), round(max(df12$result))))+ 
  annotate("text", label="Total", x=ifelse(round(max(df12$result))>=10,0.93,0.95), y=1.0, hjust="right", colour="black")+ 
  annotate("rect", xmin=0.3, xmax=0.4, ymin=0.3, ymax=0.4, alpha=0.5)+
  theme_bw()+
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.title=element_blank())+
  scale_x_continuous(breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  scale_y_continuous(breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  xlab("")+
  ylab("")
png(file="c:/Myers/Projects/Aquatic Invasive Species/2016/Annual Report/Fig.6.A.Location.png",width=1700,height=1600, res=300)
plot(Fig.6.A.Location)
dev.off()
```

```{r eval=params$GEAR_COMPARISON, echo=FALSE, message=FALSE, warning=FALSE, include=FALSE}
#RANDOMIZATION ANALYSIS FOR "RARE 5%" RICHNESS
df9<-df7
df9$Total_Catch.sum.sum<-NULL  
#DETERMINE WHICH SPECIES QUALIFY AS "RARE 5%"
rare_05<-summaryBy(Presence_Absence ~ Species_Name, data=df9, FUN=c(sum))
rare_05$Richness<-rare_05$Presence_Absence.sum/nrow(df1)
rare_05<-rare_05[with(rare_05, order(Richness)),]
rare_05$Presence_Absence.sum<-NULL #Remove unnecesary column
rare_05$YesNo<-ifelse(rare_05$Richness<=0.05,"1","0")
rare_05$Richness<-NULL #Remove unnecesary column
df9<-merge(df9, rare_05, by="Species_Name")
names(df9)[names(df9)=="YesNo"] <- "Rare_05"
df9$Rare_05<-as.numeric(df9$Rare_05)
#DETERMINE HOW MANY "RARE 5%" SPECIES WERE CAUGHT AT EACH SAMPLING LOCATION
df9$dummy<-df9$Rare_05*df9$Presence_Absence
df9<-summaryBy(dummy~OptID+Species_Name, data=df9, FUN=c(sum))
names(df9)[names(df9)=="dummy.sum"]<-"Rare_05"
df9<-merge(df9, df1, by="OptID")
df9<-df9[c("OptID", "Species_Name", "Rare_05", "Gear_Code")]
df9$Gear_Code<-as.character(df9$Gear_Code)
#SUBSET SITES ACCORDING TO GEAR TYPE
df10<-data.frame(OptID=df1$OptID, Gear_Code=as.character(df1$Gear_Code))
trawl<-subset(df10, Gear_Code=="24") #BOTTOM TRAWL
fyke<-subset(df10, Gear_Code=="2") #FYKE NET
efish<-subset(df10, Gear_Code=="6") #BOAT ELECTROSHOCKER
#PREPARE BLANK DATA FRAME FOR RANDOM DRAWS USING SAME METHODS AS TREBITZ ET AL. (2009) AND HOFFMAN ET AL. (2016)
df11<-expand.grid(n.trawl=seq(0,20,2), n.fyke=seq(0,20,2), n.efish=seq(0,20,2))
df11$sum<-rowSums(df11)
df12<-subset(df11, sum==20)
df12$sum<-NULL
df12$result<-0
vec1<-c(1:100)
#RUN LOOPS
for (i in 1:nrow(df12))
{
  for (j in 1:length(vec1))
  {
    trawl.sample<-trawl[sample(1:nrow(trawl), df12[i,1], replace=FALSE),]
    fyke.sample<-fyke[sample(1:nrow(fyke), df12[i,2], replace=FALSE),]
    efish.sample<-efish[sample(1:nrow(efish), df12[i,3], replace=FALSE),]
    sample<-rbind(trawl.sample, fyke.sample, efish.sample)
    fish<-merge(sample, df9, by="OptID")
    fish.2<-summaryBy(Rare_05~Species_Name, data=fish, FUN=c(sum))
    fish.2$Rare_05.sum<-ifelse(fish.2$Rare_05.sum>=1,1,0)
    vec1[j]<-sum(fish.2$Rare_05.sum)
  }
  df12[i,4]<-mean(vec1)
}
#GRAPH OF RANDOMIZATION RESULTS
df12$prop.trawl<-df12$n.trawl/20
df12$prop.fyke<-df12$n.fyke/20
df12$prop.efish<-df12$n.efish/20
Fig.6.B.Location<-ggplot(df12, aes(prop.trawl, prop.efish))+
  geom_tile(aes(fill = result), color='white')+
  scale_fill_gradient2(name="", low ="#0072B2", mid='white', high ="red4", space = 'rgb', guide = "colourbar", midpoint = (max(df12$result)+min(df12$result))/2, breaks=c(min(df12$result),(max(df12$result)+min(df12$result))/2, max(df12$result)), labels=c(round(min(df12$result)), round((max(df12$result)+min(df12$result))/2), round(max(df12$result))))+ 
  annotate("text", label="Rare 5%", x=ifelse(round(max(df12$result))>=10,0.93,0.95), y=1.0, hjust="right", colour="black")+ 
  annotate("rect", xmin=0.3, xmax=0.4, ymin=0.3, ymax=0.4, alpha=0.5)+
  theme_bw()+
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.title=element_blank())+
  scale_x_continuous(breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  scale_y_continuous(breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  xlab("")+
  ylab("Proportion Electrofishing Samples")
png(file="c:/Myers/Projects/Aquatic Invasive Species/2016/Annual Report/Fig.6.B.Location.png",width=1700,height=1600, res=300)
plot(Fig.6.B.Location)
dev.off()
```

```{r eval=params$GEAR_COMPARISON, echo=FALSE, message=FALSE, warning=FALSE, include=FALSE}
#RANDOMIZATION ANALYSIS FOR NON-NATIVE RICHNESS
df9<-df7
df9$Total_Catch.sum.sum<-NULL
#IDENTIFY NON-NATIVE SPECIES
df9$Non_Native<-ifelse(df9$Species_Name=="Alewife" | 
                         df9$Species_Name=="Brook silverside" | 
                         df9$Species_Name=="Brown trout" |
                         df9$Species_Name=="Chinook salmon" |
                         df9$Species_Name=="Coho salmon" |
                         df9$Species_Name=="Common carp" |
                         df9$Species_Name=="Fourspine Stickleback" |
                         df9$Species_Name=="Freshwater drum" |
                         df9$Species_Name=="Pink salmon" |
                         df9$Species_Name=="Rainbow smelt" |
                         df9$Species_Name=="Rainbow trout (Steelhead)" |
                         df9$Species_Name=="Round goby" |
                         df9$Species_Name=="Ruffe" |
                         df9$Species_Name=="Splake (brook trout x lake trout)" |
                         df9$Species_Name=="Threespine stickleback" |
                         df9$Species_Name=="Tubenose goby" |
                         df9$Species_Name=="White bass" |
                         df9$Species_Name=="White perch", 1, 0)
#DETERMINE HOW MANY NON-NATIVE SPECIES WERE CAUGHT AT EACH SAMPLING LOCATION
df9$dummy<-df9$Non_Native*df9$Presence_Absence
df9<-summaryBy(dummy~OptID+Species_Name, data=df9, FUN=c(sum))
names(df9)[names(df9)=="dummy.sum"]<-"Non_Native"
df9<-merge(df9, df1, by="OptID")
df9<-df9[c("OptID", "Species_Name", "Non_Native", "Gear_Code")]
df9$Gear_Code<-as.character(df9$Gear_Code)
#SUBSET SITES ACCORDING TO GEAR TYPE
df10<-data.frame(OptID=df1$OptID, Gear_Code=as.character(df1$Gear_Code))
trawl<-subset(df10, Gear_Code=="24") #BOTTOM TRAWL
fyke<-subset(df10, Gear_Code=="2") #FYKE NET
efish<-subset(df10, Gear_Code=="6") #BOAT ELECTROSHOCKER
#PREPARE BLANK DATA FRAME FOR RANDOM DRAWS USING SAME METHODS AS TREBITZ ET AL. (2009) AND HOFFMAN ET AL. (2016)
df11<-expand.grid(n.trawl=seq(0,20,2), n.fyke=seq(0,20,2), n.efish=seq(0,20,2))
df11$sum<-rowSums(df11)
df12<-subset(df11, sum==20)
df12$sum<-NULL
df12$result<-0
vec1<-c(1:100)
#RUN LOOPS
for (i in 1:nrow(df12))
{
  for (j in 1:length(vec1))
  {
    trawl.sample<-trawl[sample(1:nrow(trawl), df12[i,1], replace=FALSE),]
    fyke.sample<-fyke[sample(1:nrow(fyke), df12[i,2], replace=FALSE),]
    efish.sample<-efish[sample(1:nrow(efish), df12[i,3], replace=FALSE),]
    sample<-rbind(trawl.sample, fyke.sample, efish.sample)
    fish<-merge(sample, df9, by="OptID")
    fish.2<-summaryBy(Non_Native~Species_Name, data=fish, FUN=c(sum))
    fish.2$Non_Native.sum<-ifelse(fish.2$Non_Native.sum>=1,1,0)
    vec1[j]<-sum(fish.2$Non_Native.sum)
  }
  df12[i,4]<-mean(vec1)
}
#GRAPH OF RANDOMIZATION RESULTS
df12$prop.trawl<-df12$n.trawl/20
df12$prop.fyke<-df12$n.fyke/20
df12$prop.efish<-df12$n.efish/20
Fig.6.C.Location<-ggplot(df12, aes(prop.trawl, prop.efish))+
  geom_tile(aes(fill = result), color='white')+
  scale_fill_gradient2(name="", low ="#0072B2", mid='white', high ="red4", space = 'rgb', guide = "colourbar", midpoint = (max(df12$result)+min(df12$result))/2, breaks=c(min(df12$result),(max(df12$result)+min(df12$result))/2, max(df12$result)), labels=c(round(min(df12$result)), round((max(df12$result)+min(df12$result))/2), round(max(df12$result))))+ 
  annotate("text", label="Non-Native", x=ifelse(round(max(df12$result))>=10,0.93,0.95), y=1.0, hjust="right", colour="black")+ 
  annotate("rect", xmin=0.3, xmax=0.4, ymin=0.3, ymax=0.4, alpha=0.5)+
  theme_bw()+
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.title=element_blank())+
  scale_x_continuous(breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  scale_y_continuous(breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  xlab("Proportion Bottom-Trawl Samples")+
  ylab("")
png(file="c:/Myers/Projects/Aquatic Invasive Species/2016/Annual Report/Fig.6.C.Location.png",width=1700,height=1600, res=300)
plot(Fig.6.C.Location)
dev.off()
```

```{r eval=params$GEAR_COMPARISON, echo=FALSE, message=FALSE, warning=FALSE, include=FALSE}
Fig.6.Location<-grid.arrange(Fig.6.A.Location, Fig.6.B.Location, Fig.6.C.Location, ncol=1)
png(file="c:/Myers/Projects/Aquatic Invasive Species/2016/Annual Report/Fig.6.Location.png",width=1700,height=2400, res=300)
plot(Fig.6.Location)
dev.off()

```

![](Fig.6.Location.png) 

*`r if(params$GEAR_COMPARISON){"Figure 6. Estimates of total, rare 5%, and non-native species richness generated from a random selection of 20 sites using different combinations "}'*
