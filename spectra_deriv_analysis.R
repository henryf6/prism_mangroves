## ******************************************************************** ##
## mangrove_deriv_analysis.R 
##
## Author: Henry Frye
## Date Created: Aug 29 2019
##
##
## Purpose:
## Spectral derivative analysis
## ******************************************************************** ##

## ******************************************************************** ##
####Set Up####
## ******************************************************************** ##
rm(list = ls())

# Set working directory
if(Sys.info()['user']=='henryfrye') setwd("/Users/henryfrye/Dropbox/Intellectual_Endeavours/UConn/Research/phd_research/chapter_2")

# Loading Packages
library(tidyverse)
library(ggthemes)
library(reshape2)
library(hsdar)


# Load Data
black <- read.csv(file = 'mangrove_data/Leaf_spectra/blackMangroves_leaves.csv')
red <- read.csv(file = 'mangrove_data/Leaf_spectra/redMangroves_leaves.csv')
white <- read.csv(file = 'mangrove_data/Leaf_spectra/whiteMangroves_leaves.csv')

#transpose data format (run once only or re-run with removing objects)
bt <- t(black)
colnames(bt) <- bt[1,]
black <- bt[-1,]

rt <- t(red)
colnames(rt) <- rt[1,]
red <- rt[-1,]

wt <- t(white)
colnames(wt) <- wt[1,]
white <- wt[-1,]

#merge all data into one dataframe
all_spp <- rbind(black,white,red)
all_spp <- as.data.frame(all_spp)
all_spp$species <- c(rep("black",20), rep("white",20), rep("red",20) )
all_spp <- all_spp %>% dplyr::select(species,'350':'2500')

## ******************************************************************** ##
####Spectral Reading####
## ******************************************************************** ##
mangrove.spec <- speclib(as.matrix(all_spp[,2:1402]),350:1750)
SI(mangrove.spec) <- all_spp$species


# meanflt.white <- noiseFiltering(white.spec, method = "mean", p = 19)
# plot(meanflt.white)
# plot(white.spec, FUN = 1, col = "red")

mangrove.spec <- noiseFiltering(mangrove.spec, method = "mean", p = 15)
mask(mangrove.spec) <- c(1325,1500)


sp_red <- subset(mangrove.spec, V1 == "red")
sp_black <- subset(mangrove.spec, V1 == "black")
sp_white <- subset(mangrove.spec, V1 == "white")

png("mangrove_figures/species_spectra.png")
plot(sp_red, col = "darkred", ylim = c(0,1), main = "Mean Species Reflectance", lwd = 2)
plot(sp_black, col = "black",new = FALSE, lwd = 2)
plot(sp_white, col = "blue",new = FALSE, lwd = 2)
dev.off()

## ******************************************************************** ##
####Band Depth Index####
## ******************************************************************** ##
#convex hull continuum removal transformation
ch_bd.r <- transformSpeclib(sp_red, method = "ch")
ch_bd.b <- transformSpeclib(sp_black, method = "ch")
ch_bd.w <- transformSpeclib(sp_white, method = "ch")
ch_bd.all <- transformSpeclib(mangrove.spec, method = "ch")


png("mangrove_figures/band_index.png")
plot(ch_bd.r, col = "darkred", main = "Band depth", lwd = 2)
plot(ch_bd.b, col = "black", main = "Band depth", new = FALSE, lwd = 2)
plot(ch_bd.w,col = "blue", main = "Band depth", new = FALSE, lwd = 2)
abline(v = 484, col = "red")
abline(v = 667, col = "red")
#abline(v = 797, col = "red")
#abline(v = 865, col = "red")
abline(v = 976, col = "red")
abline(v = 1181, col = "red")
dev.off()
#
#getcp(ch_cline,2)

plot(ch_bd.all,col = "blue", main = "Band depth")
#explore highest band depths
mean.band<- colMeans(ch_bd.all@spectra@spectra_ma)
attributes(mean.band)$names <- c(350:1325,1500:1750)
#based on band depths, 484, 667, 976 and 1181 appear to be of interest 
#run function below
which.peaks(mean.band)


## ******************************************************************** ##
####Absorption feature extraction####
## ******************************************************************** ##
featureSelection <- specfeat(ch_bd.all, c(480,780, 975, 1150))
plot(featureSelection, fnumber = 1:4)

## Calculate properties of features
featureProp <- feature_properties(featureSelection)

sp_feat <- as.data.frame(SI(featureProp))

ggplot(sp_feat, aes(y = f480_maxwl, color = V1)) + geom_boxplot()


## ******************************************************************** ##
####Find an index####
## ******************************************************************** ##

# Calculate normalized band depth index for first feature
featureSelection_bdri <- bdri(featureSelection, 3, index = "ndbi")

## Plot result
plot(featureSelection_bdri)

## ******************************************************************** ##
####Within and out group average spectral angles for entire spectrum used.####
## ******************************************************************** ##
rad2deg <- function(rad) {(rad * 180) / (pi)}
#create empty vector of average spectral angles for each species for within group variability 
sa_red <- c()
sa_black <- c()
sa_white <- c()


for(i in 1:20){
sa_red[i] <- mean(sam(sp_red,sp_red[i,]), na.rm =TRUE)
}
mean(rad2deg(sa_red)) #3.1535
sd(rad2deg(sa_red)) # 1.209211
range(rad2deg(sa_red)) #2.372152 8.028901


for(i in 1:20){
  sa_black[i] <- mean(sam(sp_black,sp_black[i,]), na.rm =TRUE)
}
mean(rad2deg(sa_black)) #5.21
sd(rad2deg(sa_black)) #  1.764257
range(rad2deg(sa_black)) # 3.689861 9.257968

for(i in 1:20){
  sa_white[i] <- mean(sam(sp_white,sp_white[i,]), na.rm =TRUE)
}
mean(rad2deg(sa_white)) #3.150645
sd(rad2deg(sa_white)) #0.6447496
range(rad2deg(sa_white)) #2.477099 4.800206

#across group spectral angle variability 
b_r_dist <- c()
r_b_dist <- c()
b_w_dist <- c()
r_w_dist <- c()

for(i in 1:20){
b_r_dist[i] <- mean(sam(sp_red,
                sp_black[i,]),na.rm=TRUE)
}
mean(rad2deg(b_r_dist)) #6.626706
sd(rad2deg(b_r_dist)) #2.211052
range(rad2deg(b_r_dist)) #4.109806 12.244118


#just to check that inverse is the same value
for(i in 1:20){
  r_b_dist[i] <- mean(sam(sp_black,
                          sp_red[i,]),na.rm=TRUE)
}
mean(b_r_dist) #0.1156578



for(i in 1:20){
  b_w_dist[i] <- mean(sam(sp_black,
                          sp_white[i,]),na.rm=TRUE)
}
mean(rad2deg(b_w_dist)) #6.243229
sd(rad2deg(b_w_dist)) #0.9866146
range(rad2deg(b_w_dist)) #4.815639 8.227854

for(i in 1:20){
  r_w_dist[i] <- mean(sam(sp_red,
                          sp_white[i,]),na.rm=TRUE)
}
mean(rad2deg(r_w_dist)) # 3.68113
sd(rad2deg(r_w_dist)) #0.7988784
range(rad2deg(r_w_dist)) #2.453275 5.374340

## ******************************************************************** ##
####Graphic of threshholding####
## ******************************************************************** ##

angle_thresh <- data.frame("species" <- c(rep("red", 20),rep("black",20),rep("white",20)),
            "angle" <- c(rad2deg(sa_red),rad2deg(sa_black),rad2deg(sa_white)))
ggplot(angle_thresh, aes(angle, fill = species)) + 
  geom_histogram(alpha = .4, position = "identity") +
  xlab('Spectral Angle (degrees)') + 
  ylab('Count') +
  theme_tufte() +
  guides(fill =  guide_legend(title="Species")) +
  theme(text = element_text(size=12, family = "Helvetica"),legend.position = c(.8,.55))

png("mangrove_figures/species_spectra.png")
