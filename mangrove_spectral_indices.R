## ******************************************************************** ##
## mangrove_spectral_indices.R 
##
## Author: Henry Frye
## Date Created: Aug 29 2019
##
##
## Purpose:
## Comparison of major spectral indices between species
## Indices from summary table 1 of Ollinger 2011 in New Phytologist
## ******************************************************************** ##

## ******************************************************************** ##
####Set Up####
## ******************************************************************** ##
rm(list = ls())

# Set working directory
if(Sys.info()['user']=='henryfrye') setwd("/Users/henryfrye/Dropbox/Intellectual_Endeavours/UConn/Research/phd_research/chapter_2")

# Loading Packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(reshape2)

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

specindices <- data.frame(species = c(rep("black",20), rep("white",20), rep("red",20) ),
                          NDVI = (all_spp$'800' - all_spp$'680') / (all_spp$'800' + all_spp$'680'),
                          PRI = (all_spp$'531' - all_spp$'570') / (all_spp$'531' + all_spp$'570'),
                          WI = (all_spp$'900'/ all_spp$'970'))


## ******************************************************************** ##
####Spectral index comparison between species####
## ******************************************************************** ##

#NDVI
NDVI <- ggplot(specindices, aes(x = species, y = NDVI)) + 
  geom_boxplot(notch = TRUE) + 
  xlab("Species") +
  ylab("Normalized Difference Vegetation Index") +
  theme_classic() +
  theme(text = element_text(size=20, family = "Helvetica")) 
NDVI
ggsave("mangrove_figures/NDVI_species.png", NDVI, dpi = 300)


#PRI
PRI <- ggplot(specindices, aes(x = species, y = PRI)) + 
  geom_boxplot(notch = TRUE) + 
  xlab("Species") +
  ylab("Photochemical Reflectance Index") +
  theme_classic() +
  theme(text = element_text(size=20, family = "Helvetica"))
PRI
ggsave("mangrove_figures/PRI_species.png", PRI, dpi = 300)

#WI
WI <- ggplot(specindices, aes(x = species, y = WI)) +
  geom_boxplot(notch = TRUE) + 
  xlab("Species") +
  ylab("Water Index") +
  theme_classic() + 
  theme(text = element_text(size=20, family = "Helvetica")) 
WI
ggsave("mangrove_figures/WI_species.png",WI, dpi = 300)

