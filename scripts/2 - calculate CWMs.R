## Calculate plot-level CWMs
## Author: Emily H
## Created: November 27, 2023
## Last edited: April 18, 2023

#install.packages("tidyverse")
#install.packages("FD")
#install.packages("lme4)
#install.packages("lmerTest")
#install.packages("car")
#install.packages("nnet")

library(FD)
library(lme4)
library(lmerTest)
library(car)
library(nnet)
library(tidyverse)

setwd("~/Documents/GitHub/Analyses/RRI Analyses 2024")

#### import data ####
## import July 2021 data
kin21 <- read.csv("data/July 2021 % cover.csv")
kin21[kin21 == "t"] <- '1'
kin21 <- kin21 %>%
  mutate_at(c(6:76), as.numeric) %>%
  dplyr::select(!c(Date.Collected, Block, Waypoint, Treatment, Bare.Ground, Litter, moss, Cow.pat, X)) %>%
  rename(Hel.hoo2 = Heli.hoo) %>%
  mutate(Hel.hoo1 = Hel.hoo + Hel.hoo2,
         Ery.che1 = Ery.che + rough.mustard) %>% #collapse double observations
  dplyr::select(!c(Hel.hoo, Hel.hoo2, fuzz.auricle.grass, slender.wheat.grass, Ery.che, rough.mustard)) %>%
  dplyr::select(!(starts_with("un"))) %>%
  rename(Ago.gla = Tradub.like,
         Hel.hoo = Hel.hoo1,
         Ery.che = Ery.che1,
         Pas.smi = Pas.mi) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  column_to_rownames("Plot") %>%
  dplyr::select(.,order(colnames(.))) #put columns in alphabetical order
str(kin21)

## import treatment data 
trts <- read.csv("data/Plot treatments.csv") %>%
  dplyr::select(!X)
str(trts)

## import 2023 richness data
no.spp23 <- read_rds("output/richness23.rds")

## import trait database spp average traits
traits <- read_rds("data/Cahill trait database species averages.rds") %>% # trait data: species are rows and traits are columns
  rownames_to_column(var = "Species")
traits$Species <- paste0(substr(traits$Species, start = 1, stop = 3),
                         ".",
                         substr(traits$Species, start = 4, stop = 6)) #add "." into spp codes
str(traits)

#### Calculate CWMs ####
traits <- traits %>%
  dplyr::select(Species, Max.Height, SLA.Leaf, SRL, RTD) %>% #select traits I want in CWM calculations
  column_to_rownames(var = "Species") %>%
  drop_na() #drop NA values so dbFD will run
weights <- kin21 %>% # rows are plots and species are columns
  select_if(colSums(.) != 0) #remove species with 0 total abundance across all plots
traits <- subset(traits, rownames(traits) %in% colnames(weights)) # Subset traits to only species present in abundance data
traits <- traits[order(match(rownames(traits), colnames(weights))), ] # Make sure the two data frames are in the same order

df <- data.frame(Spp = colnames(weights),
                 Trait = colnames(weights) %in% rownames(traits))
drops <- as.character(subset(df, Trait=="FALSE")$Spp) # There is no trait data for these taxa.
weights<-weights[ , !(names(weights) %in% drops)] #drop spp without trait data from abundance df

fdiv_results <- dbFD(traits, weights, stand.x=TRUE, corr="none") # Calculate diversity. "stand.x=TRUE" standardizes traits to mean 0 and unit variance.
##note: I had to drop spp that did not have a trait value for this to work
## can I revisit with a more complete dataset? It should work in theory!

qual.Fric <- fdiv_results$qual.FRic
fdiv <- data.frame( # Extract FDiv metrics to use in downstream analyses.
  richness = fdiv_results$nbsp,
  FRic = fdiv_results$FRic,
  FEve = fdiv_results$FEve,
  FDiv = fdiv_results$FDiv,
  FDis = fdiv_results$FDis,
  RaoQ = fdiv_results$RaoQ
)

CWM <- fdiv_results$CWM %>% # Extract a separate data frame of community-weighted means
  mutate(Plot = 1:204, .before = 1) %>%
  rename(SLA = SLA.Leaf)

## check RTD for outliers
hist(CWM$RTD)

## export results 
write_rds(CWM, "output/July 2021 CWMs.rds")
