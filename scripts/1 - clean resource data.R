## clean 2021 soil moisture and ramet density and 2023 productivity data
## Author: Emily H
## Created: December 22, 2023
## Last edited: April 18, 2023

#install.packages("tidyverse")

library(tidyverse)

setwd("~/Documents/GitHub/Analyses/RRI Analyses 2024")

#### import data ####
## import treatment data 
trts <- read.csv("data/Plot treatments.csv") %>%
  dplyr::select(!X) %>%
  mutate(Light = str_replace_all(Light, fixed("Tie-backs"), "Ambient")) ## reassign tie=back plots to the ambient treatment
str(trts)
write_rds(trts, "output/updated_trts.rds")

## soil moisture
moisture21 <- read.csv("data/May 2021 Soil Moisture.csv") %>%
  mutate_at(c('Plot'), as.numeric) %>%
  select(Plot, Gravimetric.soil.moisture) %>%
  rename(Soil.moisture = Gravimetric.soil.moisture)
str(moisture21)

## productivity
productivity23 <- read.csv("data/August 2023 aboveground biomass & litter weights.csv") %>%
  mutate(Productivity = Aboveground.biomass..g./(0.25*0.5)) %>% #divide biomass values by surface area of clipping quadrat (25x50 cm)
  select(Plot, Productivity) %>%
  mutate_at(c('Plot', 'Productivity'), as.numeric) %>%
  filter(Plot != 73) #remove outlier to improve normality of residuals
str(productivity23)
hist(productivity23$Productivity)

## ramet density
density21 <- read.csv("data/August 2021 ramet density.csv") %>%
  mutate_at(c('Plot', 'Ramet.density'), as.numeric)
str(density21)

## richness
## import July 2023 data
richness23 <- read.csv("data/July 2023 % cover.csv") %>%
  mutate_at(c(7:84), as.numeric) %>%
  select(!(c(Date.Collected, Waypoint, Block, Plot.1, Treatment, Bare.Ground, Litter, Cow.pat, Rock))) %>%
  select(!(starts_with("Unk"))) %>%
  rename(Rub.ida2 = Arb.ida) %>%
  mutate(Rub.ida1 = Rub.ida + Rub.ida2) %>% #collapse double observations
  select(!c(Rub.ida, Rub.ida2)) %>%
  rename(Lat.ven = Other.lathyrus,
         Rub.ida = Rub.ida1,
         Cir.und = Cir.umb) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  mutate(Richness = apply(.[2:72] > 0, 1, sum)) %>% #calculate spp richness
  select(Plot, Richness) %>%
  merge(productivity23, by = "Plot")
str(richness23)

## import N data 
N21 <- read.csv("data/2021 Soil N values.csv") %>%
  mutate(NH4.change = August.NH4.g - May.NH4.g) %>% #calculate change over season
  mutate(NO3.change = August.NO3.g - May.NO3.g) %>% #calcualte change over season
  dplyr::select(Plot, NH4.change, NO3.change)
str(N21)
#check correlation of NH4 and NO3
cor(N21[2:3])

## import PAR data
PAR21 <- read.csv("data/July 2021 PAR.csv") %>%
  dplyr::select(Plot, Above.Canopy.PAR, Below.Canopy.PAR) %>%
  merge(trts, by = "Plot")
str(PAR21)

#### prep PAR data ####
par.noshade <- PAR21 %>%
  filter(Light != "Shade") %>% #separate out shade trts
  unite("Treatment", Light:Nutrients:Thin, sep= "/") #reunite trts in one column
str(par.noshade)

#calculate average above-canopy par by block
block.avg.ABCpar <- aggregate(par.noshade$Above.Canopy.PAR, list(par.noshade$Block), mean)
block.avg.ABCpar <- block.avg.ABCpar %>% rename(Avg.Above.Canopy.PAR = x,
                                                Block = Group.1)
str(block.avg.ABCpar)

#replace shade pars with block avg pars
PAR21.new <- PAR21 %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="1", 752.125)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="2", 707.750)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="3", 807.375)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="4", 678.625)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="5", 745.625)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="6", 534.250)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="7", 700.875)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="8", 840.250)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="9", 1208.500)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="10", 1145.875)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="11", 982.375)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="12", 988.500)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="13", 1150.000)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="14", 1368.625)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="15", 1317.250)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="16", 1169.250)) %>%
  mutate(Treatment=replace(Above.Canopy.PAR, Light=="Shade"& Block=="17", 1109.875))

#add new column for light penetration
PAR21.new <- PAR21.new %>%
  dplyr::mutate(True.light.penetration = Below.Canopy.PAR/Above.Canopy.PAR) %>%
  dplyr::select(Plot, True.light.penetration)
str(PAR21.new)

#### merge dfs ####
resources <- trts %>%
  merge(no.spp, by = "Plot") %>% # add richness and productivity data 
  merge(PAR21.new, by = "Plot")  %>% # add PAR data
  merge(N21, by = "Plot")  %>% # add N data
  filter(Light != "Tie-backs") #drop tie-back treatment (ineffective treatment)
str(resources)

#### write files ####
##richness
write_rds(richness23, "output/richness23.rds")

##ramet density
write_rds(density21, "output/density21.rds")

##soil moisture
write_rds(moisture21, "output/moisture21.rds")

## PAR
write_rds(PAR21.new, "output/PAR21 new.rds")

## all resources
write_rds(resources, "output/cleaned 2021 resources.rds")


