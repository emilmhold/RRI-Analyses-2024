## clean 2021 soil moisture and ramet density and 2023 productivity data
## Author: Emily H
## Created: December 22, 2023
## Last edited: April 18, 2023

#install.packages("tidyverse")

library(tidyverse)

setwd("~/Documents/GitHub/Analyses/RRI Analyses 2024")

#### import data ####
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

#### write files ####
##richness
write_rds(richness23, "output/richness23.rds")

##ramet density
write_rds(density21, "output/density21.rds")

##soil moisture
write_rds(moisture21, "output/moisture21.rds")
