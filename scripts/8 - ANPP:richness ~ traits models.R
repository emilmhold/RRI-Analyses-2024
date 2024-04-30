## run ANPP/richness ~ traits models
## Author: Emily H
## Created: April 18, 2023
## Last edited: April 18, 2023

#install.packages("tidyverse")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("car")
#install.packages("broom.mixed")

library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(broom.mixed)

setwd("~/Documents/GitHub/Analyses/RRI Analyses 2024")

options(scipen = 999) #turn off scientific notation

#### import data ####
##treatments
trts <- read_rds("output/updated_trts.rds")

## CWMs
CWM <- read_rds("output/July 2021 CWMs.rds")

## richness and ANPP
richness23 <- read_rds(file = "output/richness23.rds") %>%
  merge(trts, by = "Plot") %>%
  merge(CWM, by = "Plot")
str(richness23)

#### models ####
## richness
richness.CWM.lmm <- lmer(Richness ~ Max.Height*SLA*SRL*RTD + (1|Block), data = richness23)
richness.anova <- Anova(richness.CWM.lmm) %>% # gather results into a table for easy export
  rename(p.value = "Pr(>Chisq)") %>%
  rownames_to_column(var = "term")
str(richness.anova)
richness.results <- tidy(richness.CWM.lmm, "fixed") %>% # gather results into table
  mutate(EGS = "Richness", .before = 1) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(EGS, term, estimate, std.error) %>%
  full_join(richness.anova, by = "term")
str(richness.CWM.lmm)
qqnorm(resid(richness.CWM.lmm))
qqline(resid(richness.CWM.lmm))
ks.test(resid(richness.CWM.lmm),"pnorm",mean=mean(resid(richness.CWM.lmm)),sd=sd(resid(richness.CWM.lmm))) 

## productivity
productivity.CWM.lmm <- lmer(Productivity ~ Max.Height*SLA*SRL*RTD + (1|Block), data = richness23)
productivity.anova <- Anova(productivity.CWM.lmm) %>% # gather results into a table for easy export
  rename(p.value = "Pr(>Chisq)") %>%
  rownames_to_column(var = "term")
productivity.results <- tidy(productivity.CWM.lmm, "fixed") %>% # gather results into table
  mutate(EGS = "Productivity", .before = 1) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(EGS, term, estimate, std.error) %>%
  full_join(productivity.anova, by = "term")
qqnorm(resid(productivity.CWM.lmm))
qqline(resid(productivity.CWM.lmm))
ks.test(resid(productivity.CWM.lmm),"pnorm",mean=mean(resid(productivity.CWM.lmm)),sd=sd(resid(productivity.CWM.lmm))) 

#### combine and export tables ####
EGS.CWM.results <- rbind(richness.results, productivity.results) %>%
  remove_rownames() %>%
  write_csv(file = "tables/EGS ~ CWM lmm results.csv")

#### figures ####
EGS_results.forplot <- richness23 %>%
  pivot_longer(cols = c(Max.Height, SLA, SRL, RTD), names_to = "Trait", values_to = "CWM.value") %>%
  pivot_longer(cols = c(Richness, Productivity), names_to = "EGS", values_to = "EGS.value")
str(EGS_results.forplot)

EGS.plot <- ggplot(EGS_results.forplot, aes(x = CWM.value, y = EGS.value)) +
  geom_point() + 
  facet_wrap(~Trait + EGS, scales = 'free', ncol = 2) +
  xlab("CWM trait value") +
  ylab("EGS value") + 
  geom_smooth(data = subset(EGS_results.forplot, EGS == "Productivity" & Trait == "Max.Height"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(EGS_results.forplot, EGS == "Productivity" & Trait == "SRL"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(EGS_results.forplot, EGS == "Productivity" & Trait == "RTD"), method = "lm", se = FALSE) +
  theme_classic(base_size = 18)
EGS.plot
ggsave(filename = "EGS by CWMs.png", 
       EGS.plot,
       path = "figures/",
       width = 16,
       height = 20,
       units = "in"
)

