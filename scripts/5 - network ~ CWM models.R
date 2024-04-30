## run network ~ CWM models
## Author: Emily H
## Created: April 18, 2023
## Last edited: April 21, 2023

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
cor(CWM[2:5])

## networks
results_df <- read_rds("output/plot-wise network metrics.rds") %>%
  merge(CWM, by = "Plot") %>%
  merge(trts, by = "Plot")

#### models ####
## connectance
connectance.CWM.lmm <- lmer(Connectance ~ Max.Height*SLA*SRL*RTD + (1|Block), data = results_df)
conncectance.anova <- Anova(connectance.CWM.lmm) %>% # gather results into a table for easy export
  rename(p.value = "Pr(>Chisq)") %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  dplyr::select(Chisq, Df, p.value)
str(conncectance.anova)
connectance.results <- tidy(connectance.CWM.lmm, "fixed") %>% # gather results into table
  mutate(Network.metric = "Connectance", .before = 1) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(Network.metric, term, estimate, std.error) %>%
  cbind(conncectance.anova) %>% # cbind anova table
  as.data.frame()
str(connectance.results)
qqnorm(resid(connectance.CWM.lmm))
qqline(resid(connectance.CWM.lmm))
ks.test(resid(connectance.CWM.lmm),"pnorm",mean=mean(resid(connectance.CWM.lmm)),sd=sd(resid(connectance.CWM.lmm))) 

## nestedness
nestedness.CWM.lmm <- lmer(Nestedness ~ Max.Height*SLA*SRL*RTD + (1|Block), data = results_df)
nestedness.anova <- Anova(nestedness.CWM.lmm) %>% 
  rename(p.value = "Pr(>Chisq)") %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  dplyr::select(Chisq, Df, p.value)
nestedness.results <- tidy(nestedness.CWM.lmm, "fixed") %>%
  mutate(Network.metric = "Nestedness", .before = 1) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(Network.metric, term, estimate, std.error) %>%
  cbind(nestedness.anova) %>%
  as.data.frame()
qqnorm(resid(nestedness.CWM.lmm))
qqline(resid(nestedness.CWM.lmm))
ks.test(resid(nestedness.CWM.lmm),"pnorm",mean=mean(resid(nestedness.CWM.lmm)),sd=sd(resid(nestedness.CWM.lmm))) 

## clustering
clustering.CWM.lmm <- lmer(Clustering ~ Max.Height*SLA*SRL*RTD + (1|Block), data = results_df)
clustering.anova <- Anova(clustering.CWM.lmm) %>% 
  rename(p.value = "Pr(>Chisq)") %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  dplyr::select(Chisq, Df, p.value)
clustering.results <- tidy(clustering.CWM.lmm, "fixed") %>%
  mutate(Network.metric = "Clustering", .before = 1) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(Network.metric, term, estimate, std.error) %>%
  cbind(clustering.anova) %>%
  as.data.frame()
qqnorm(resid(clustering.CWM.lmm))
qqline(resid(clustering.CWM.lmm))
ks.test(resid(clustering.CWM.lmm),"pnorm",mean=mean(resid(clustering.CWM.lmm)),sd=sd(resid(clustering.CWM.lmm))) 

#### combine and export tables ####
network.CWM.results <- rbind(connectance.results, nestedness.results, clustering.results) %>%
  remove_rownames() %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  mutate("estimate ± standard error" = paste(estimate, std.error, sep = " ± ")) %>%
  rename(CWM.trait = term) %>%
  dplyr::select(Network.metric, CWM.trait, "estimate ± standard error", Chisq, Df, p.value)
write_csv(network.CWM.results, file = "tables/Network ~ CWM lmm results.csv")
str(network.CWM.results)

#### figures ####
CWM_results.forplot <- results_df %>%
  pivot_longer(cols = c(Max.Height, SLA, SRL, RTD), names_to = "Trait", values_to = "CWM.value") %>%
  pivot_longer(cols = c(Connectance, Nestedness, Clustering), names_to = "Network.metric", values_to = "Network.value") %>%
  mutate(Trait = str_replace_all(Trait, fixed("Max.Height"), "Maximum height (cm)")) %>%
  mutate(Trait = str_replace_all(Trait, fixed("SLA"), "SLA (cm2/g)")) %>%
  mutate(Trait = str_replace_all(Trait, fixed("SRL"), "SRL (cm/mg)")) %>%
  mutate(Trait = str_replace_all(Trait, fixed("RTD"), "RTD (mg/cm2)"))
str(CWM_results.forplot)

CWM.plot <- ggplot(CWM_results.forplot, aes(x = CWM.value, y = Network.value)) +
  geom_point() + 
  facet_wrap(~Trait + Network.metric, scales = 'free', ncol = 3) +
  ylab("Network metric") +
    xlab("CWM trait value") + 
  geom_smooth(data = subset(CWM_results.forplot, Trait == "Maximum height (cm)" & Network.metric == "Connectance"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(CWM_results.forplot, Trait == "SLA (cm2/g)" & Network.metric == "Connectance"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(CWM_results.forplot, Trait == "Maximum height (cm)" & Network.metric == "Nestedness"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(CWM_results.forplot, Trait == "SLA (cm2/g)" & Network.metric == "Nestedness"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(CWM_results.forplot, Trait == "Maximum height (cm)" & Network.metric == "Clustering"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(CWM_results.forplot, Trait == "SLA (cm2/g)" & Network.metric == "Clustering"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(CWM_results.forplot, Trait == "RTD (mg/cm2)" & Network.metric == "Clustering"), method = "lm", se = FALSE) +
  theme_classic(base_size = 18)
CWM.plot
ggsave(filename = "Network metrics by CWM.png", 
       CWM.plot,
       path = "figures/",
       width = 16,
       height = 16,
       units = "in"
)
