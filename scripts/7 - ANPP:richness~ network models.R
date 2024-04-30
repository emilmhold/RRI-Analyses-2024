## run ANPP/richness ~ network models
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

## richness and ANPP
richness23 <- read_rds(file = "output/richness23.rds")

## networks
results_richness_df <- read_rds("output/plot-wise network metrics.rds") %>%
  merge(trts, by = "Plot") %>%
  merge(richness23, by = "Plot")
str(results_richness_df)

#### models ####
## richness
richness.network.lmm <- lmer(Richness ~ Connectance*Nestedness*Clustering + (1|Block), data = results_richness_df)
richness.anova <- Anova(richness.network.lmm) %>%  # gather results into a table for easy export
  rename(p.value = "Pr(>Chisq)") %>%
  rownames_to_column(var = "term")
richness.results <- tidy(richness.network.lmm, "fixed") %>% # gather results into table
  mutate(EGS = "Species richness", .before = 1) %>%
  mutate(term = str_replace_all(term, fixed("(Intercept)"), "Nestedness")) %>%
  dplyr::select(EGS, term, estimate, std.error) %>%
  full_join(richness.anova, by = "term")
str(richness.results)
qqnorm(resid(richness.network.lmm))
qqline(resid(richness.network.lmm))
ks.test(resid(richness.network.lmm),"pnorm",mean=mean(resid(richness.network.lmm)),sd=sd(resid(richness.network.lmm))) 

## productivity
productivity.network.lmm <- lmer(Productivity ~ Connectance*Nestedness*Clustering + (1|Block), data = results_richness_df)
productivity.anova <- Anova(productivity.network.lmm) %>%  # gather results into a table for easy export
  rename(p.value = "Pr(>Chisq)") %>%
  rownames_to_column(var = "term")
productivity.results <- tidy(productivity.network.lmm, "fixed") %>% # gather results into table
  mutate(EGS = "Productivity", .before = 1) %>%
  mutate(term = str_replace_all(term, fixed("(Intercept)"), "Nestedness")) %>%
  dplyr::select(EGS, term, estimate, std.error) %>%
  full_join(productivity.anova, by = "term")
str(productivity.results)
qqnorm(resid(productivity.network.lmm))
qqline(resid(productivity.network.lmm))
ks.test(resid(productivity.network.lmm),"pnorm",mean=mean(resid(productivity.network.lmm)),sd=sd(resid(productivity.network.lmm))) 

#### combine and export tables ####
EGS.network.results <- rbind(richness.results, productivity.results) %>%
  remove_rownames() %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  mutate("Estimate ± standard error" = paste(estimate, std.error, sep = " ± ")) %>%
  rename("Network metric" = term) %>%
  mutate(EGS = str_replace_all(EGS, fixed("Productivity"), "ANPP")) %>%
  dplyr::select(EGS, "Network metric", "Estimate ± standard error", Chisq, Df, p.value) %>%
  write_csv(file = "tables/EGS ~ network lmm results.csv")

#### figures ####
EGS_results.forplot <- results_richness_df %>%
  pivot_longer(cols = c(Connectance, Nestedness, Clustering), names_to = "Network.metric", values_to = "Network.value") %>%
  pivot_longer(cols = c(Richness, Productivity), names_to = "EGS", values_to = "EGS.value") %>%
  mutate(EGS = str_replace_all(EGS, fixed("Productivity"), "ANPP (g/m2)"))
str(EGS_results.forplot)

EGS.plot <- ggplot(EGS_results.forplot, aes(x = Network.value, y = EGS.value)) +
  geom_point() + 
  facet_wrap(~Network.metric + EGS, scales = 'free', ncol = 2) +
  xlab("Network metric") +
  ylab("EGS variable") + 
  geom_smooth(data = subset(EGS_results.forplot, EGS == "Richness" & Network.metric == "Connectance"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(EGS_results.forplot, EGS == "Richness" & Network.metric == "Clustering"), method = "lm", se = FALSE) +
  geom_smooth(data = subset(EGS_results.forplot, EGS == "Richness" & Network.metric == "Nestedness"), method = "lm", se = FALSE) +
  theme_classic(base_size = 18)
EGS.plot
ggsave(filename = "EGS by network metrics.png", 
       EGS.plot,
       path = "figures/",
       width = 16,
       height = 16,
       units = "in"
)

