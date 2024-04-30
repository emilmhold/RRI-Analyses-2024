## run network ~ resource models
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
## resources
resources <- read_rds("output/cleaned 2021 resources.rds")
str(resources)
cor(resources[9:10]) #NH4 and NO3 are pretty correlated so I'll use NH4

## networks
results_resources_df <- read_rds("output/plot-wise network metrics.rds") %>%
  merge(resources, by = "Plot")
str(results_resources_df)

#### models ####
##connectance
connectance.resources.lmm <- lmer(Connectance ~ True.light.penetration*NH4.change + (1|Block), data = results_resources_df)
conncectance.anova <- Anova(connectance.resources.lmm) %>% # gather results into a table for easy export
  rename(p.value = "Pr(>Chisq)") %>%
  dplyr::select(Chisq, Df, p.value)
connectance.results <- tidy(connectance.resources.lmm, "fixed") %>% # gather results into table
  mutate(Network.metric = "Connectance", .before = 1) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(Network.metric, term, estimate, std.error) %>%
  cbind(conncectance.anova) # cbind anova table
qqnorm(resid(connectance.resources.lmm))
qqline(resid(connectance.resources.lmm))
ks.test(resid(connectance.resources.lmm),"pnorm",mean=mean(resid(connectance.resources.lmm)),sd=sd(resid(connectance.resources.lmm))) 

## nestedness
nestedness.resources.lmm <- lmer(Nestedness ~ True.light.penetration*NH4.change + (1|Block), data = results_resources_df)
nestedness.anova <- Anova(nestedness.resources.lmm) %>% 
  rename(p.value = "Pr(>Chisq)") %>%
  dplyr::select(Chisq, Df, p.value)
nestedness.results <- tidy(nestedness.resources.lmm, "fixed") %>%
  mutate(Network.metric = "Nestedness", .before = 1) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(Network.metric, term, estimate, std.error) %>%
  cbind(nestedness.anova)
qqnorm(resid(nestedness.resources.lmm))
qqline(resid(nestedness.resources.lmm))
ks.test(resid(nestedness.resources.lmm),"pnorm",mean=mean(resid(nestedness.resources.lmm)),sd=sd(resid(nestedness.resources.lmm))) 

## clustering
clustering.resources.lmm <- lmer(Clustering ~ True.light.penetration*NH4.change + (1|Block), data = results_resources_df)
clustering.anova <- Anova(clustering.resources.lmm) %>% 
  rename(p.value = "Pr(>Chisq)") %>%
  dplyr::select(Chisq, Df, p.value)
clustering.results <- tidy(clustering.resources.lmm, "fixed") %>%
  mutate(Network.metric = "Clustering", .before = 1) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(Network.metric, term, estimate, std.error) %>%
  cbind(clustering.anova)
qqnorm(resid(clustering.resources.lmm))
qqline(resid(clustering.resources.lmm))
ks.test(resid(clustering.resources.lmm),"pnorm",mean=mean(resid(clustering.resources.lmm)),sd=sd(resid(clustering.resources.lmm))) 

#### combine and export tables ####
network.resource.results <- rbind(connectance.results, nestedness.results, clustering.results) %>%
  remove_rownames() %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  mutate("Estimate ± standard error" = paste(estimate, std.error, sep = " ± ")) %>%
  rename(Resource = term) %>%
  mutate(Resource = str_replace_all(Resource, fixed("True.light.penetration"), "Light penetration")) %>%
  mutate(Resource = str_replace_all(Resource, fixed("NH4.change"), "Soil NH4")) %>%
  dplyr::select(Network.metric, Resource, "Estimate ± standard error", Chisq, Df, p.value)
write_csv(network.resource.results, file = "tables/Network ~ resources lmm results.csv")

#### figures ####
resource_results.forplot <- results_resources_df %>%
  pivot_longer(cols = c(Connectance, Nestedness, Clustering), names_to = "Network.metric", values_to = "Network.value") %>%
  pivot_longer(cols = c(True.light.penetration, NH4.change), names_to = "Resource", values_to = "Resource.value") %>%
  mutate(Resource = str_replace_all(Resource, fixed("True.light.penetration"), "Light penetration (PPFD)")) %>%
  mutate(Resource = str_replace_all(Resource, fixed("NH4.change"), "Soil NH4 (mg/L)"))
str(resource_results.forplot)

resource.plot <- ggplot(resource_results.forplot, aes(x = Resource.value, y = Network.value)) +
  geom_point() + 
  facet_wrap(~Network.metric + Resource, scales = 'free', ncol = 2) +
  ylab("Network metric") +
  xlab("Resource value") + theme_classic(base_size = 18)
resource.plot
ggsave(filename = "Network metrics by resources.png", 
       resource.plot,
       path = "figures/",
       width = 16,
       height = 16,
       units = "in"
)

