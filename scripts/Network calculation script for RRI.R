## make co-occurrence networks
## Author: Emily H
## Created: November 30, 2023
## Last edited: April 21, 2024

#install.packages("cooccur")
#install.packages("tidyverse")
#install.packages("igraph")
#install.packages("vegan")
#install.packages("viridis")

library(cooccur)
library(tidyverse)
library(igraph)
library(vegan)
library(viridis)

setwd("~/Documents/GitHub/Analyses/RRI Analyses 2024")

#### import data ####
kin <- read.csv("data/July 2021 % cover.csv")
kin[kin == "t"] <- '1'
kin <- kin %>%
  mutate_at(c(6:76), as.numeric) %>%
  select(!c(Date.Collected, Block, Waypoint, Treatment, Bare.Ground, Litter, moss, Cow.pat, X)) %>%
  rename(Hel.hoo2 = Heli.hoo) %>%
  mutate(Hel.hoo1 = Hel.hoo + Hel.hoo2,
         Ery.che1 = Ery.che + rough.mustard) %>% #collapse double observations
  select(!c(Hel.hoo, Hel.hoo2, fuzz.auricle.grass, slender.wheat.grass, Ery.che, rough.mustard)) %>%
  select(!(starts_with("un"))) %>%
  rename(Ago.gla = Tradub.like,
         Hel.hoo = Hel.hoo1,
         Ery.che = Ery.che1,
         Pas.smi = Pas.mi) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  select(.,order(colnames(.))) #put columns in alphabetical order
str(kin) #60 species
## 48 species would be 80% of the total 
## only 24 spp are predicted/found to have >1 cooccurrence in communities

#import treatment data
trts <- read.csv("data/Plot treatments.csv") %>%
  select(!X)
str(trts)

#prep data for networks
kin.net <- trts %>%
  merge(kin, by = "Plot") %>%
  filter(Light == "Ambient" & Nutrients == "No fert" & Thin == "Not thinned") %>% #select control plots
  select(!c(Block, Light, Nutrients, Thin)) %>% #remove metadata
  column_to_rownames(var = "Plot") %>%
  t() %>%
  decostand(method = "pa") #convert to presence/absence
str(kin.net)


#### create cooccurrence matrices ####
cooccur.kin <- cooccur(kin.net, type = "spp_site", thresh = TRUE, spp_names = TRUE) #note: thresh = TRUE means spp pairs 
#expected to have less than one co-occurrence are filtered from the dataset
summary(cooccur.kin)
# open figure 
#png("figures/July 2021 kin cooccurrence matrix.png", width = 1000, height = 1000)
plot(cooccur.kin) 
#close the file
#dev.off()
result.kin <- prob.table(cooccur.kin)

spp_names <- c("Ach.mil") %>% #start with Ach.mil which is missing from sp2_name vector
  append(unique(result.kin$sp2_name)) #extract spp names for network vertices
spp_names

#### come back and check this ####
result.kin$sign<-NA #create column to denote sign of cooccurrences
result.kin$sign[result.kin$p_lt<result.kin$p_gt]<-"N" # negative cooccurrences
result.kin$sign[result.kin$p_gt<result.kin$p_lt]<-"P" # positive cooccurrences

write_rds(result.kin, "output/aug 2021 kin cooccurrence matrix.rds")

#result.kin <- result.kin %>%
#select(sp1_name,sp2_name,sig)%>%
#filter(!is.na(sig)) %>%
#mutate(colour = if_else(sig == "P", "gray80", "red")) 

#count how many positive vs. negative cooccurrences
#cooccur.count <- result.kin %>%
#group_by(sp1_name) %>%
#summarise(count = length(sp1_name))
#cooccur.count

## 7 spp @ alpha = 0.05
## 10 spp @ alpha = 0.1
## 16 spp @ alpha = 0.3
## 22 spp @ alpha = 0.5 
## 24 spp @ alpha = 0.6

#### make networks ####
## in network
result.kin_for.graph <- result.kin %>% select(sp1_name, sp2_name, prob_cooccur, sign)
net.kin <- graph_from_data_frame(d = result.kin_for.graph, directed = FALSE)
E(net.kin)$weight <- result.kin$prob_cooccur #set edge weights as probability of spp cooccurrences
net.kin[sparse = FALSE] #prints adjacency matrix with 0s and 1s

net.kin
plot(net.kin)


col.vir<-viridis(4)

cnodes <- rep("gray80",vcount(net.kin)) #node colours: repeat gray80 for the count of vertices in net.kin
#cnodes[nodes$provenance=='introduced']<- col.vir[1] #purple
#cnodes[nodes$provenance=='Both']<-col.vir[4] #yellow
#V(net)$label <- nodes$USDA_code
E(net.kin)$width <- E(net.kin)$weight*10 #set edge width
ecol <- rep("gray80", ecount(net.kin)) #edges colours: repeat gray 80 for the count of edges in net.kin
ecol[(result.kin$sign=='P')] <- col.vir[3] #if positive cooccurrence, then make edge green
ecol[(result.kin$sign=='N')] <- col.vir[2] #if negative cooccurrence, then make edge blue

## make and save figure
# open figure 
#png("figures/July 2021 kin network_controls_weighted.png", width = 1000, height = 1000)
#plot
plot.igraph(net.kin,edge.color=ecol,vertex.color=cnodes, vertex.frame.color=cnodes,
            vertex.label.color="black", #vertex.label=V(net)$USDA_code,
            layout=layout.circle(net.kin),vertex.size=10, vertex.label.family= "Helvetica",
            vertex.label.cex = 1.5, main = "July 2021 Kinsella network (control plots only)",
            cex.main = 2) 
#legend
legend("topright", legend=c("positive", "negative"),seg.len=1,
       col=c(col.vir[3], col.vir[2]), lty=1, cex=2, lwd=4,
       title="Co-occurence", box.lty=0)
#close the file
#dev.off()

#### create plot-level networks and calculate network metrics ####
#### import data ####
## import species richness and productivity data
no.spp <- read_rds("output/richness23.rds")

#### calculate edge weights ####
cooccur.count <- result.kin %>%
  group_by(sp1_name) %>%
  mutate(weighted.sign = if_else(sign == "P", 1, -1)) %>%
  mutate(weighted.prob = prob_cooccur*weighted.sign) %>%
  select(sp1_name, sp2_name, weighted.prob, sign) %>%
  na.omit()
str(cooccur.count)

## find which cover species are not in cooccurrence df and remove
df <- data.frame(Spp = colnames(kin_percent),
                 interactions = colnames(kin_percent) %in% cooccur.count$sp1_name)
drops <- as.character(subset(df, interactions=="FALSE")$Spp) # There is no cooccurrence data for these taxa.
drops

#remove spp without cooccurrence data
kin_percent2<-kin_percent[ , !(names(kin_percent) %in% drops)]
str(kin_percent2) ## 23 spp - it worked!

#### change kin percent from wide to long ####
kin_percent_long <- kin_percent2 %>%
  rownames_to_column(var = "Plot") %>%
  pivot_longer(cols = 2:24, 
               cols_vary = "slowest", 
               names_to = "sp1_name",
               values_to = "sp1_abundance") %>% 
  filter(sp1_abundance > 0) #remove zeroes

## create sp2_name column and sp2_name abundance
#create df to draw from for sp2_abundance
kin_percent_long3 <- kin_percent2 %>%
  rownames_to_column(var = "Plot") %>%
  pivot_longer(cols = 2:24, 
               cols_vary = "slowest", 
               names_to = "sp2_name",
               values_to = "sp2_abundance") #spread data
str(kin_percent_long3)

kin_percent_long2 <- kin_percent_long %>%
  full_join(cooccur.count, by = "sp1_name", multiple = "all") %>% #merge with cooccur.count to input weighted prob values
  merge(kin_percent_long3, by = c("Plot","sp2_name")) %>%
  filter(sp2_abundance > 0) %>% #remove zeroes
  dplyr::select(Plot, sp1_name, sp1_abundance, sp2_name, sp2_abundance, weighted.prob, sign) #reorder columns
str(kin_percent_long2)

#### calculate new edge weights for graphs ####
plot.edges <- kin_percent_long2 %>%
  mutate(plot.weighted.edge = sp1_abundance*sp2_abundance*weighted.prob) %>% 
  #calculate each edge as the product of the two species' abundances*the edge weight from the kinsella network
  dplyr::select(Plot, sp1_name, sp2_name, plot.weighted.edge)

#### write for loop to create networks and calculate metrics for each plot ####
# write nestedness function
calculate_nestedness <- function(net) {
  adj_matrix <- as_adjacency_matrix(net, sparse = FALSE)
  nestedness_rows <- apply(adj_matrix, 1, function(row) sum(row != 0) / length(row))
  nestedness_cols <- apply(adj_matrix, 2, function(col) sum(col != 0) / length(col))
  nestedness <- (sum(nestedness_rows) + sum(nestedness_cols)) / (2 * min(length(nestedness_rows), length(nestedness_cols)))
  return(nestedness)
}

## trying out for loop functions with plot 2
plot2.df <- plot.edges %>% filter(Plot == "2") %>% dplyr::select(-Plot) #make df of plot1 data
net <- graph_from_data_frame(d = plot2.df, directed = FALSE) #make network
E(net)$weight <- plot2.df$plot.weighted.edge #set edge weights as probability of spp cooccurrences
clustering <- transitivity(net, type = "global")
net[sparse = FALSE] #adjacency matrix

##create empty df
results_df = data.frame()

for (i in 1:204) {
  # Filter data for the current plot
  plot_df <- plot.edges %>% 
    filter(Plot == as.character(i)) %>% 
    select(-Plot)
  
  # Make network from plot data
  net <- graph_from_data_frame(d = plot_df, directed = FALSE)
  
  # Set edge weights as probability of species co-occurrences
  E(net)$weight <- plot_df$plot.weighted.edge
  
  # Calculate connectance
  connectance <- ecount(net) / vcount(net)^2 
  
  # Calculate nestedness
  nestedness <- calculate_nestedness(net)
  
  # Calculate clustering coefficients
  clustering <- transitivity(net, type = "global")
  
  # Append results to the dataframe
  results_df <- bind_rows(results_df, data.frame(Plot = i, Connectance = connectance, Nestedness = nestedness, Clustering = clustering))
}

##export df 
str(results_df)
write_rds(results_df, "output/plot-wise network metrics.rds")

#### figures ####
## change data into long format
results_df_new <- results_df %>%
  pivot_longer(cols = c(Connectance, Nestedness, Clustering), names_to = "Network.metric", values_to = "Network.value")
str(results_df_new)  

network.histograms <- ggplot(data = results_df_new, aes(x = Network.value)) +
  geom_histogram() +
  facet_wrap(~Network.metric, scales = 'free', ncol = 3) +
  xlab("Network metric") + 
  ylab("Frequency") +
  theme_classic(base_size = 18)
network.histograms
ggsave(filename = "Network metric histograms.png", 
       network.histograms,
       path = "figures/",
       width = 16,
       height = 16,
       units = "in"
)
