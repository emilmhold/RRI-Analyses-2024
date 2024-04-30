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
png("figures/July 2021 kin network_controls_weighted.png", width = 1000, height = 1000)
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
dev.off()


