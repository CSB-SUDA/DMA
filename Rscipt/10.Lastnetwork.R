library(igraph)
library(tidyverse)
library(dplyr)
source("EM.R")

parameters <- read_csv(file="data/parameters.csv")
parameters <- parameters[parameters$N0NodeModule!=parameters$NplusNodeModule,]
dataHubs <- read_csv(file="data/dateHubs.csv")
parameters <- parameters[parameters$node1%in%dataHubs$hubs & parameters$node2%in%dataHubs$hubs,]
parameters$survivalRelated <- parameters$N0pvalue<0.05 | parameters$Npluspvalue<0.05 | parameters$allpvalue<0.05
parameters$survivalGene <-  parameters$allpvalue<0.05

g <- graph_from_data_frame(parameters,directed = F)

parameters <- parameters %>% add_column(edge.betweenness_N0cor = edge_betweenness(g,weights = abs(parameters$N0cor))) %>%
              add_column(edge.betweenness_Npluscor = edge_betweenness(g,weights = abs(parameters$Npluscor)))
parameters$squareEdgeBetweeness <- parameters$edge.betweenness_N0cor/2+parameters$edge.betweenness_Npluscor/2
write.csv(parameters %>% dplyr::select(node1,node2,squareEdgeBetweeness) %>% arrange(desc(squareEdgeBetweeness)),file="edgeBetw_weight.csv",row.names = F,quote = F)

parameters$edge.betweennes <- edge_betweenness(g)
write.csv(parameters %>% dplyr::select(node1,node2,edge.betweennes) %>% arrange(desc(edge.betweennes)),file="edgeBetw_noweight.csv",row.names = F,quote = F)

N0tfc <- parameters %>% dplyr::select(node1,node2) %>% TFC(weights=abs(parameters$N0cororigin))
Nplustfc <- parameters %>% dplyr::select(node1,node2) %>% TFC(weights=abs(parameters$Npluscororigin))
tfcweight <- merge(N0tfc,Nplustfc,by= c("name1","name2"))
tfcweight$tfc.weight <- (tfcweight[,3]+tfcweight[,4])/2
write.csv(tfcweight %>% dplyr::select(name1,name2,tfc.weight) %>% arrange(desc(tfc.weight)),file="TFC_weight.csv",row.names = F,quote = F)

tfcnoweight <- parameters %>% dplyr::select(node1,node2) %>% TFC
write.csv(tfcnoweight %>% arrange(desc(TFCs)),file="TFC_noweight.csv",row.names = F,quote = F)
