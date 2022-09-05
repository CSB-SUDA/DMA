library(igraph)
library(tidyverse)
library(dplyr)
source("11.EM.R")

parameters <- read_csv(file="data/parameters.csv")
parameters <- parameters[parameters$N0NodeModule!=parameters$NplusNodeModule,]
dataHubs <- read_csv(file="data/dateHubs.csv")
parameters <- parameters[parameters$node1%in%dataHubs$hubs & parameters$node2%in%dataHubs$hubs,]

g <- graph_from_data_frame(parameters,directed = F)

tfcnoweight <- parameters %>% dplyr::select(node1,node2) %>% TFC
write.csv(tfcnoweight %>% arrange(desc(TFCs)),file="TFC_noweight.csv",row.names = F,quote = F)
