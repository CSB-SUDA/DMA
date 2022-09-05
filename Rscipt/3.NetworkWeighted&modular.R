library(igraph)
library(reshape2)
library(ggpubr)
library(ggplot2)

TCGAexprs <- read.csv(file="data/BRCA_log2TPM_DEGmatrix.csv",
                                         row.names = 1,header = T,check.names = F)
#=====
stringdb <- read.delim(file="data/string_interactions_short_0_4.tsv",header=T)
stringdb <- stringdb[,c(1,2)]
colnames(stringdb) <- c("node1","node2")
clinical <- read.csv(file="data/BRCA_clinical.csv",header=T,row.names = 1)
stringdb <- stringdb[(stringdb$node1%in%rownames(TCGAexprs)) & (stringdb$node2%in%rownames(TCGAexprs)),]
#=====
for(linei in 1:nrow(stringdb)){
    cat(paste0(linei,"\r"))
    stringdb$N0cor[linei] <- cor(as.numeric(TCGAexprs[stringdb[linei,1],clinical$N0N=="N0"]),
                    as.numeric(TCGAexprs[stringdb[linei,2],clinical$N0N=="N0"]))
    stringdb$N0cororigin[linei] <- stringdb$N0cor[linei]
    stringdb$Npluscor[linei] <- cor(as.numeric(TCGAexprs[stringdb[linei,1],clinical$N0N!="N0"]),
                          as.numeric(TCGAexprs[stringdb[linei,2],clinical$N0N!="N0"]))
    stringdb$Npluscororigin[linei] <- stringdb$Npluscor[linei]
}
stringdb$N0cor <- 0.5*log((1+stringdb$N0cor)/(1-stringdb$N0cor))
stringdb$Npluscor <- 0.5*log((1+stringdb$Npluscor)/(1-stringdb$Npluscor))
write.csv(stringdb,file="stringdb.csv",quote = F,row.names = F)

nodes <- unique(c(stringdb$node1,stringdb$node2))
net_pc<-graph_from_data_frame(d=stringdb,vertices=nodes,directed=F)

set.seed(1)
ecN0 <- multilevel.community(net_pc,weights = abs(stringdb$N0cor))
set.seed(1)
ecNplus <- multilevel.community(net_pc,weights = abs(stringdb$Npluscor))
length(ecN0)
save(ecN0,ecNplus,file = "ecNplus-ecN0.RData")
length(ecNplus)
#================================================================
ecN0Cluster <- NULL
for(i in 1:length(ecN0)){
    if(length(ecN0[[i]])>=0)
    ecN0Cluster <- rbind(ecN0Cluster,cbind(ecN0[[i]],paste0("cluster",i)))
}
write.csv(ecN0Cluster,file="ecN0Clusterlg10.csv",quote = F,row.names = F)
ecNplusCluster <- NULL
for(i in 1:length(ecNplus)){
    if(length(ecNplus[[i]])>=0)
    ecNplusCluster <- rbind(ecNplusCluster,cbind(ecNplus[[i]],paste0("cluster",i)))
}
write.csv(ecNplusCluster,file="ecNplusClusterlg10.csv",quote = F,row.names = F)
#====================================================================
library(networkD3)

N0data <- ecN0Cluster
colnames(N0data) <- c("gene","N0")
Nplusdata <- ecNplusCluster
colnames(Nplusdata) <- c("gene","Nplus")
data <- merge(N0data,Nplusdata,by="gene",all = T)
data$N0 <- paste0("N0_",data$N0)
data$Nplus <- paste0("Nplus_",data$Nplus)
Nodes <- unique(c(data$N0,data$Nplus))
Nodes <- data.frame(name=Nodes)
data$N0 <- match(data$N0,Nodes$name)-1
data$Nplus <- match(data$Nplus,Nodes$name)-1
data <- data[,2:3]
links <- data.frame(reshape2::melt(table(data)))
links <- links[links$value!=0,]

p<- sankeyNetwork(Links = links, Nodes = Nodes, Source = 'N0',
              Target = 'Nplus', Value = 'value', NodeID = 'name',
              units = 'TWh', fontSize = 12, nodeWidth = 120)
saveNetwork(p,"sankey.html")
