library("igraph")
library("ggplot2")
library("ggpubr")

stringdb <- read.csv(file="data/stringdb.csv",header = T)
g <- graph_from_data_frame(stringdb, directed=F)
g1 <- data.frame(table(degree(g)))
sum(g1[21:nrow(g1),2])/length(V(g))
#0.1279683

degree_g <- degree(g)
hubs <- names(degree_g[degree_g>20])
hubTable <- NULL
for(hubi in hubs){
  hubline <- c(hubi,mean(stringdb[stringdb$node1%in%hubi | stringdb$node2%in%hubi,4]),
               mean(stringdb[stringdb$node1%in%hubi | stringdb$node2%in%hubi,6]))
  hubTable <- rbind(hubTable,hubline)
}
hubTable <- data.frame(hubTable)
hubTable[,2] <- as.numeric(as.character(hubTable[,2]))
hubTable[,3] <- as.numeric(as.character(hubTable[,3]))
hubTable$avpcc <- (hubTable[,2]+hubTable[,3])/2
colnames(hubTable) <- c("hubs","N0PCCmean","NplusPCCmean","avpcc")

partHubs <- data.frame(hubTable[hubTable$avpcc>median(hubTable$avpcc),])
partHubs <- partHubs[order(partHubs$avpcc,decreasing = T),]
dateHubs <- data.frame(hubTable[hubTable$avpcc<median(hubTable$avpcc),])
dateHubs <- dateHubs[order(dateHubs$avpcc,decreasing = T),]
write.csv(dateHubs,file = "dateHubs.csv",row.names = F)
write.csv(partHubs,file = "partHubs.csv",row.names = F)
#===============================================================
boxplotTable <- data.frame(rbind(cbind(partHubs[,4],"PCChigh"),cbind(dateHubs[,4],"PCClow")))
boxplotTable$X1 <- as.numeric(boxplotTable$X1)
colnames(boxplotTable) <- c("value","group")

p <- ggviolin(boxplotTable , "group", "value", fill = "group",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              add = "boxplot", add.params = list(fill = "white"))+
  geom_signif(comparisons = list(c("PCChigh", "PCClow")),y_position = 0.71)+
  ylab("PCC")
ggsave(p,filename = "groupPCC.pdf",width = 4,height = 4)

boxplotTable <- data.frame(rbind(cbind(degree_g[partHubs[,1]],"PCChigh"),cbind(degree_g[dateHubs[,1]],"PCClow")))
boxplotTable$X1 <- as.numeric(boxplotTable$X1)
colnames(boxplotTable) <- c("value","group")

p<- ggviolin(boxplotTable , "group", "value", fill = "group",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"))+
  geom_signif(comparisons = list(c("PCChigh", "PCClow")),y_position = 150)+
  ylab("degree")
ggsave(p,filename = "groupdegree.pdf",width = 4,height = 4)
#====================================================
#<&
parameters <- read.csv(file="data/parameters.csv",header = T)
datahubs <- read.csv(file="dateHubs.csv",header = T)
parthubs <- read.csv(file="partHubs.csv",header = T)

diffModuleFre<- parameters[parameters$N0NodeModule!=parameters$NplusNodeModule,]
diffModuleFre <- data.frame(table(c(diffModuleFre$node1,diffModuleFre$node2)))
rownames(diffModuleFre) <- diffModuleFre[,1]

PCClow <- (na.omit(diffModuleFre[datahubs$hubs,]))[,2]
PCChigh <- (na.omit(diffModuleFre[parthubs$hubs,]))[,2]

boxplot <- data.frame(group = c(rep("Part hub",length(PCChigh)),rep("Date hub",length(PCClow))),value=c(PCChigh,PCClow))

p<- ggviolin(data = boxplot,x="group",y="value", fill = "group",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = c("jitter","boxplot"), add.params = list(fill = "white"))+
  theme_bw()+ylim(c(-5,50))+
  geom_signif(comparisons = list(c("Part hub", "Date hub")),y_position = 45,test = wilcox.test)+
  ylab("the number of intermodular edges") +
  theme(axis.text.y=element_text(vjust=1,size=15,face = "bold"),
        axis.text.x=element_text(vjust=1,size=15,face = "bold"),
        axis.title.x = element_text(vjust=1,size=15,face = "bold"),
        axis.title.y = element_text(vjust=1,size=15,face = "bold"))


ggsave(p,filename = "BothChangeEdgDegree.pdf",width = 6.5,height = 5.5)

rm(list=ls())
parameters <- read.csv(file="data/parameters.csv",header = T)
datahubs <- read.csv(file="dateHubs.csv",header = T)
parthubs <- read.csv(file="partHubs.csv",header = T)
diffModuleFre<- parameters[parameters$N0NodeModule==parameters$NplusNodeModule,]
diffModuleFre <- data.frame(table(c(diffModuleFre$node1,diffModuleFre$node2)))
rownames(diffModuleFre) <- diffModuleFre[,1]

PCClow <- (na.omit(diffModuleFre[datahubs$hubs,]))[,2]
PCChigh <- (na.omit(diffModuleFre[parthubs$hubs,]))[,2]

boxplot <- data.frame(group = c(rep("PCChigh",length(PCChigh)),rep("PCClow",length(PCClow))),value=c(PCChigh,PCClow))
ggviolin(data = boxplot,x="group",y="value")

p<- ggviolin(data = boxplot,x="group",y="value", fill = "group",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = c("jitter","boxplot"), add.params = list(fill = "white"))+
  geom_signif(comparisons = list(c("PCChigh", "PCClow")),y_position = 100,test = wilcox.test)+
  ylab("the number of changing edge")
p
ggsave(p,filename = "BothKeepEdgDegree.pdf",width = 6,height = 6)
