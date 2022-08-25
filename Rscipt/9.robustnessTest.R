library("ggpubr")

rm(list=ls())
parameters <- read.csv(file="data/parameters.csv",header = T)
sameModule <- parameters[parameters$N0NodeModule==parameters$NplusNodeModule,]
diffModule<- parameters[parameters$N0NodeModule!=parameters$NplusNodeModule,]
datahubs <- read.csv(file="data/dateHubs.csv",header = T)
parthubs <- read.csv(file="data/partHubs.csv",header = T)
sameModule <- sameModule[sameModule$node1%in%datahubs$hubs | sameModule$node2%in%datahubs$hubs,]
diffModule <- diffModule[diffModule$node1%in%datahubs$hubs | diffModule$node2%in%datahubs$hubs,]

diffModule[,4] <- abs(diffModule[,4])
parameters[,4] <- abs(parameters[,4])
PCClowdiffsample <- NULL
for(i in 1:100000){
  cat(paste0(i,"\r"))
  PCClowdiffsample <- c(PCClowdiffsample,mean(parameters[sample(1:nrow(parameters),nrow(diffModule),replace = F),4]))
}
Freplot <- data.frame(table(cut(PCClowdiffsample,
                                breaks = seq(mean(diffModule[,4])-0.01,
                                             max(PCClowdiffsample),
                                             (max(PCClowdiffsample)-mean(diffModule[,4]))/100))))

colnames(Freplot) <- c("Interval","Frequency")
ggplot(data = Freplot,aes(x=Interval,y=Frequency))+
  geom_col(width=0.9)+theme_bw()+
  geom_vline(aes(xintercept=8),colour="#990000",linetype="dashed",size=1)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(breaks=c("(0.225,0.227]","(0.268,0.269]","(0.3,0.301]","(0.328,0.329]"))+
  annotate("text", x = 18, y = 4000, label = as.character(round(mean(diffModule[,4]),5))) +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12))

ggplot(data = Freplot,aes(x=Interval,y=Frequency,group = 1)) +
  geom_line(size=0.6) + theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(aes(xintercept=8),colour="#990000",size=0.8)+
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12))

min(PCClowdiffsample)
mean(diffModule[,4])

#=======
diffModule[,6] <- abs(diffModule[,6])
parameters[,6] <- abs(parameters[,6])
PCClowdiffsample <- NULL
for(i in 1:100000){
  cat(paste0(i,"\r"))
  PCClowdiffsample <- c(PCClowdiffsample,mean(parameters[sample(1:nrow(parameters),nrow(diffModule),replace = F),6]))
}
Freplot <- data.frame(table(cut(PCClowdiffsample,
                                breaks = seq(mean(diffModule[,6])-0.01,
                                             max(PCClowdiffsample),
                                             (max(PCClowdiffsample)-mean(diffModule[,6]))/100))))

colnames(Freplot) <- c("Interval","Frequency")
ggplot(data = Freplot,aes(x=Interval,y=Frequency))+
  geom_col(width=0.9)+theme_bw()+
  geom_vline(aes(xintercept=8),colour="#990000",linetype="dashed",size=1)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(breaks=c("(0.208,0.209]","(0.248,0.249]","(0.278,0.279]","(0.303,0.304]"))+
  annotate("text", x = 18, y = 4000, label = as.character(round(mean(diffModule[,6]),5))) +
  theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"),
        axis.text.x=element_text(vjust=1,size=10,face = "bold"),
        axis.title.x = element_text(vjust=1,size=10,face = "bold"),
        axis.title.y = element_text(vjust=1,size=10,face = "bold"))

min(PCClowdiffsample)
mean(diffModule[,6])

ggplot(data = Freplot,aes(x=Interval,y=Frequency,group = 1)) +
  geom_line(size=0.6) + theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(aes(xintercept=8),colour="#990000",size=0.8)+
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12))
