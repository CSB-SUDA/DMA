checkEdgeIn <- function(network = interModuleEdge,ouredge = interModuleEdgei){
    for(i in 1:nrow(network)){
        if(network[i,1] > network[i,2]){
            tempedge <- network[i,1]
            network[i,1] <- network[i,2]
            network[i,2] <- tempedge
        }
    }
    for(i in 1:nrow(ouredge)){
        if(ouredge[i,1] > ouredge[i,2]){
            tempedge <- ouredge[i,1]
            ouredge[i,1] <- ouredge[i,2]
            ouredge[i,2] <- tempedge
        }
    }
    network$Label <- paste0(network$node1,network$node2)
    ouredge$Label <- paste0(ouredge$node1,ouredge$node2)
    return(network$Label%in%ouredge$Label)
}
#===========================================================
library("igraph")
library("ggplot2")
sampleTime <- 200

stringdb <- read.csv(file="data/parameters.csv",header = T)
interModuleEdge <- stringdb[(stringdb$N0NodeModule != stringdb$NplusNodeModule),]
intraModuleEdge <- stringdb[(stringdb$N0NodeModule == stringdb$NplusNodeModule),]

N0lineinterTable <- NULL
N0lineintraTabel <- NULL
for(i in 1:(round(nrow(interModuleEdge)/20))){
    cat(paste0("\n",i,"<=======>\n"))
    lineinter <- NULL
    lineintra <- NULL
    for(j in 1:sampleTime){
        cat(paste0(j,"\r"))
        sampleV <- sample(1:nrow(interModuleEdge),size = i*20,replace = F)
        interModuleEdgei <- interModuleEdge[sampleV,]
        sampleV <- sample(1:nrow(intraModuleEdge),size = i*20,replace = F)
        intraModuleEdgei <- intraModuleEdge[sampleV,]
        
        intertempdb <- stringdb[!checkEdgeIn(stringdb,interModuleEdgei),]
        intratempdb <- stringdb[!checkEdgeIn(stringdb,intraModuleEdgei),]
        interg <- graph_from_data_frame(intertempdb, directed=F)
        intrag <- graph_from_data_frame(intratempdb, directed=F)
        lineinter <- c(lineinter,mean(betweenness(interg,weights = abs(intertempdb$N0cor))))
        lineintra <- c(lineintra,mean(betweenness(intrag,weights = abs(intratempdb$N0cor))))
    }
    N0lineinterTable <- rbind(N0lineinterTable,lineinter)
    N0lineintraTabel <- rbind(N0lineintraTabel,lineintra)
    cat(paste0(mean(lineinter),"<===>",mean(lineintra)))
}
colnames(N0lineintraTabel) <- paste0("sample",1:sampleTime)
rownames(N0lineintraTabel) <- paste0(1:(round(nrow(interModuleEdge)/20)))
colnames(N0lineinterTable) <- paste0("sample",1:sampleTime)
rownames(N0lineinterTable) <- paste0(1:(round(nrow(interModuleEdge)/20)))
write.csv(N0lineintraTabel,file="N0betweennessintraTable.csv",row.names = F)
write.csv(N0lineinterTable,file="N0betweennessinterTable.csv",row.names = F)

lineintraTable <- read.csv(file="N0betweennessintraTable.csv",header = T)
lineinterTabel <- read.csv(file="N0betweennessinterTable.csv",header = T)
betweTable <- rbind(cbind(rowMeans(lineintraTable),rowMeans(lineinterTabel)))
colnames(betweTable) <- c("intramodule edge","intermodule edge")
rownames(betweTable) <- 1:nrow(betweTable)
betweTable <- reshape2::melt(betweTable)
colnames(betweTable)[2] <- "hubTypes"
betweTable$hubTypes <- factor(betweTable$hubTypes,levels = c("intermodule edge","intramodule edge"))
sdymin <- c(rowMeans(lineintraTable) - apply(lineintraTable,1,sd)/sqrt(ncol(lineintraTable)),rowMeans(lineinterTabel) - apply(lineinterTabel,1,sd)/sqrt(ncol(lineinterTabel)))
sdymax <- c(rowMeans(lineintraTable) + apply(lineintraTable,1,sd)/sqrt(ncol(lineintraTable)),rowMeans(lineinterTabel) + apply(lineinterTabel,1,sd)/sqrt(ncol(lineinterTabel)))

p<- ggplot(data = betweTable,aes(x = Var1,y=value,group = hubTypes,colour=hubTypes))+
  geom_line(size=1)+ geom_errorbar(aes(ymin = sdymin,ymax = sdymax, width=1))+
  geom_point(shape=16)+scale_color_manual(values = c('#BC3C28','#66AAD2'))+
  theme_classic()+xlab("")+
  ylab("")+
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x =element_text(size=14), axis.text.y=element_text(size=14))
ggsave(p,filename = "betweennessN0edge1.pdf",width = 6,height = 4)
#=========================================================
#closeness
library("igraph")
library("ggplot2")
sampleTime <- 200

stringdb <- read.csv(file="data/parameters.csv",header = T)
interModuleEdge <- stringdb[(stringdb$N0NodeModule != stringdb$NplusNodeModule),]
intraModuleEdge <- stringdb[(stringdb$N0NodeModule == stringdb$NplusNodeModule),]

N0lineinterTable <- NULL
N0lineintraTabel <- NULL
for(i in 1:(round(nrow(interModuleEdge)/20))){
  cat(paste0("\n",i,"<=======>\n"))
  lineinter <- NULL
  lineintra <- NULL
  for(j in 1:sampleTime){
    cat(paste0(j,"\r"))
    sampleV <- sample(1:nrow(interModuleEdge),size = i*20,replace = F)
    interModuleEdgei <- interModuleEdge[sampleV,]
    sampleV <- sample(1:nrow(intraModuleEdge),size = i*20,replace = F)
    intraModuleEdgei <- intraModuleEdge[sampleV,]
    
    intertempdb <- stringdb[!checkEdgeIn(stringdb,interModuleEdgei),]
    intratempdb <- stringdb[!checkEdgeIn(stringdb,intraModuleEdgei),]
    interg <- graph_from_data_frame(intertempdb, directed=F)
    intrag <- graph_from_data_frame(intratempdb, directed=F)
    lineinter <- c(lineinter,mean(closeness(interg,weights = abs(intertempdb$N0cor),normalized = T)))
    lineintra <- c(lineintra,mean(closeness(intrag,weights = abs(intratempdb$N0cor),normalized = T)))
  }
  N0lineinterTable <- rbind(N0lineinterTable,lineinter)
  N0lineintraTabel <- rbind(N0lineintraTabel,lineintra)
  cat(paste0(mean(lineinter),"<===>",mean(lineintra)))
}
colnames(N0lineintraTabel) <- paste0("sample",1:sampleTime)
rownames(N0lineintraTabel) <- paste0(1:(round(nrow(interModuleEdge)/20)))
colnames(N0lineinterTable) <- paste0("sample",1:sampleTime)
rownames(N0lineinterTable) <- paste0(1:(round(nrow(interModuleEdge)/20)))
write.csv(N0lineintraTabel,file="N0closenessintraTable.csv",row.names = F)
write.csv(N0lineinterTable,file="N0closenessinterTable.csv",row.names = F)

lineintraTable <- read.csv(file="N0closenessintraTable.csv",header = T)
lineintraTable <- 1/lineintraTable
lineinterTabel <- read.csv(file="N0closenessinterTable.csv",header = T)
lineinterTabel <- 1/lineinterTabel
betweTable <- rbind(cbind(rowMeans(lineintraTable),rowMeans(lineinterTabel)))
colnames(betweTable) <- c("intramodule edge","intermodule edge")
rownames(betweTable) <- 1:nrow(betweTable)
betweTable <- reshape2::melt(betweTable)
colnames(betweTable)[2] <- "hubTypes"
betweTable$hubTypes <- factor(betweTable$hubTypes,levels = c("intermodule edge","intramodule edge"))
sdymin <- c(rowMeans(lineintraTable) - apply(lineintraTable,1,sd)/sqrt(ncol(lineintraTable)),rowMeans(lineinterTabel) - apply(lineinterTabel,1,sd)/sqrt(ncol(lineinterTabel)))
sdymax <- c(rowMeans(lineintraTable) + apply(lineintraTable,1,sd)/sqrt(ncol(lineintraTable)),rowMeans(lineinterTabel) + apply(lineinterTabel,1,sd)/sqrt(ncol(lineinterTabel)))


p<- ggplot(data = betweTable,aes(x = Var1,y=value,group = hubTypes,colour=hubTypes))+
  geom_line(size=1)+ geom_errorbar(aes(ymin = sdymin,ymax = sdymax, width=1))+
  geom_point(shape=16)+scale_color_manual(values = c('#BC3C28','#66AAD2'))+
  theme_bw()+xlab("")+
  ylab("")+theme_classic()+
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x =element_text(size=14), axis.text.y=element_text(size=14))

ggsave(p,filename = "shortestPathWayN0edge1.pdf",width = 6,height = 4)

