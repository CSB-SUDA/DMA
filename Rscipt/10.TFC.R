#' Measure the role of interactions as an information medium and Identify modules from network
#'
#' @title EM function
#' @param edges data.frame, two columns,for exmaple:gene symbol A, gene symbol B.
#' @import igraph GOSemSim
#' @export EM
#' @author Fan Wang

library(igraph)
library(GOSemSim)
TFC <- function(edges,weights=NULL) {
  #check the format of parameters
  if(ncol(edges) <2)
  {
    stop("Param 'edges' input error!
         Please input two columns!
         For example:protein-protein interactions")
  }
  options(scipen=200)
  ecount<-nrow(edges)
  #Create the network
  g = graph_from_data_frame(edges, directed=FALSE)
  nodes=as.data.frame(rbind(as.matrix(edges[,1]),as.matrix(edges[,2])))
  nodes=unique(nodes)
  #Calculate the edge betweenness
  bet<-edge_betweenness(g, e = E(g), directed = F,weights=weights)
  bet<-cbind(edges,bet)
  colnames(bet)=c("name1","name2","bet")
  #Calculate functional similarity
  hsGO <- godata('org.Hs.eg.db', keytype = "SYMBOL",
                  ont="BP", computeIC=FALSE)
  funcsim<-mgeneSim(nodes[,1], semData=hsGO, measure="Wang",
                  combine="BMA", verbose=FALSE)
  upsim <- upper.tri(funcsim)
  mergesim<-data.frame(row = rownames(funcsim)[row(funcsim)[upsim]],
                    column = rownames(funcsim)[col(funcsim)[upsim]],
                    cor =(funcsim)[upsim] )
  temp=mergesim
  temp[,1]=mergesim[,2]
  temp[,2]=mergesim[,1]
  finalsim<-rbind(mergesim,temp)
  colnames(finalsim)=c("name1","name2","cor")
  colnames(edges)=c("name1","name2")
  GO<-merge(finalsim,edges,by=c("name1","name2"),sort=F)
  #TFC of edges without functional similarity to be zero
  nasim=dplyr::anti_join(edges[,1:2], GO[,1:2],by = c("name1", "name2"))
  if(nrow(nasim)!=0){
    nasim$cor=0
    GO=rbind(GO,nasim)
  }
  #Obtain topological-functional connection score
  NET=merge(GO,bet,by=c("name1","name2"))
  scaleA<-NET$bet-min(NET$bet)
  scaleB<-max(NET$bet)-min(NET$bet)
  #Edge betweenness normalization
  NET$scaleTN <- scaleA/scaleB
  TF<-NET$scaleTN+NET$cor
  NET$TFC  <- 100*TF/(2-TF)
  TFC=NET[,1:2]
  TFC$TFCs=NET$TFC
  #Top TFC score interactions distribution
  #Zoom in to see details
  TFCplot<-TFC[order(TFC$TFCs,decreasing=T),]	
  return(TFCplot) 
}
