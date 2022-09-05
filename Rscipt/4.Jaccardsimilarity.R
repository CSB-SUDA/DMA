ecN0Cluster <- read.csv(file = "data/ecN0Clusterlg10.csv")
ecNplusCluster <- read.csv(file = "data/ecNplusClusterlg10.csv")
ecN0 <- unique(ecN0Cluster$V2)
ecNplus <- unique(ecNplusCluster$V2)


jtable <- NULL
jtableRowName <- NULL
jtableColName <- NULL
for(i in 1:length(ecN0)){
    jtableline <- NULL
    for(j in 1:length(ecNplus)){
      N0gene <- ecN0Cluster[ecN0Cluster$V2==ecN0[i],1]
      Nplusgene <- ecNplusCluster[ecNplusCluster$V2==ecNplus[j],1]
      if((length(N0gene) >=10) & (length(Nplusgene) >=10)){
        jeffonvalue <- length(intersect(N0gene,Nplusgene))/
          length(union(Nplusgene,N0gene))
        jtableline <- c(jtableline,jeffonvalue)
        if(j==1)jtableRowName <- c(jtableRowName,paste0("M",i,"(",length(N0gene),")"))
        if(i==1)jtableColName <- c(jtableColName,paste0("M",j,"(",length(Nplusgene),")"))
      } 
    }
    jtable <- rbind(jtable,jtableline) 
}
rownames(jtable) <- jtableRowName
colnames(jtable) <- jtableColName

dim(jtable)
jtableLable <- jtable
pdf(file="Jaccardsimilarity.pdf",height = 5,width = 6)
pheatmap::pheatmap(jtable,display_numbers = round(jtableLable,2),cluster_rows = T,cluster_cols = T)
dev.off()