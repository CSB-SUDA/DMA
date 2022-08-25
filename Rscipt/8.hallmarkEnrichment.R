library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")

data <- read.delim(file = "Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv",header = T)
data <- data[,c(4,1)]
data <- unique(data)
BgRatio <- read.csv(file = "genecode-TCGA-ensemblID-geneName.csv",header = T)
BgRatio <- unique(BgRatio[BgRatio$V2 == "protein_coding",3])

data1 <- data[data$HALLMARK == "proliferative signalling" |
                data$HALLMARK == "change of cellular energetics" |
                data$HALLMARK == "escaping programmed cell death" |
                data$HALLMARK == "genome instability and mutations" |
                data$HALLMARK == "angiogenesis" |
                data$HALLMARK == "invasion and metastasis" |
                data$HALLMARK == "tumour promoting inflammation" |
                data$HALLMARK == "cell replicative immortality" |
                data$HALLMARK == "escaping immune response to cancer" |
                data$HALLMARK == "suppression of growth",]
colnames(data1) <- c("TERM","GENE")

TERM2GENE <- data1
partHub <- read.csv(file="data/partHubs.csv",header = T)
dateHub <- read.csv(file="data/dateHubs.csv",header = T)
write.csv(TERM2GENE,file = "TERM2GENE.csv")

TERM2GENE <- rbind(TERM2GENE,data.frame(TERM = "BgRatio",GENE=BgRatio))

partekk <- enricher(gene= partHub$hubs,
                    pvalueCutoff = 1,
                    TERM2GENE = TERM2GENE,
                    qvalueCutoff = 1)
dotplot(partekk,showCategory=10,color = 'pvalue')
partKEGGData <- data.frame(partekk)
rownames(partKEGGData) <- partKEGGData[,2]
write.csv(partKEGGData,file = "partKEGGData.csv")

dateekk <- enricher(gene= dateHub$hubs,
                    pvalueCutoff = 1,
                    TERM2GENE = TERM2GENE,
                    qvalueCutoff = 1)
dotplot(dateekk,showCategory=10,color = 'pvalue')
dateKEGGData <- data.frame(dateekk)
rownames(dateKEGGData) <-dateKEGGData[,2]
write.csv(dateKEGGData,file = "dateKEGGData.csv")

partKEGGDataOurs <- partKEGGData
dateKEGGDataOurs <- dateKEGGData
write.csv(partKEGGDataOurs,file = "partGOCCData.csv",quote = F)
write.csv(dateKEGGDataOurs,file = "dateGOCCData.csv",quote = F)

allKEGGID <- union(partKEGGDataOurs[,2],dateKEGGDataOurs[,2])

partKEGGData <- partKEGGData[allKEGGID,]
rownames(partKEGGData) <- allKEGGID
partKEGGData$Description <- allKEGGID
partKEGGData$p.adjust[is.na(partKEGGData$p.adjust)] <- 1
dateKEGGData <- dateKEGGData[allKEGGID,]
rownames(dateKEGGData) <- allKEGGID
dateKEGGData$Description <- allKEGGID
dateKEGGData$p.adjust[is.na(dateKEGGData$p.adjust)] <- 1
partKEGGData <- partKEGGData[order(partKEGGData$p.adjust,decreasing = F),]

barplotdata <- data.frame(Description = c(rownames(partKEGGData),rownames(dateKEGGData)),
                          pvalue = c(partKEGGData$pvalue,dateKEGGData$pvalue),
                          group = c(rep("partHub",nrow(partKEGGData)),rep("dateHub",nrow(dateKEGGData))))
barplotdata$pvalue[is.na(barplotdata$pvalue)] <- 1

barplotdata$pvalue <- (-log10(barplotdata$pvalue))
barplotdata$pvalue <- ifelse(barplotdata$group == "dateHub", 
                             barplotdata$pvalue * -1,
                             barplotdata$pvalue * 1)

dateKEGGData <- dateKEGGData[order(dateKEGGData$pvalue),]

barplotdata$Description <- factor(barplotdata$Description,levels=rev(dateKEGGData$Description))

face = "bold"
axis.text.color = "black"
axis.text.font.size = 8
axis.title.font.size = 8
legend.text.font.size = 8
p1 <- ggplot(data = barplotdata) + geom_col(width = 0.7,aes(x = factor(Description), 
                                                            y = pvalue, 
                                                            fill = group)) + 
  scale_fill_manual(values=c(dateHub = "#BC3C28", partHub = "#0072B5"))+
  scale_y_continuous(breaks = seq(from = -10, to = 10,
                                  by = 2),
                     labels = c(seq(10, 0, -2), 
                                seq(2, 10, 2))) + 
  coord_flip() + theme_bw()+ ylim(c(-max(abs(barplotdata$pvalue)),max(abs(barplotdata$pvalue))))+
  geom_hline(aes(yintercept=log10(0.05)),size=0.7,linetype="dashed",colour = "#0072B5")+
  geom_hline(aes(yintercept=-log10(0.05)),size=0.7,linetype="dashed",colour = "#BC3C28") +
  theme(axis.text=element_text(face = face, color = axis.text.color),
        axis.text.y = element_text(face = face, color=axis.text.color,size = axis.text.font.size),
        axis.text.x = element_text(face = face, color=axis.text.color,size = axis.text.font.size),
        axis.title.x = element_text(face = face,size = axis.title.font.size),
        axis.title.y = element_text(face = face,size = axis.title.font.size),
        axis.title  = element_text(face = face,size = axis.title.font.size),
        title = element_text(size= axis.title.font.size,face=face),
        legend.text =element_text(face = face,size = legend.text.font.size),
        legend.title=element_blank())+ylab("-log10(p value)")+xlab("pathway name")
ggsave(p1,filename = "hallmark.pdf",width = 8,height = 2)    
#==========================================================
DEG <- read.csv(file = "data/BRCA_log2TPM_DEGmatrix.csv",header = T,row.names = 1)

resTable <- data.frame(gene=rep("",1000),count=rep("",1000),pvalue=rep("",1000))
for(i in 1:1000){
  cat(paste0(i,"\r"))
  tempGene <- rownames(DEG)[sample(1:nrow(DEG),size = nrow(dateHub))]
  tempkk <- enricher(gene= tempGene,
                     pvalueCutoff = 1,
                     TERM2GENE = TERM2GENE,
                     qvalueCutoff = 1)
  tempData <- data.frame(tempkk)
  if(ncol(tempData)==0){
    resTable[i,1]<-paste0(tempGene,collapse = "-")
    resTable[i,2]<-paste0(rep(0,10),collapse = "-")
    resTable[i,3]<-paste0(rep(1,10),collapse = "-")
    next
  }
  tempData <- tempData[rownames(dateKEGGData),]
  rownames(tempData) <- rownames(dateKEGGData)
  tempData$pvalue[is.na(tempData$pvalue)] <- 1
  tempData$pvalue[tempData$Count==1] <- 1
  tempData$Count[is.na(tempData$Count)] <- 0
  resTable[i,1]<-paste0(tempGene,collapse = "-")
  resTable[i,2]<-paste0(tempData$Count,collapse = "-")
  resTable[i,3]<-paste0(tempData$pvalue,collapse = "-")
}

