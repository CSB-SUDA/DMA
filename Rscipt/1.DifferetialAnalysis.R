library(ggplot2)
library(stringr)
library(limma)
library(ggrepel)
library(GeoTcgaData)

ostime <- read.table(file="./data/TCGA-BRCA.survival.tsv",fill=T,header=T,row.names = 1)
ENSG2SYMBOL <- read.csv(file="./data/protein_coding_genecode_v22.csv",header = F)
groupData <- read.delim(file="./data/TCGA-BRCA.GDC_phenotype.tsv.gz",header=T,fill=T,row.names = 1)
HiSeqV2 <- read.table(file="./data/TCGA-BRCA.htseq_fpkm.tsv.gz",header=T,row.names = 1,check.names = F)

HiSeqV2 <- 2^HiSeqV2-1
HiSeqV2 <- fpkmToTpm_matrix(HiSeqV2)
HiSeqV2 <- log2(HiSeqV2+1)

ENSG2SYMBOL <- ENSG2SYMBOL[ENSG2SYMBOL$V2=="protein_coding",]
commonGene <- intersect(rownames(HiSeqV2),ENSG2SYMBOL$V1)
HiSeqV2 <- HiSeqV2[commonGene,]
ENSG2SYMBOL <- ENSG2SYMBOL[match(commonGene,ENSG2SYMBOL$V1),]
rownames(HiSeqV2) <- ENSG2SYMBOL$V3
HiSeqV2 <- aggregate(HiSeqV2,list(ENSG2SYMBOL$V3),mean)
rownames(HiSeqV2) <- HiSeqV2[,1]
HiSeqV2 <- HiSeqV2[,-1]
#==================================================

HiSeqV2 <- HiSeqV2[,substr(colnames(HiSeqV2),14,15)!=11]
commonSample <- intersect(colnames(HiSeqV2),rownames(groupData))
commonSample <- intersect(commonSample,rownames(ostime))
groupData <- groupData[commonSample,]
ostime <- ostime[commonSample,]
HiSeqV2 <- HiSeqV2[,commonSample]
#==================================================

groupData$pathologic_N <- substr(groupData$pathologic_N,1,2)
grouplist <- groupData$pathologic_N
HiSeqV2 <- HiSeqV2[,((grouplist!="") & (grouplist!="NX" ))]
groupData <- groupData[((grouplist!="") & (grouplist!="NX" )),]
groupData$pathologic_N[groupData$pathologic_N!="N0"]="Nplus"
table(groupData$pathologic_N)
ostime <- ostime[rownames(groupData),]
HiSeqV2 <- HiSeqV2[,rownames(groupData)]
#==================================================

m.mad <- apply(HiSeqV2,1,mad)
HiSeqV2 <- HiSeqV2[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
#==================================================

group_list <- as.factor(groupData$pathologic_N)

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(HiSeqV2)

contrast.matrix<-makeContrasts("Nplus-N0",levels = design)

fit <- lmFit(HiSeqV2,design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)

DESeq2_DEG = na.omit(tempOutput)
#=====================================================
write.csv(DESeq2_DEG,file=paste0("BRCA_log2TPM_DEG.csv"),quote = F,row.names = T)
write.csv(HiSeqV2,file=paste0("BRCA_log2TPM_matrix.csv"),quote = F,row.names = T)

nrDEG=DESeq2_DEG[,c(1,4)]
colnames(nrDEG)=c('logFoldChange','pvalue')
nrDEG$pvalue <- as.numeric(nrDEG$pvalue)
nrDEG$logFoldChange <- as.numeric(nrDEG$logFoldChange)
nrDEG$label <- "nosig"
nrDEG$label[nrDEG$pvalue<0.01 & nrDEG$logFoldChange>0] <- "up"
nrDEG$label[(nrDEG$pvalue<0.01) & (nrDEG$logFoldChange < 0)] <- "down"
DEGs <- rownames(nrDEG[nrDEG$pvalue < 0.01,])
nrDEG$pvalue <- -log10(nrDEG$pvalue)

Volcano <- function(nrDEG){
    nrDEG <- nrDEG[order(nrDEG$pvalue,decreasing = T),]
    nrDEG$text <- rownames(nrDEG)
    nrDEG$text[11:nrow(nrDEG)]<- ""
    ggplot(nrDEG,aes(x=logFoldChange,y=pvalue,color = label))+
        geom_point(size=1.3,alpha = .6)+
        scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
        geom_hline(yintercept=-log10(0.01),linetype=4,size=0.6)+
        geom_vline(xintercept=0,linetype=4,size=0.6)+
        theme(panel.grid =element_blank())+theme_classic()+
        xlim(-max(nrDEG$logFoldChange),max(nrDEG$logFoldChange))+xlab("log2(fold change)")+
        ylab("-log10(p value)") + geom_label_repel(data=nrDEG,aes(
            x=logFoldChange,y=pvalue,
            label = text), color="black",max.overlaps=100)
}
p <- Volcano(nrDEG)
ggsave(p,filename = "volo.pdf",width = 8,height = 6)

write.csv(HiSeqV2[DEGs,],file=paste0("BRCA_log2TPM_DEGmatrix.csv"),quote = F,row.names = T)
groupData <- data.frame(rowname=rownames(groupData),N0N=groupData[,"pathologic_N"],
                        OS = ostime$OS,OS.time = ostime$OS.time)
groupData$OS.time <- groupData$OS.time/365
write.csv(groupData,file=paste0("BRCA_clinical.csv"),quote = F,row.names = F)

