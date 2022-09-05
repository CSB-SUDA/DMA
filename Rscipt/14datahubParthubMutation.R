library(maftools)
library(ggplot2)

CLINICALdata <- read.csv(file="data/BRCA_clinical.csv",header = T)
dateGene <- read.csv(file="data/dateHubs.csv",header = T)
partGene <- read.csv(file="data/partHubs.csv",header = T)
lamlTable <- read.delim(file="data/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz",comment.char = "#")
lamlTable$Tumor_Sample_Barcode <- substr(lamlTable$Tumor_Sample_Barcode,1,16)
lamlTable <- lamlTable[lamlTable$Tumor_Sample_Barcode%in%CLINICALdata$rowname,]
write.table(lamlTable,file="laml.maf",quote = F,row.names = F,sep="\t")

laml = read.maf(maf = "laml.maf")

lamlgene <- data.frame(laml@gene.summary)
rownames(lamlgene) <- lamlgene[,1]
lamlgene <- lamlgene[,-1]

datelaml <- lamlgene[dateGene$hubs,]
datelaml <- na.omit(datelaml)
datelaml$rowSums <- rowSums(datelaml)
datelaml <- datelaml[order(datelaml$rowSums,decreasing = T),]
datelaml <- datelaml[1:30,]
partlaml <- lamlgene[partGene$hubs,]
partlaml <- na.omit(partlaml)
partlaml$rowSums <- rowSums(partlaml)
partlaml <- partlaml[order(partlaml$rowSums,decreasing = T),]
partlaml <- partlaml[1:30,]

pdf(file="datelaml.pdf",width = 10,height = 8)
oncoplot(laml,top = 20,removeNonMutated = T,genes=rownames(datelaml),gene_mar = 6.5)
dev.off()

pdf(file="partlaml.pdf",width = 10,height = 8)
oncoplot(laml,top = 20,removeNonMutated = T,genes = rownames(partlaml))
dev.off()