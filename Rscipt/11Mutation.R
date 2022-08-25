library(maftools)
library(ggplot2)

CLINICALdata <- read.csv(file="data/BRCA_clinical.csv",header = T)
dateGene <- read.csv(file="data/dateHubs.csv",header = T)
partGene <- read.csv(file="data/partHubs.csv",header = T)

lamlTable <- read.delim(file="data/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz",comment.char = "#")
lamlTable$Tumor_Sample_Barcode <- substr(lamlTable$Tumor_Sample_Barcode,1,16)

N0laml <- lamlTable[lamlTable$Tumor_Sample_Barcode%in%CLINICALdata[CLINICALdata$N0N=="N0",1],]
Npluslaml <- lamlTable[lamlTable$Tumor_Sample_Barcode%in%CLINICALdata[CLINICALdata$N0N=="Nplus",1],]
write.table(N0laml,file="N0laml.maf",quote = F,row.names = F,sep="\t")
write.table(Npluslaml,file="Npluslaml.maf",quote = F,row.names = F,sep="\t")
#===================================================================
SELECTGene <- dateGene$hubs

lamlN0 = read.maf(maf = "N0laml.maf")
lamlplus = read.maf(maf = "Npluslaml.maf")

lamlN0gene <- data.frame(lamlN0@gene.summary)
rownames(lamlN0gene) <- lamlN0gene[,1]
lamlN0gene <- lamlN0gene[,-1]
lamlplusgene <- data.frame(lamlplus@gene.summary)
rownames(lamlplusgene) <- lamlplusgene[,1]
lamlplusgene <- lamlplusgene[,-1]
commonGene <- intersect(rownames(lamlN0gene),rownames(lamlplusgene))
SELECTGene <- intersect(commonGene,SELECTGene)
lamlN0gene <- lamlN0gene[SELECTGene,1:9]
lamlplusgene <- lamlplusgene[SELECTGene,1:9]
lamlN0gene <- lamlN0gene/rowSums(lamlN0gene)
lamlplusgene <- lamlplusgene/rowSums(lamlplusgene)
rownames(lamlN0gene) <- paste0("N0.",rownames(lamlN0gene))
rownames(lamlplusgene) <- paste0("Nplus.",rownames(lamlplusgene))

#pdf(file="summaryN0.pdf",width = 10,height = 14)
#plotmafSummary(maf = lamlN0, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#dev.off()
pdf(file="datelamlN0.pdf",width = 9,height = 12)
oncoplot(lamlN0,top = 20,removeNonMutated = T,genes = SELECTGene)
dev.off()

#pdf(file="summaryNplus.pdf",width = 8,height = 6)
#plotmafSummary(maf = lamlplus, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#dev.off()
pdf(file="datelamlNplus.pdf",width = 9,height = 12)
oncoplot(lamlplus,top = 20,removeNonMutated = T,genes = SELECTGene)
dev.off()
#=============================================
SELECTGene <- partGene$hubs

lamlN0 = read.maf(maf = "N0laml.maf")
lamlplus = read.maf(maf = "Npluslaml.maf")

lamlN0gene <- data.frame(lamlN0@gene.summary)
rownames(lamlN0gene) <- lamlN0gene[,1]
lamlN0gene <- lamlN0gene[,-1]
lamlplusgene <- data.frame(lamlplus@gene.summary)
rownames(lamlplusgene) <- lamlplusgene[,1]
lamlplusgene <- lamlplusgene[,-1]
commonGene <- intersect(rownames(lamlN0gene),rownames(lamlplusgene))
SELECTGene <- intersect(commonGene,SELECTGene)
lamlN0gene <- lamlN0gene[SELECTGene,1:9]
lamlplusgene <- lamlplusgene[SELECTGene,1:9]
lamlN0gene <- lamlN0gene/rowSums(lamlN0gene)
lamlplusgene <- lamlplusgene/rowSums(lamlplusgene)
rownames(lamlN0gene) <- paste0("N0.",rownames(lamlN0gene))
rownames(lamlplusgene) <- paste0("Nplus.",rownames(lamlplusgene))


pdf(file="partlamlN0.pdf",width = 9,height = 12)
oncoplot(lamlN0,top = 20,removeNonMutated = T,genes = SELECTGene)
dev.off()

pdf(file="partlamlNplus.pdf",width = 9,height = 12)
oncoplot(lamlplus,top = 20,removeNonMutated = T,genes = SELECTGene)
dev.off()
#==============================================================
if(T){
    mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                     axis.title = element_text(size = 12,color ="black"), 
                     axis.text = element_text(size= 12,color = "black"),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     axis.text.x = element_text(angle = 55, hjust = 1 ),
                     panel.grid=element_blank(),
                     legend.position = "top",
                     legend.text = element_text(size= 12),
                     legend.title= element_text(size= 12)
    ) }

counti <- 4

for(i in 1:counti){
    data_m <- t(rbind(lamlN0gene[i,],lamlplusgene[i,]))
    data_m <- reshape2::melt(data_m)
    p <- ggplot(data_m, aes(x=Var1, y=value))
    p <- p + geom_bar(stat="identity", position="dodge", aes(fill=Var2))+mytheme+
        xlab("variant class")+ylab("frequency")+scale_fill_brewer(palette = 'Set1')
    a<- paste0("p",i,"<-p")
    eval(parse(text = a))
}

cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = SELECTGene[1:counti])

