library("clusterProfiler")
library("enrichplot")

kegg_gmt <- read.gmt("./data/c2.cp.kegg.v7.5.1.symbols.gmt")
geneList <- read.csv(file = "./data/BRCA_log2TPM_DEG.csv",header = T,row.names = 1)
geneList <- geneList[order(geneList$logFC,decreasing = T),]

genelist <- geneList$logFC
names(genelist) <- rownames(geneList)
gsea <- GSEA(genelist,
             TERM2GENE = kegg_gmt,pAdjustMethod = "BH",pvalueCutoff = 1,eps=0) #GSEA????
gseadata <- data.frame(gsea)
gseadata$p.adjust <- p.adjust(gseadata$pvalue,method = "BH")

num <- which(rownames(gseadata)=="KEGG_CALCIUM_SIGNALING_PATHWAY")
pdf("KEGG_CALCIUM_SIGNALING_PATHWAY",width = 10,height = 8)
gseaplot2(gsea,num,title = paste0(gseadata$Description[num]),base_size = 16)
dev.off()
