library("survival")
library("survminer")
library("RColorBrewer")

TCGAclinical <- read.csv(file="./data/BRCA_clinical.csv",header = T,row.names = 1)
survivalTable <- data.frame(feature = TCGAclinical$N0N,
                            OS = TCGAclinical$OS,
                            OS.time = TCGAclinical$OS.time)
survivalTable$label <- survivalTable$feature
sfit <- survfit(Surv(OS.time,OS)~label,data=survivalTable)

p<- ggsurvplot(sfit, 
           data=survivalTable,
           ggtheme =theme_classic(),
           palette = c("#0070B4","#BA3521"),
           pval = TRUE)

ggsave(p$plot,filename = "survival.pdf",width = 4.5,height = 4)
