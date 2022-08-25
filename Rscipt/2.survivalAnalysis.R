library("survival")
library("survminer")
library("RColorBrewer")

TCGAclinical <- read.csv(file="BRCA_clinical.csv",header = T,row.names = 1)
survivalTable <- data.frame(feature = TCGAclinical$N0N,
                            OS = TCGAclinical$OS,
                            OS.time = TCGAclinical$OS.time)
survivalTable$label <- survivalTable$feature
survivalTable$feature[1:200] <- "middle"
sfit <- survfit(Surv(OS.time,OS)~label,data=survivalTable)
ggsurvplot(sfit, 
           data=survivalTable,
           pval = TRUE,
           risk.table=TRUE,
           legend.labs = c("N0", "Nplus"),
           cumevents.height = 0.2,
           ncensor.plot.height = 0.2,
           ncensor.plot = TRUE,
           ggtheme=theme_bw())

p<- ggsurvplot(sfit, 
           data=survivalTable,ggtheme =theme_classic(),pval = TRUE)
p$plot
