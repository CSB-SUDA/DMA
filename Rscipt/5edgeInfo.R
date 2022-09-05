library("survival")
library("survminer")
library("poolr")


stringdb <- read.csv(file="./data/stringdb.csv",header = T)
ecNplus <- read.csv(file="data/ecNplusClusterlg10.csv",header = T)
ecN0 <- read.csv(file="data/ecN0Clusterlg10.csv",header = T)
DGEmatrix <-read.csv(file="data/BRCA_log2TPM_DEGmatrix.csv",header = T,check.names = F,row.names = 1)
clinicaldata <- read.csv(file="data/BRCA_clinical.csv",header = T,row.names = 1)
clinicaldataN0 <- clinicaldata[clinicaldata$N0N=="N0",]
clinicaldataNplus <- clinicaldata[clinicaldata$N0N=="Nplus",]
DGEmatrixN0 <- DGEmatrix[,rownames(clinicaldataN0)]
DGEmatrixNplus <- DGEmatrix[,rownames(clinicaldataNplus)]
#========================================================================
for(linei in 1:nrow(stringdb)){
    cat(paste0(linei,"\r"))
    node1 <- stringdb[linei,1]
    node2 <- stringdb[linei,2]
    node1 <- as.numeric(DGEmatrixN0[node1,])
    node2 <- as.numeric(DGEmatrixN0[node2,])
    if(mean(node1)>mean(node2)){
        tempnode <- node1
        node1 <- node2
        node2 <- tempnode
    }
    clinicaldataN0$score <- (node2+1)/(node1+1)
    clinicaldataN0$label <- "low"
    clinicaldataN0$label[clinicaldataN0$score>median(clinicaldataN0$score)] <- "high"
    sdf<-survdiff(Surv(OS.time,OS)~label,data=clinicaldataN0)
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    stringdb$N0pvalue[linei] <- p.val
}
for(linei in 1:nrow(stringdb)){
    cat(paste0(linei,"\r"))
    node1 <- stringdb[linei,1]
    node2 <- stringdb[linei,2]
    node1 <- as.numeric(DGEmatrixNplus[node1,])
    node2 <- as.numeric(DGEmatrixNplus[node2,])
    if(mean(node1)>mean(node2)){
        tempnode <- node1
        node1 <- node2
        node2 <- tempnode
    }
    clinicaldataNplus$score <- (node2+1)/(node1+1)
    clinicaldataNplus$label <- "low"
    clinicaldataNplus$label[clinicaldataNplus$score>median(clinicaldataNplus$score)] <- "high"
    sdf<-survdiff(Surv(OS.time,OS)~label,data=clinicaldataNplus)
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    stringdb$Npluspvalue[linei] <- p.val
}
for(linei in 1:nrow(stringdb)){
    cat(paste0(linei,"\r"))
    node1 <- stringdb[linei,1]
    node2 <- stringdb[linei,2]
    node1 <- as.numeric(DGEmatrix[node1,])
    node2 <- as.numeric(DGEmatrix[node2,])
    if(mean(node1)>mean(node2)){
        tempnode <- node1
        node1 <- node2
        node2 <- tempnode
    }
    clinicaldata$score <- (node2+1)/(node1+1)
    clinicaldata$label <- "low"
    clinicaldata$label[clinicaldata$score>median(clinicaldata$score)] <- "high"
    sdf<-survdiff(Surv(OS.time,OS)~label,data=clinicaldata)
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    stringdb$allpvalue[linei] <- p.val
}

#============================================================
for(linei in 1:nrow(stringdb)){
    cat(paste0(linei,"\r"))
    node1 <- stringdb[linei,1]
    node2 <- stringdb[linei,2]
    stringdb$N0NodeModule[linei] <- ecN0[match(node1,ecN0[,1]),2]==ecN0[match(node2,ecN0[,1]),2]
    stringdb$NplusNodeModule[linei] <- ecNplus[match(node1,ecNplus[,1]),2]==ecNplus[match(node2,ecNplus[,1]),2]
}
#===========================================================
write.csv(stringdb,file="parameters.csv",quote = F,row.names = F)
