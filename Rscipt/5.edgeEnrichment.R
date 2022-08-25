library(biomaRt)
library(clusterProfiler)
#======>
Rectome <- read.delim(file = "data/reactome.homo_sapiens.interactions.psi-mitab.txt",header = T)
RectomePathway <- read.delim(file = "data/ReactomePathways.txt",header = F)
Rectome <- Rectome[(grepl("uniprotkb",Rectome$X.ID.s..interactor.A) & grepl("uniprotkb",Rectome$ID.s..interactor.B)),]
Rectome[,1] <- gsub(".*:","",Rectome[,1])
Rectome[,2] <- gsub(".*:","",Rectome[,2])
for(i in 1:nrow(Rectome)){
  if(Rectome[i,1]>Rectome[i,2]){
    tempVar <- Rectome[i,1]
    Rectome[i,1] <- Rectome[i,2]
    Rectome[i,2] <- tempVar
  }
}
Rectome$RectomeEnrichmentID <- paste0(Rectome[,1],"-",Rectome[,2])
TERM2GENE <- NULL
for(i in 1:nrow(Rectome)){
  cat(paste0("dealing to get TERM2GENE ..... NO.",i,"\r"))
  tempdata <- data.frame(unlist(strsplit(Rectome[i,28],split = "\\|")),Rectome$RectomeEnrichmentID[i])
  TERM2GENE <- rbind(TERM2GENE,tempdata)
}
colnames(TERM2GENE) <- c("term","gene")
TERM2GENE <- TERM2GENE[grepl("pathway",TERM2GENE$term),]
TERM2GENE <- unique(TERM2GENE)
TERM2GENE$term <- substr(TERM2GENE$term,9,100)
TERM2NAME <- data.frame(term=RectomePathway$V1,name=RectomePathway$V2)
TERM2NAME <- TERM2NAME[grepl("R-HSA",TERM2NAME$term),]
#=======>
parameterALL <- read.csv(file="data/parameters.csv",header = T)
parameters <- parameterALL[parameterALL$N0NodeModule==parameterALL$NplusNodeModule,]
my_hgnc_symbol_id <- unique(c(parameters[,1],parameters[,2]))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
hg_symbols<- getBM(attributes=c('hgnc_symbol',"uniprotswissprot"), 
                   filters= 'hgnc_symbol', 
                   values = my_hgnc_symbol_id,
                   mart = mart)
hg_symbols <- hg_symbols[hg_symbols$uniprotswissprot!="",]
hg_symbols <- unique(hg_symbols)
parameters[,1] <- hg_symbols[match(parameters[,1],hg_symbols$hgnc_symbol),2]
parameters[,2] <- hg_symbols[match(parameters[,2],hg_symbols$hgnc_symbol),2]
parameters <- na.omit(parameters)
for(i in 1:nrow(parameters)){
  if(parameters[i,1]>parameters[i,2]){
    tempVar <- parameters[i,1]
    parameters[i,1] <- parameters[i,2]
    parameters[i,2] <- tempVar
  }
}
queryEnrichenrID <- paste0(parameters[,1],"-",parameters[,2])

p1 <- enricher(queryEnrichenrID,TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME)

#write.csv(data.frame(p1),file="intraEdge.csv",row.names = F)
dotplot(p1)


parameterALL <- read.csv(file="data/parameters.csv",header = T)
parameters <- parameterALL[parameterALL$N0NodeModule!=parameterALL$NplusNodeModule,]
my_hgnc_symbol_id <- unique(c(parameters[,1],parameters[,2]))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
as <- listAttributes(mart)
hg_symbols<- getBM(attributes=c('hgnc_symbol',"uniprotswissprot"), 
                   filters= 'hgnc_symbol', 
                   values = my_hgnc_symbol_id,
                   mart = mart)
hg_symbols <- hg_symbols[hg_symbols$uniprotswissprot!="",]
hg_symbols <- unique(hg_symbols)
parameters[,1] <- hg_symbols[match(parameters[,1],hg_symbols$hgnc_symbol),2]
parameters[,2] <- hg_symbols[match(parameters[,2],hg_symbols$hgnc_symbol),2]
parameters <- na.omit(parameters)
for(i in 1:nrow(parameters)){
  if(parameters[i,1]>parameters[i,2]){
    tempVar <- parameters[i,1]
    parameters[i,1] <- parameters[i,2]
    parameters[i,2] <- tempVar
  }
}
queryEnrichenrID <- paste0(parameters[,1],"-",parameters[,2])


p2 <- enricher(queryEnrichenrID,TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME)

write.csv(data.frame(p2),file="interEdge.csv",row.names = F)
dotplot(p2)
