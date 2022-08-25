library(igraph)
library(ggplot2)
library(ggrepel)
library(biomaRt)

parameterALL <- read.csv(file="data/parameters.csv",header = T)
parameters <- parameterALL[parameterALL$N0NodeModule==parameterALL$NplusNodeModule,]
my_hgnc_symbol_id <- unique(c(parameters[,1],parameters[,2]))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
as <- listAttributes(mart)
hg_symbols<- getBM(attributes=c('hgnc_symbol',"uniprotswissprot"), 
                   filters= 'hgnc_symbol', 
                   values = my_hgnc_symbol_id,
                   mart = mart)
hg_symbols <- hg_symbols[hg_symbols$uniprotswissprot!="",]
hg_symbols <- unique(hg_symbols)

intra <- read.csv(file="data/intraEdge.csv",header = T)
intragene <- unique(unlist(strsplit(paste0(intra$geneID[1:9],collapse = "/"),split = "[-/]")))
intragene <- hg_symbols[hg_symbols$uniprotswissprot%in%intragene,1]

inter <- read.csv(file="data/interEdge.csv",header = T)
intergene <- unique(unlist(strsplit(paste0(inter$geneID[1:9],collapse = "/"),split = "[-/]")))
intergene <- hg_symbols[hg_symbols$uniprotswissprot%in%intergene,1]

dateHubs <- read.csv(file="data/dateHubs.csv",header = T)
partHubs <- read.csv(file="data/partHubs.csv",header = T)



intersect(intergene,dateHubs$hubs)
intersect(intragene,dateHubs$hubs)
intersect(intergene,partHubs$hubs)
intersect(intragene,partHubs$hubs)

df <- matrix(c(length(intersect(intergene,partHubs$hubs)), length(intersect(intergene,dateHubs$hubs)), 
               length(intersect(intragene,partHubs$hubs)), length(intersect(intragene,dateHubs$hubs))), nr = 2,
             dimnames = list(c("intergene", "intragene"), c( "partHubs", "dateHubs" )) )

fisher.test(df)

g <- graph_from_data_frame(parameterALL[,1:2])
g_nodes <- get.vertex.attribute(g)
gtable<- data.frame(rowname = g_nodes,degree = degree(g),betweenness = betweenness(g))

gtabledate <- gtable[gtable$name%in%intergene,]
gtabledate$label <- "No sig"
gtabledate$label[gtabledate$name%in%dateHubs$hubs] <- "Date hub"
gtabledate$label[gtabledate$name%in%partHubs$hubs] <- "Part hub"
gtabledate$text <- ""
gtabledate$text[gtabledate$label!="No sig"] <- gtabledate$name[gtabledate$label!="No sig"]


ggplot(gtabledate,aes(x=degree,y=betweenness,color = label))+
  geom_point(size=2.5,alpha = .6)+theme_bw()+
  scale_color_manual(values =c("#BC3C28","grey","#0072B5"))+
  geom_label_repel(data=gtabledate,aes(x=degree,y=betweenness,
    label = text),size=2,color="black",max.overlaps=20) + theme_classic()+
  theme(axis.text.y=element_text(vjust=1,size=12,face = "bold"),
        axis.text.x=element_text(vjust=1,size=12,face = "bold"),
        axis.title.x = element_text(vjust=1,size=12,face = "bold"),
        axis.title.y = element_text(vjust=1,size=12,face = "bold"))

#=============
g <- graph_from_data_frame(parameterALL[,1:2])
g_nodes <- get.vertex.attribute(g)
gtable<- data.frame(rowname = g_nodes,degree = degree(g),betweenness = betweenness(g))

gtabledate <- gtable[gtable$name%in%intragene,]
gtabledate$label <- "No sig"
gtabledate$label[gtabledate$name%in%dateHubs$hubs] <- "Date hub"
gtabledate$label[gtabledate$name%in%partHubs$hubs] <- "Part hub"
gtabledate$text <- ""
gtabledate$text[gtabledate$label!="No sig"] <- gtabledate$name[gtabledate$label!="No sig"]

ggplot(gtabledate,aes(x=degree,y=betweenness,color = label))+
  geom_point(size=2.5,alpha = .6)+ theme_classic()+
  scale_color_manual(values =c("#BC3C28","grey","#0072B5"))+
  geom_label_repel(data=gtabledate,aes(x=degree,y=betweenness,
                                       label = text),size=2,color="black",max.overlaps=20) +
  theme(axis.text.y=element_text(vjust=1,size=12,face = "bold"),
        axis.text.x=element_text(vjust=1,size=12,face = "bold"),
        axis.title.x = element_text(vjust=1,size=12,face = "bold"),
        axis.title.y = element_text(vjust=1,size=12,face = "bold"))

