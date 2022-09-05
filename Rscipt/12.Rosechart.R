library(ggplot2)
my_data <- read.csv('./data/TFC_noweight.csv',header = T, stringsAsFactors = F)

my_data <- my_data[1:30,]
my_data$Interaction <- paste0(my_data$name1,"-",my_data$name2)

my_data$RGB <- colorRampPalette(c("blue","red"))(nrow(my_data))

colour <- as.vector(t(my_data$RGB))
my_data$Interaction <- factor(my_data$Interaction,levels = my_data$Interaction)


g <- ggplot(my_data,aes(x=Interaction,y=TFCs))+
  geom_bar(aes(fill=factor(TFCs)),width=1,stat = 'identity')+
  scale_fill_manual(values = colour)+
  coord_polar(theta = 'x',start = 0,direction = 1)+
  ylim(0,180)+
  theme(legend.position = 'none',axis.text.x=element_blank())

ggsave(filename = "Rosechart.pdf",g,width = 8,height = 8)
