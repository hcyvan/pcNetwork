library(ggpubr)
library(cowplot)
library(dplyr)

#############################
#############################
draw <- function(data,order,color,comparison){
  ggboxplot(data,x = "group", y = "value",order=order,
            ylab = 'Protein Expression Level', xlab = '',
            color = "group",shape = "group",palette =color,
            add = "jitter", add.params = list(fill = "white"),ggtheme = theme_pubr())+
    theme(axis.text=element_text(size=rel(1.2)),
          axis.text.x = element_text(size=rel(0),angle = 45),
          legend.text= element_text(size=rel(1)),
          legend.title=element_blank(),
          axis.title=element_text(size=rel(1.2)))+
    stat_compare_means(comparisons = comparison,label = "p.signif")
}
############################

data <- read.delim('./reports/thesis/origin/77829.csv',sep = ',', stringsAsFactors = FALSE)
colnames(data)<-c('id','tissue','gleason', 'value')
data <- mutate(data, gleason=as.integer(gleason), value=as.numeric(value))%>%
  filter(!is.na(value))%>%
  mutate(group2=ifelse(tissue=='癌',ifelse(gleason>7, 'Gleason>7','Gleason<=7'),'ParaCancer'))%>%
  mutate(group1=ifelse(tissue=='癌','Cancer','ParaCancer'))

data<-mutate(data,group=group2)
p1<-draw(data,order=c('ParaCancer', 'Gleason<=7', 'Gleason>7'),
     color=c("#00AFBB", "#E7B800", "#FC4E07"),
     comparison = list(c("Gleason<=7", "ParaCancer"),
                       c("Gleason>7", "ParaCancer"),
                       c("Gleason>7", "Gleason<=7")))
data<-mutate(data,group=group1)
p2<-draw(data,order=c('ParaCancer', 'Cancer'),
         color = c("#00AFBB", "#FC4E07"),
         comparison = list(c("Cancer", "ParaCancer")))

plot_grid(p2,p1)
############################################################
data <- read.delim('./reports/thesis/origin/78561.csv',sep = ',', stringsAsFactors = FALSE)
colnames(data)<-c('id','tissue','gleason', 'value')
data <- mutate(data, gleason=as.integer(gleason), value=as.numeric(value))%>%
  filter(!is.na(value))%>%
  mutate(group2=ifelse(tissue=='癌',ifelse(gleason>7, 'Gleason>7','Gleason<=7'),'ParaCancer'))%>%
  mutate(group1=ifelse(tissue=='癌','Cancer','ParaCancer'))

data<-mutate(data,group=group2)
p1<-draw(data,order=c('ParaCancer', 'Gleason<=7', 'Gleason>7'),
         color=c("#00AFBB", "#E7B800", "#FC4E07"),
         comparison = list(c("Gleason<=7", "ParaCancer"),
                           c("Gleason>7", "ParaCancer"),
                           c("Gleason>7", "Gleason<=7")))
data<-mutate(data,group=group1)
p2<-draw(data,order=c('ParaCancer', 'Cancer'),
         color = c("#00AFBB", "#FC4E07"),
         comparison = list(c("Cancer", "ParaCancer")))

plot_grid(p2,p1)