library(pcProfile)
library(data.table)
library(dplyr)
source('./R/lib.R')

# The Diff Score is a transformation of the p-value that provides directionality to the p-value based on the difference between the average signal in the reference group vs. the comparison group. The formula is: DiffScore = 10*sgn(µcond-µref)*log10p; For a p-value of 0.05, DiffScore = ± 13; For a p-value of 0.01, DiffScore = ± 22; For a p-value of 0.001, DiffScore = ± 33 The p-value column is hidden by default. To display this column, use the Column Chooser.
###################################
E <- function(a,b) {
  a<-log(a[2:10]/a[1],2)
  b<-log(b[2:10]/b[1],2)
  z<- c(1/3,2/3,1,2,4,8,16,24,48)
  mw <- function(x,y,cutoff=log(1.5,2)) {
    ifelse(abs(x)>cutoff, sign(x),0)*ifelse(abs(y)>cutoff, sign(y),0)
  }
  g<-abs(mw(a,b))
  h<-abs(mw(a[1:8],b[2:9]))/(z[2:9]-z[1:8])
  p<-ifelse(mw(a[1:8],b[1:8])*mw(a[2:9],b[2:9])==-1,-1,0)
  g[is.na(g)]<-0
  h[is.na(h)]<-0
  p[is.na(p)]<-0
  sum(1,g,h,p)
}

TRS<-function(pairs,exp,cutoff=1.28){
  total<-length(unique(as.vector(pairs$tf)))
  i<-0
  tft<-do.call(rbind,lapply(split(pairs,as.vector(pairs$tf)), function(sub){
    i<<-i+1
    tf <- sub$tf[1]
    cat(paste0('process: ',tf,'    ',i,'/',total,'\n'))
    g1<-unlist(exp[symbol==tf][,-c(1,2)])
    a<-do.call(rbind,lapply(split(sub, seq(nrow(sub))), function(x){
      g2<-unlist(exp[symbol==x$symbol][,-c(1,2)])
      if (sum(g1)==0|sum(g2)==0){
        NA
      }else{
        test <- cor.test(g1,g2)
        e<-E(g1,g2)
        c(x$tf,x$symbol,e,test$estimate,test$p.value)
      }
    }))
    if(all(is.na(a))) {
      NA
    } else {
      a
    }
  }))
  tft<-as.data.frame(tft)
  colnames(tft)<-c('tf','target','e','r','p.value')
  ret<-tft%>%filter(!is.na(tf))%>%
    mutate(e=as.numeric(as.vector(e)),
           r=as.numeric(as.vector(r)),
           p.value=as.numeric(as.vector(p.value)))%>%
    mutate(FDR=p.adjust(p.value),
           trs=e*log((1+r)/(1-r))/2)
  filter(ret,abs(trs)>cutoff)
}
######################################
data("tf2gene.gtrd")
tf2gene<-as.data.table(tf2gene.gtrd)
tf2gene[,symbol:=pf.ensembl2symbol(tf2gene$gene)]

tf2gene2<-do.call(rbind,lapply(split(tf2gene, as.vector(tf2gene$tf)),function(x){
  x[N>=quantile(x$N,0.85)&symbol!=""]
}))

exp<-as.data.table(readRDS('./data/beachip.rds')$value)
########################################## AR target
#-----------Gene
tfs1.gene <- tf2gene2[tf=='AR']$symbol
length(tfs1.gene) #5083
level1.gene<-tf2gene2[tf=='AR'&symbol%in%tfs1.gene]
trs1gene<- TRS(level1.gene,exp=exp,cutoff = 2) # 60
dim(trs1gene)
lapply(dge, function(x){
  length(intersect(trs1gene$target,x))
})
#-----------TF
tfs1 <- tf2gene2[tf=='AR'&symbol%in%unique(tf2gene2$tf)]$symbol
tfs1 #276
level1<-tf2gene2[tf=='AR'&symbol%in%tfs1]

trs1<- TRS(level1,exp=exp,cutoff = 2) # 60
trs1.export<-data.frame(TF=trs1$target,TRS=round(trs1$trs,2))%>%arrange(desc(abs(TRS)))
write.csv(trs1.export,'reports/thesis/trs1.export.csv',row.names = FALSE)
######################################### AR-target target
level2<-tf2gene2[tf%in%unique(trs1$target)]
# trs2<- TRS(level2,exp=exp,cutoff = 2) # 281512
# 
# saveRDS(trs2, './data/trs2.rds')
trs2 <- readRDS('./data/trs2.rds')
diff.ar<-pf.get.diff.ar()
early<-diff.ar$early
late<-diff.ar$late
inter<-diff.ar$inter
early.only<-setdiff(early,inter)
late.only<-setdiff(late,inter)
length(early.only)
length(late.only)
dim(filter(trs2,target%in%early.only))
dim(filter(trs2,target%in%late.only))
#---------------Gene
trs2gene<-unique(as.vector(trs2$target))
lapply(dge, function(x){
  length(intersect(trs2gene,x))
})
lapply(dge, function(x){
  length(x)
})

#-------------plot
AR1<-sapply(dge, function(x){
  length(intersect(trs1gene$target,x))
})[-1]
AR2<-sapply(dge, function(x){
  length(intersect(trs2gene,x))
})[-1]
all<-sapply(dge, function(x){
  length(x)
})[-1]
labels<-c('20min', '40min', '1h', '2h', '4h', '8h', '16h', '24h', '48h')
time.labels=factor(labels, levels=labels)
df.ar1<-data.frame(group='AR',time=time.labels,data=AR1)
df.ar2<-data.frame(group='Secondary',time=time.labels,data=AR2)
df.all<-data.frame(group='ALL',time=time.labels,data=all)
df <- rbind(df.ar1,df.ar2,df.all)
ggplot(data=df, aes(x=time, y=data, group=group)) +
  geom_line(aes(linetype=group))+
  geom_point(aes(shape=group))+
  labs(x=NULL,y='Differential gene number')+
  theme_classic(base_size = 20)+scale_color_grey()

###############################3
trs2.early.only<-filter(trs2,target%in%early.only)
trs2.early.n<-sapply(split(trs2.early.only,as.vector(trs2.early.only$tf)),function(x){nrow(x)})
trs2.late.only<-filter(trs2,target%in%late.only)
trs2.late.n<-sapply(split(trs2.late.only,as.vector(trs2.late.only$tf)),function(x){nrow(x)})
mean(trs2.late.n)
mean(trs2.early.n)
t.test(trs2.early.n,trs2.late.n)

tf2<-as.vector(trs1.export$TF)
tf2.early<-trs2.early.n[match(tf2,names(trs2.early.n))]
tf2.late<-trs2.late.n[match(tf2,names(trs2.late.n))]
tf2.early[is.na(tf2.early)]<-0
tf2.late[is.na(tf2.late)]<-0
tf.el<-data.frame(TF=tf2,Early=tf2.early,Late=tf2.late)
write.csv(tf.el,'reports/thesis/ar.tf2.early.late.csv',row.names = FALSE)
#############################
trs2.early.only<-filter(trs2,target%in%early.only)

a<-rbind(data.frame(tf=trs1$tf,target=trs1$target,level='1'),
         data.frame(tf=trs2.early.only$tf,target=trs2.early.only$target,level='2'))
nodes<-unique(c(as.vector(trs1$target),as.vector(trs2.early.only$tf)))
b<-data.frame(node=nodes,
              type='tf',
              stringsAsFactors=TRUE)
b<-rbind(b,data.frame(node='AR',type='tf1'))
write.csv(a, 'reports/thesis/cyto/ar.early.edge.csv',row.names = FALSE,quote=FALSE)
write.csv(b, 'reports/thesis/cyto/ar.early.node.csv',row.names = FALSE,quote = FALSE)

sort(sapply(split(a,as.vector(a$tf)),function(x){nrow(x)}),decreasing=TRUE)
