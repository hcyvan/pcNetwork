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
tfs1 <- tf2gene2[tf=='AR'&symbol%in%unique(tf2gene2$tf)]$symbol
tfs1 #276
level1<-tf2gene2[tf=='AR'&symbol%in%tfs1]

trs1<- TRS(level1,exp=exp,cutoff = 2) # 60
######################################### AR-target target
level2<-tf2gene2[tf%in%unique(trs1$target)]
trs2<- TRS(level2,exp=exp,cutoff = 2) # 281512

saveRDS(trs2, './data/trs2.rds')
trs2 <- readRDS('./data/trs2.rds')


