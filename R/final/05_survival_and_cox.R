library(dplyr)
library(stringr)
library(ggpubr)
library(cowplot)
library(survival)
library(survminer)

source('./R/lib.R')
#################### clinical data
clinical1.raw <- read.delim('./data/tcga-prad/nationwidechildrens.org_clinical_patient_prad.txt', stringsAsFactors=FALSE)[-c(1,2),]
clin.df <- dplyr::select(clinical1.raw,
                   case=bcr_patient_barcode,
                   gleason=gleason_score,
                   new_tumor=days_to_biochemical_recurrence_first,
                   #new_tumor=new_tumor_event_dx_indicator,
                   last=last_contact_days_to)%>%
    mutate(new_tumor=suppressWarnings(as.numeric(new_tumor)),
           last=suppressWarnings(as.numeric(last)))
df<-data.frame(case=as.vector(clin.df$case),
               time=ifelse(is.na(clin.df$new_tumor), clin.df$last,clin.df$new_tumor),
               status=ifelse(is.na(clin.df$new_tumor), 1,2),
               stringsAsFactors = FALSE)

clinical2.raw <- read.delim('./data/tcga-prad/clinical/clinical.tsv', stringsAsFactors=FALSE)
clin.o <- dplyr::select(clinical2.raw,
                    case=submitter_id,
                    death=days_to_death,
                    last=days_to_last_follow_up)%>%
    mutate(death=suppressWarnings(as.numeric(death)),
           last=suppressWarnings(as.numeric(last)))

overall<-data.frame(case=as.vector(clin.o$case),
               time=ifelse(is.na(clin.o$death), clin.o$last,clin.o$death),
               status=ifelse(is.na(clin.o$death), 1,2),
               stringsAsFactors = FALSE)
#################################
diff<-pf.get.diff()

get.surv.data<-function(gene,type='os'){
  if (type=='os') {
    surv<-overall
  }else if(type=='dfs'){
    surv<-df
  }else {
    stop()
  }
  data<-pf.filter.zfpkm(gene)
  if (length(gene)==1) {
    gzfpkm<-data.frame(data)
    cases<-names(data)
  } else {
    gzfpkm <- data.frame(t(data))
    cases<-rownames(t(data))
  }
  colnames(gzfpkm)<-gene
  case.1<-as.vector(sapply(str_split(cases,'\\.'),function(x){paste(x[1:3], collapse = '-')}))
  case.2<-as.vector(sapply(str_split(cases,'\\.'),function(x){x[4]}))
  gzfpkm.2<-mutate(gzfpkm,case=case.1, flag=case.2)%>%
    filter(!grepl('11A|11B',flag))%>%dplyr::select(-flag)
  right_join(surv,gzfpkm.2,by=c('case'='case'))
}

coxph.2 <- function(gene,type=type){
  surv.data<-get.surv.data(gene,type=type)%>%dplyr::select(-case)
  model<-summary(coxph(Surv(time, status) ~ ., data=surv.data))
  ret <- c(model$logtest[3],as.vector(model$coefficients[,c(1,5)]))
  names(ret) <- NULL
  ret
}

#################################################
############# Univariate
#################################################
#---------- dfs
i<-0;total<-length(diff$GeneID)
dfs.uni<-t(sapply(diff$GeneID, function(x){
  i<<-i+1
  print(paste0(i,'/',total))
  coxph.2(x,'dfs')
}))
du<-data.frame(gene=rownames(dfs.uni),
                    p=dfs.uni[,1],
                    FDR=p.adjust(dfs.uni[,1],method = 'BH'),
                    b=dfs.uni[,2],
                    ph=exp(dfs.uni[,2]),
                    ph.p=dfs.uni[,3],
                    ph.FDR=p.adjust(dfs.uni[,3],method = 'BH'),
                    stringsAsFactors = FALSE)%>%arrange(ph.FDR)
saveRDS(du, file = './data/coxph.dfs.nui.rds')
#---------- os
i<-0;total<-length(diff$GeneID)
os.uni<-t(sapply(diff$GeneID, function(x){
  i<<-i+1
  print(paste0(i,'/',total))
  coxph.2(x,'os')
}))
ou<-data.frame(gene=rownames(os.uni),
               p=os.uni[,1],
               FDR=p.adjust(os.uni[,1],method = 'BH'),
               b=os.uni[,2],
               ph=exp(os.uni[,2]),
               ph.p=os.uni[,3],
               ph.FDR=p.adjust(os.uni[,3],method = 'BH'),
               stringsAsFactors = FALSE)%>%arrange(ph.FDR)
saveRDS(ou, file = './data/coxph.os.nui.rds')
# --------------------------------
du<-readRDS(file = './data/coxph.dfs.nui.rds')
ou<-readRDS(file = './data/coxph.os.nui.rds')
dug<-filter(du, FDR<0.01,ph.FDR<0.01)
oug<-filter(ou, p<0.01,ph.p<0.01)
saveRDS(dug, file = './data/coxph.dfs.nui.dug.rds')
saveRDS(oug, file = './data/coxph.os.nui.oug.rds')
# -------------------------thesis

dug.thesis <- data.frame(ID=dug$gene,
                        Symobl=pf.ensembl2symbol(dug$gene),
                         HR=round(dug$ph,2),
                         logFC=round(diff[match(dug$gene,diff$GeneID),'logFC'],2),
                         stringsAsFactors = FALSE)%>%arrange(desc(HR))
write.csv(dug.thesis, './reports/thesis/coxph.fds.nui.csv',row.names = FALSE)

oug.thesis <- data.frame(ID=oug$gene,
                         Symobl=pf.ensembl2symbol(oug$gene),
                         HR=round(oug$ph,2),
                         logFC=round(diff[match(oug$gene,diff$GeneID),'logFC'],2),
                         stringsAsFactors = FALSE)%>%arrange(desc(HR))
write.csv(oug.thesis, './reports/thesis/coxph.os.nui.csv',row.names = FALSE)

iug<-intersect(oug$gene,dug$gene)
iug.thesis<-data.frame(
  ID=pf.ensembl2symbol(iug),
  Type=ifelse(pf.ensembl2biotype(iug)%in%pv.pcg, 'PCG', 'lncNRA'),
  logFC=round(diff$logFC[match(iug,diff$GeneID)],3),
  OS_HR=round(oug$ph[match(iug,oug$gene)],3),
  DFS_HR=round(dug$ph[match(iug,dug$gene)],3)
)%>%arrange(desc(abs(logFC)))
write.csv(iug.thesis, './reports/thesis/coxph.iu.nui.csv',row.names = FALSE)

#---------------------- GO


#----------------------------------------------
oug.a <- data.frame(ID=oug$gene,
                         Symobl=pf.ensembl2symbol(oug$gene),
                         type=pf.ensembl2biotype(oug$gene),
                         HR=round(oug$ph,2),
                         logFC=round(diff[match(oug$gene,diff$GeneID),'logFC'],2),
                         stringsAsFactors = FALSE)%>%arrange(desc(HR))
dim(filter(oug.a, type%in%pv.lncRNA))
dim(filter(oug.a, type%in%pv.pcg))
table(data.frame(oug.a$HR>1))
table(data.frame(oug.a$logFC>0))
table(data.frame(oug.a$HR>1, oug.a$logFC>0))

dug.a <- data.frame(ID=dug$gene,
                    Symobl=pf.ensembl2symbol(dug$gene),
                    type=pf.ensembl2biotype(dug$gene),
                    HR=round(dug$ph,2),
                    logFC=round(diff[match(dug$gene,diff$GeneID),'logFC'],2),
                    stringsAsFactors = FALSE)%>%arrange(desc(HR))
dim(filter(dug.a, type%in%pv.lncRNA))
dim(filter(dug.a, type%in%pv.pcg))
table(data.frame(dug.a$HR>1))
table(data.frame(dug.a$logFC>0))
table(data.frame(dug.a$HR>1, dug.a$logFC>0))

#################################################
############# Multivariate
#################################################
genes<-oug$gene
om.multi.l<-list()
total<-length(genes)
i<-0
a<-sapply(genes,function(x){
  i<<-i+1
  print(paste0(i,'/',total))
  j<-0
  sapply(genes, function(y){
    j<<-j+1
    if (i<j) {
      om.multi.l[[paste0(x,'-',y)]]<<-c(x,y,coxph.2(c(x,y),type='os'))
    }
  })
})
omm<-data.frame(do.call(rbind,om.multi.l),row.names = NULL)

om<-data.frame(g1=as.vector(omm[,1]),
                     g2=as.vector(omm[,2]),
                     p=as.numeric(as.vector(omm[,3])),
                     FDR=p.adjust(as.numeric(as.vector(omm[,3])),
                                  method = 'BH'),
                     b1=as.numeric(as.vector(omm[,4])),
                     b2=as.numeric(as.vector(omm[,5])),
                     ph1=exp(as.numeric(as.vector(omm[,4]))),
                     ph2=exp(as.numeric(as.vector(omm[,5]))),
                     ph1.p=as.numeric(as.vector(omm[,6])),
                     ph2.p=as.numeric(as.vector(omm[,7])),
                     ph1.FDR=p.adjust(as.numeric(as.vector(omm[,6])),
                                      method = 'BH'),
                     ph2.FDR=p.adjust(as.numeric(as.vector(omm[,7])),
                                      method = 'BH'),
                     row.names = NULL,
                     stringsAsFactors = FALSE)%>%dplyr::arrange(FDR)

saveRDS(om, file = './data/coxph.os.mulit.rds')

filter(om, FDR<0.01,ph1.FDR<0.05,ph2.FDR<0.05)
#------------------
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

 #############################3
getDataSurv <- function(clinical, gene='', event='death') {
  clinical <- data.frame(clinical,
                         time=ifelse(is.na(clinical[,event]), clinical$last,clinical[,event]),
                         status=ifelse(is.na(clinical[,event]), 1,2))
  left_join(get.fpkm.anno(), clinical, by=c('case'='case'))%>%
    data.frame(fpkm=fpkm[,gene], stringsAsFactors = F) %>%
    filter(!type%in%c('11A','11B')) %>%
    distinct(case, .keep_all = T) %>%
    mutate(x.mean=ifelse(fpkm>mean(fpkm), 'high', 'low')) %>%
    mutate(x.median=ifelse(fpkm>median(fpkm), 'high', 'low'))
}
##################################333

fpkm.anno <- left_join(get.fpkm.anno(),
                       clinical1, by=c('case'='case'))
fpkm.anno <- fpkm.anno %>% mutate(group=ifelse(gleason>=7, 'high', 'low'),
                                  group.2=ifelse(type%in%c('11A','11B'),'normal','tumor'))
fpkm.anno$group[fpkm.anno$type%in%c('11A','11B')] <- 'normal'


getDataSurv <- function(clinical, gene='', event='death') {
    clinical <- data.frame(clinical,
                           time=ifelse(is.na(clinical[,event]), clinical$last,clinical[,event]),
                           status=ifelse(is.na(clinical[,event]), 1,2))
    left_join(get.fpkm.anno(), clinical, by=c('case'='case'))%>%
        data.frame(fpkm=fpkm[,gene], stringsAsFactors = F) %>%
        filter(!type%in%c('11A','11B')) %>%
        distinct(case, .keep_all = T) %>%
        mutate(x.mean=ifelse(fpkm>mean(fpkm), 'high', 'low')) %>%
        mutate(x.median=ifelse(fpkm>median(fpkm), 'high', 'low'))
}
################################################################### Plot
draw.barplot.survival <- function(gene, symbol='', outdir='./reports/candidate/') {
    final<-data.frame(fpkm.anno, fpkm=fpkm[,gene], stringsAsFactors = F)

    data.surv1 <- getDataSurv(clinical1, gene=gene, event='new_tumor')
    data.surv2 <- getDataSurv(clinical2, gene=gene, event='death')
    km.dfs.media <- survfit(Surv(time, status) ~ x.median, data=data.surv1)
    km.os.media <- survfit(Surv(time, status) ~ x.median, data=data.surv2)
    ## plot
    p1 <- ggboxplot(final,
                    x = "group", y = "fpkm",
                    order=c('normal', 'low', 'high'),
                    ylab = 'FPKM', xlab = F,
                    color = "group",
                    shape = "group",
                    palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                    add = "jitter", add.params = list(fill = "white"))+
        stat_compare_means(comparisons = list(c("low", "normal"), c("high", "normal"), c("low", "high")),
                           label = "p.signif")
    p1 <- ggpar(p1, legend.title = 'Tissue')
    p2 <- ggboxplot(final,
                    x = "group.2", y = "fpkm",
                    order=c('normal', 'tumor'),
                    ylab = 'FPKM',xlab = F,
                    color = "group.2",
                    shape = "group.2",
                    palette =c("#00AFBB", "#FC4E07"),
                    add = "jitter",
                    add.params = list(fill = "white"))+
        stat_compare_means(comparisons = list(c("tumor", "normal")), label = "p.signif")
    p2 <- ggpar(p2, legend.title = 'Tissue')

    p3 <- ggsurvplot(km.dfs.media, xlab='Disease Free Survival Probability', ylab=element_blank())$plot
    p4 <- ggsurvplot(km.os.media, xlab='Overall Survival Probability', ylab=element_blank())$plot
    p.all <- plot_grid(p1, p2, p3, p4, ncol = 2, labels=c('A', 'B', 'C', 'D'))
    pic.name=paste0(outdir, ifelse(symbol=='', gene, symbol), '.jpg')
    ggsave(p.all, file=pic.name, width = 10, height = 8)
}

load('./cache/lncTP.0.x.rda')
genes <- helper.get.candidate()
genes <- genes%>%filter(fd.p<0.05|fd.p.m<0.05|o.p<0.05|o.p.m<0.05|fd.cox.p<0.05|o.cox.p<0.05) # <========================= This is the final result!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#genes <- helper.get.lncRNA.PCG()
i<-0
lapply(split(genes, seq(nrow(genes))), function(x){
    print(i)
    i<<-i+1
    #lncTP.filterd<-filterByX(lncTP.0.3, x$GeneID)
    #print(lncTP.filterd)
    #print(x)
    draw.barplot.survival(gene = x$GeneID, symbol = x$symbol)
})->tmp
################################################ Univariate Regression
survDiffGene <- function(gene) {
    data.surv1 <- getDataSurv(clinical1, gene=gene, event='new_tumor')
    data.surv2 <- getDataSurv(clinical2, gene=gene, event='death')
    survdiff.p <- function(diff) {
        pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
    }
    fds.p.mean <- survdiff.p(survdiff(Surv(time, status) ~ x.mean, data=data.surv1))
    fds.p.median <- survdiff.p(survdiff(Surv(time, status) ~ x.median, data=data.surv1))
    os.p.mean <- survdiff.p(survdiff(Surv(time, status) ~ x.mean, data=data.surv2))
    os.p.median <- survdiff.p(survdiff(Surv(time, status) ~ x.median, data=data.surv2))

    fd.cox.coefs <- summary(coxph(Surv(time, status) ~ fpkm, data=data.surv1))$coefficients
    o.cox.coefs <- summary(coxph(Surv(time, status) ~ fpkm, data=data.surv2))$coefficients

    c(fds.p.mean, fds.p.median, os.p.mean, os.p.median,
      fd.cox.coefs[2], fd.cox.coefs[5], o.cox.coefs[2], o.cox.coefs[5])
}

genes <- helper.get.lncRNA.PCG()
i<-0
lapply(split(genes, seq(nrow(genes))), function(x){
  print(i)
  i<<-i+1
  survDiffGene(gene = x$GeneID)
})->gene.surv.p
surv.p.data <- as.data.frame(t(matrix(unlist(gene.surv.p), byrow = F, ncol= length(gene.surv.p))))
colnames(surv.p.data) <- c('fds.mean','fds.median','os.mean', 'os.median', 'fd.cox.PH', 'fd.cox.p.value', 'o.cox.PH', 'o.cox.p.value')
surv.p.data <- data.frame(gene=genes$GeneID, symbol=genes$symbol, type=genes$GeneType,
                          surv.p.data,
                          fd.km.p.adjust=p.adjust(surv.p.data$fds.median,method='BH'),
                          o.km.p.adjust=p.adjust(surv.p.data$os.median,method='BH'),
                          fd.cox.p.adjust=p.adjust(surv.p.data$fd.cox.p.value, method='BH'),
                          o.cox.p.adjust=p.adjust(surv.p.data$o.cox.p.value, method='BH'),
                          stringsAsFactors = F)

surv.p.data <- arrange(surv.p.data, fds.mean)
write.csv(surv.p.data, file = './reports/surv.p.data.csv')

### stat
fd.cox.gene<-list()
fd.cox <- surv.p.data%>%arrange(fd.cox.p.adjust)%>%filter(fd.cox.p.adjust<0.05)%>%
  select(gene=gene,symbol=symbol,type=type, PH=fd.cox.PH,p.value=fd.cox.p.value,FDR=fd.cox.p.adjust)
write.csv(fd.cox, file = './reports/fd.cox.csv')
fd.cox.gene[['lncRNA']]<-(fd.cox%>%filter(type%in%config$lncRNA))$gene
fd.cox.gene[['pcg']]<-(fd.cox%>%filter(type%in%config$PCGs))$gene

o.cox.gene<-list()
o.cox <- surv.p.data%>%arrange(o.cox.p.adjust)%>%filter(o.cox.p.adjust<0.05)%>%
  select(gene=gene,symbol=symbol,type=type, PH=o.cox.PH,p.value=o.cox.p.value,FDR=o.cox.p.adjust)
write.csv(o.cox, file = './reports/o.cox.csv')
o.cox.gene[['lncRNA']]<-(o.cox%>%filter(type%in%config$lncRNA))$gene
o.cox.gene[['pcg']]<-(o.cox%>%filter(type%in%config$PCGs))$gene

fd.km.gene<-list()
fd.km <- surv.p.data%>%filter(fd.km.p.adjust<0.05)
fd.km.gene[['lncRNA']]<-(fd.km%>%filter(type%in%config$lncRNA))$gene
fd.km.gene[['pcg']]<-(fd.km%>%filter(type%in%config$PCGs))$gene

o.km.gene<-list()
o.km <- surv.p.data%>%filter(o.km.p.adjust<0.05)
o.km.gene[['lncRNA']]<-(o.km%>%filter(type%in%config$lncRNA))$gene
o.km.gene[['pcg']]<-(o.km%>%filter(type%in%config$PCGs))$gene

surv <- list(fd.cox=fd.cox.gene,
             o.cox=o.cox.gene,
             fd.km=fd.km.gene,
             o.km=o.km.gene)

saveRDS(surv, file='./cache/candidate.surv.rds')


#################################################################################### Multivariate Regression
# surv.p.data%>%arrange(cox.p.value)%>%filter(cox.p.value <0.05)->genes
# geneList=genes$gene
# geneList=geneList[1:5]
#
# mutate(clinical, time=ifelse(is.na(new_tumor), last, new_tumor), status=ifelse(is.na(new_tumor),1,2)) -> clinical.fix
# mutate(clinical2, time=ifelse(is.na(death), last, death), status=ifelse(is.na(death),1,2)) -> clinical2.fix
# get.data.surv <- function(clini.fix) {
#   if (length(geneList)==1) {
#     tmp <- left_join(get.fpkm.anno(), clini.fix, by=c('case'='case'))%>%
#       data.frame(geneList=fpkm[,geneList], stringsAsFactors = F)
#     names(tmp)[ncol(tmp)] <- geneList
#     tmp
#   } else {
#     tmp <- left_join(get.fpkm.anno(), clini.fix, by=c('case'='case'))%>%
#       data.frame(fpkm[,geneList], stringsAsFactors = F)
#   }
#   tmp %>%distinct(case, .keep_all = T)
# }
# data.surv <- get.data.surv(clinical.fix)
# data.surv2 <- get.data.surv(clinical2.fix)
#
# fit <- coxph(as.formula(paste('Surv(time, status)~',str_c(geneList, collapse = '+'))), data=data.surv)
# summary(fit)
#
#
# ##
# gs=7
# gene <- 'ENSG00000080493'  # SLC45A2
# gene2 <- 'ENSG00000113763'  # UNC5A
# filter(diff.cor.pairs, Var1==gene, Var2==gene2)
# # SLC45A2, UNC5A
