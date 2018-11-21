library(dplyr)
library(ggpubr)
library(cowplot)

source('./lib/globals.R')
source('./lib/helpers.R')


clinical.raw <- read.delim('./data/tcga-prad/nationwidechildrens.org_clinical_patient_prad.txt', stringsAsFactors = F)[-c(1,2),]
clinical <- select(clinical.raw, case=bcr_patient_barcode, gleason=gleason_score,new_tumor=new_tumor_event_dx_indicator, last=last_contact_days_to)%>%
            mutate(new_tumor=suppressWarnings(as.numeric(new_tumor)), last=suppressWarnings(as.numeric(last)))
rownames(clinical) <- NULL

clinical2.raw <- read.delim('./data/tcga-prad/clinical/clinical.tsv', stringsAsFactors = F)
clinical2 <- select(clinical2.raw, case=submitter_id, death=days_to_death, last=days_to_last_follow_up)%>%
            mutate(death=suppressWarnings(as.numeric(death)), last=suppressWarnings(as.numeric(last)))
rownames(clinical2) <- NULL


# data.fpkm <- helper.get.fpkm.count()
saveRDS(data.fpkm, './cache/data.fpkm.rds')
data.fpkm <- readRDS('./cache/data.fpkm.rds')
fpkm <- t(data.fpkm)

####################3
get.fpkm.anno <- function() {
  case.matrix <- str_split(rownames(fpkm), '\\.', simplify = T)
  data.frame(case=apply(case.matrix[,1:3], 1, function(x){str_c(x, collapse = '-')}), type=case.matrix[,4], stringsAsFactors = F)
}

################################################################### Plot
draw.barplot.survival <- function(gene, gs=7, symbol='') {
  ################################################################### Bar Plot
  fpkm.anno <- left_join(get.fpkm.anno(), clinical, by=c('case'='case'))
  fpkm.anno <- fpkm.anno %>% mutate(group=ifelse(gleason>=gs, 'high', 'low'), group.2=ifelse(type%in%c('11A','11B'),'normal','tumor'))
  fpkm.anno$group[fpkm.anno$type%in%c('11A','11B')] <- 'normal'
  
  final<-data.frame(fpkm.anno, fpkm=fpkm[,gene], stringsAsFactors = F)
  
  
  p1 <- ggboxplot(final,
                  x = "group", y = "fpkm",
                  order=c('normal', 'low', 'high'),
                  ylab = 'FPKM', xlab = F,
                  color = "group",
                  shape = "group",
                  palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                  add = "jitter", add.params = list(fill = "white"))+
    stat_compare_means(comparisons = list(c("low", "normal"), c("high", "normal"), c("low", "high")), label = "p.signif")
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
  ############################################################################# survival curve
  
  mutate(clinical, time=ifelse(is.na(new_tumor), last, new_tumor), status=ifelse(is.na(new_tumor),1,2)) -> clinical.fix
  mutate(clinical2, time=ifelse(is.na(death), last, death), status=ifelse(is.na(death),1,2)) -> clinical2.fix
  
  get.data.surv <- function(clini.fix) {
    left_join(get.fpkm.anno(), clini.fix, by=c('case'='case'))%>%
      data.frame(fpkm=fpkm[,gene], stringsAsFactors = F) %>%
      filter(!type%in%c('11A','11B')) %>%
      distinct(case, .keep_all = T) %>%
      mutate(x.mean=ifelse(fpkm>mean(fpkm), 'high', 'low')) %>%
      mutate(x.median=ifelse(fpkm>median(fpkm), 'high', 'low'))
  }
  data.surv <- get.data.surv(clinical.fix)
  data.surv2 <- get.data.surv(clinical2.fix)
  
  
  km.dfs.mean <- survfit(Surv(time, status) ~ x.mean, data=data.surv)
  km.dfs.media <- survfit(Surv(time, status) ~ x.median, data=data.surv)
  
  km.os.mean <- survfit(Surv(time, status) ~ x.mean, data=data.surv2)
  km.os.media <- survfit(Surv(time, status) ~ x.median, data=data.surv2)
  
  p3 <- ggsurvplot(km.dfs.media, xlab='Disease Free Survival Probability', ylab=element_blank())$plot
  p4 <- ggsurvplot(km.os.media, xlab='Overall Survival Probability', ylab=element_blank())$plot
  
  p.all <- plot_grid(p1, p2, p3, p4, ncol = 2, labels=c('A', 'B', 'C', 'D'))
  
  pic.name=paste0('./reports/survplot/',ifelse(symbol=='', gene, symbol), '.jpg')
  ggsave(p.all, file=pic.name, width = 10, height = 8)
}

# genes <- helper.get.lncRNA.PCG()
# i<-0
# lapply(split(genes, seq(nrow(genes))), function(x){
#   print(i)
#   i<<-i+1
#   draw.barplot.survival(gene = x$GeneID, symbol = x$symbol)
# })->a


#################################################################################### Univariate Regression
survdiff.gene <- function(gene) {
  mutate(clinical, time=ifelse(is.na(new_tumor), last, new_tumor), status=ifelse(is.na(new_tumor),1,2)) -> clinical.fix
  mutate(clinical2, time=ifelse(is.na(death), last, death), status=ifelse(is.na(death),1,2)) -> clinical2.fix
  get.data.surv <- function(clini.fix) {
    left_join(get.fpkm.anno(), clini.fix, by=c('case'='case'))%>%
      data.frame(fpkm=fpkm[,gene], stringsAsFactors = F) %>%
      filter(!type%in%c('11A','11B')) %>%
      distinct(case, .keep_all = T) %>%
      mutate(x.mean=ifelse(fpkm>mean(fpkm), 'high', 'low')) %>%
      mutate(x.median=ifelse(fpkm>median(fpkm), 'high', 'low'))
  }
  data.surv <- get.data.surv(clinical.fix)
  data.surv2 <- get.data.surv(clinical2.fix)
  survdiff.p <- function(diff) {
    pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
  }
  fds.p.mean <- survdiff.p(survdiff(Surv(time, status) ~ x.mean, data=data.surv))
  fds.p.median <- survdiff.p(survdiff(Surv(time, status) ~ x.median, data=data.surv))
  os.p.mean <- survdiff.p(survdiff(Surv(time, status) ~ x.mean, data=data.surv2))
  os.p.median <- survdiff.p(survdiff(Surv(time, status) ~ x.median, data=data.surv2))
  
  fd.cox.coefs <- summary(coxph(Surv(time, status) ~ fpkm, data=data.surv))$coefficients
  o.cox.coefs <- summary(coxph(Surv(time, status) ~ fpkm, data=data.surv2))$coefficients
  
  c(fds.p.mean, fds.p.median, os.p.mean, os.p.median, fd.cox.coefs[2], fd.cox.coefs[5], o.cox.coefs[2], o.cox.coefs[5])
}

genes <- helper.get.lncRNA.PCG()
i<-0
lapply(split(genes, seq(nrow(genes))), function(x){
  print(i)
  i<<-i+1
  survdiff.gene(gene = x$GeneID)
})->gene.surv.p
surv.p.data <- as.data.frame(t(matrix(unlist(gene.surv.p), byrow = F, ncol= length(gene.surv.p))))
colnames(surv.p.data) <- c('fds.mean','fds.median','os.mean', 'os.median', 'fd.cox.PH', 'fd.cox.p.value', 'fd.cox.PH', 'o.cox.p.value')
surv.p.data <- data.frame(gene=genes$GeneID, symbol=genes$symbol, surv.p.data,
                          fd.cox.p.adjust=p.adjust(surv.p.data$fd.cox.p.value, method='BH'),
                          o.cox.p.adjust=p.adjust(surv.p.data$o.cox.p.value, method='BH'),
                          stringsAsFactors = F)
surv.p.data <- arrange(surv.p.data, fds.mean)
write.csv(surv.p.data, file = './reports/surv.p.data.csv')

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


