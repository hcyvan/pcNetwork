library(pcProfile)
library(affy)
# library(affyQCReport)
# library(simpleaffy)
# library(affyPLM)
library(annotate)
library(hgu95c.db)
library(hgu95b.db)
library(hgu95av2.db)


# idir <- './data/AR/AR-chip/GSE6606_Primary'
# idir <- './data/AR/AR-chip/GSE6604_Normal'
# idir <- './data/AR/AR-chip/GSE6605_Metastatic'
# cel.p <- list.celfiles(path=idir,full.names=TRUE);
# for (cel in cel.p) {
#   data.p = ReadAffy(filenames =cel)
#   cdf <- data.p@cdfName
#   print(cdf)
#   if (!dir.exists(file.path(idir, cdf))) {
#     dir.create(file.path(idir, cdf))
#   }
#   file.rename(cel, file.path(idir, cdf, basename(cel)))
# }

# data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6604_Normal/HG_U95Av2')
# sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
# saqc = qc(data)
# plot(saqc)
# dataPLM = fitPLM(data)
# boxplot(dataPLM, main="NUSE", outline = FALSE, col="lightblue", las=3, whisklty=0, staplelty=0)
# Mbox(dataPLM, main="RLE", outline = FALSE, col="mistyrose", las=3, whisklty=0, staplelty=0)
# datarma = rma(data)
# dataexp <- exprs(datarma)
# rownames(dataexp)<-getSYMBOL(rownames(dataexp), 'hgu95av2')

########################### Normal + Primary + Metastatic
cel1 <- list.celfiles(path='./data/AR/AR-chip/GSE6604_Normal/HG_U95Av2',full.names=TRUE)
cel2 <- list.celfiles(path='./data/AR/AR-chip/GSE6606_Primary/HG_U95Av2',full.names=TRUE)
cel3 <- list.celfiles(path='./data/AR/AR-chip/GSE6605_Metastatic/HG_U95Av2',full.names=TRUE)
cel <- c(cel1,cel2,cel3)
data <- ReadAffy(filenames = cel)
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpNPM1 <- exprs(datarma)
rownames(dataexpNPM1)<-getSYMBOL(rownames(dataexpNPM1), 'hgu95av2')
dataexpNPM1<-dataexpNPM1[!is.na(rownames(dataexpNPM1)),]
saveRDS(dataexpNPM1, './data/dataexpNPM1.rds')
##########################################################
library(limma)
dataexpNPM1 <- readRDS('./data/dataexpNPM1.rds')
group <- c(rep('N',18),rep('P',65),rep('M',25))
design <- model.matrix(~factor(group))


fit <- lmFit(dataexpNPM1,design)
fit2 <- contrasts.fit(fit,coef=2)
fit2 <- eBayes(fit2)
diff2<-topTable(fit2,n=10000)%>%filter(adj.P.Val<0.05,abs(logFC)>=log(1.5,2))

fit3 <- contrasts.fit(fit,coef=3)
fit3 <- eBayes(fit3)
diff3<-topTable(fit3,n=10000)%>%filter(adj.P.Val<0.05,abs(logFC)>=log(1.5,2))


##############################################################
###########################################
########################## Normal
data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6604_Normal/HG_U95Av2')
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpN1 <- exprs(datarma)
rownames(dataexpN1)<-getSYMBOL(rownames(dataexpN1), 'hgu95av2')

data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6604_Normal/HG_U95B')
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpN2 <- exprs(datarma)
rownames(dataexpN2)<-getSYMBOL(rownames(dataexpN2), 'hgu95b')

data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6604_Normal/HG_U95C')
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpN3 <- exprs(datarma)
rownames(dataexpN3)<-getSYMBOL(rownames(dataexpN3), 'hgu95c')

########################## Primary
data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6606_Primary/HG_U95Av2')
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpP1 <- exprs(datarma)
rownames(dataexpP1)<-getSYMBOL(rownames(dataexpP1), 'hgu95av2')

data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6606_Primary/HG_U95B')
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpP2 <- exprs(datarma)
rownames(dataexpP2)<-getSYMBOL(rownames(dataexpP2), 'hgu95b')

data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6606_Primary/HG_U95C')
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpP3 <- exprs(datarma)
rownames(dataexpP3)<-getSYMBOL(rownames(dataexpP3), 'hgu95c')

########################### Metastatic
data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6605_Metastatic/HG_U95Av2')
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpM1 <- exprs(datarma)
rownames(dataexpM1)<-getSYMBOL(rownames(dataexpM1), 'hgu95av2')

data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6605_Metastatic/HG_U95B')
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpM2 <- exprs(datarma)
rownames(dataexpM2)<-getSYMBOL(rownames(dataexpM2), 'hgu95b')

data<-ReadAffy(celfile.path = './data/AR/AR-chip/GSE6605_Metastatic/HG_U95C')
sampleNames(data) = sub("\\.CEL.gz$", "", sampleNames(data))
datarma = rma(data)
dataexpM3 <- exprs(datarma)
rownames(dataexpM3)<-getSYMBOL(rownames(dataexpM3), 'hgu95c')

saveRDS(dataexpN1, './data/dataexpN1.rds')
saveRDS(dataexpN2, './data/dataexpN2.rds')
saveRDS(dataexpN3, './data/dataexpN3.rds')
saveRDS(dataexpP1, './data/dataexpP1.rds')
saveRDS(dataexpP2, './data/dataexpP2.rds')
saveRDS(dataexpP3, './data/dataexpP3.rds')
saveRDS(dataexpM1, './data/dataexpM1.rds')
saveRDS(dataexpM2, './data/dataexpM2.rds')
saveRDS(dataexpM3, './data/dataexpM3.rds')
##########################################################
real <- readRDS('./data/trs2.rds')
dataexpP1<-readRDS('./data/dataexpP1.rds')
dataexpM1<-readRDS('./data/dataexpM1.rds')
dataexpP2<-readRDS('./data/dataexpP2.rds')
dataexpM2<-readRDS('./data/dataexpM2.rds')
dataexpP3<-readRDS('./data/dataexpP3.rds')
dataexpM3<-readRDS('./data/dataexpM3.rds')

tfs <- unique(c(as.vector(real$tf),'AR'))
targets <- unique(as.vector(real$target))
genes <- union(tfs, targets)


datageneP1<- dataexpP1[match(genes, rownames(dataexpP1)),]
datageneP1<-datageneP1[!is.na(rownames(datageneP1)),]
datageneM1<- dataexpM1[match(genes, rownames(dataexpM1)),]
datageneM1<-datageneM1[!is.na(rownames(datageneM1)),]
datageneP2<- dataexpP2[match(genes, rownames(dataexpP2)),]
datageneP2<-datageneP2[!is.na(rownames(datageneP2)),]
datageneM2<- dataexpM2[match(genes, rownames(dataexpM2)),]
datageneM2<-datageneM2[!is.na(rownames(datageneM2)),]
datageneP3<- dataexpP3[match(genes, rownames(dataexpP3)),]
datageneP3<-datageneP3[!is.na(rownames(datageneP3)),]
datageneM3<- dataexpM3[match(genes, rownames(dataexpM3)),]
datageneM3<-datageneM3[!is.na(rownames(datageneM3)),]

corP1 <- multicor(datageneP1,rds='./data/dataexpP1.cor.4.rds')
corM1 <- multicor(datageneM1,rds='./data/dataexpM1.cor.4.rds')
corP2 <- multicor(datageneP2,rds='./data/dataexpP2.cor.rds')
corM2 <- multicor(datageneM2,rds='./data/dataexpM2.cor.rds')
corP3 <- multicor(datageneP3,rds='./data/dataexpP3.cor.rds')
corM3 <- multicor(datageneM3,rds='./data/dataexpM3.cor.rds')

corP<-corP1
corM<-corM1

# filter(corP,FDR<0.05,r>0.7)%>%filter(v1=='IRF1'|v2=='IRF1')
# filter(corM,FDR<0.05,r>0.7)%>%filter(v1=='TP53'|v2=='TP53')
# filter(corM,FDR<0.05,r>0.7)%>%filter(v1=='IRF1'|v2=='IRF1')
filter(corM,FDR<0.05,abs(r)>0.7)%>%filter(v1=='AR'|v2=='AR')

finalP<-filter(corP1,FDR<0.05,abs(r)>0.7)
finalM<-filter(corM1,FDR<0.05,abs(r)>0.7)


a<-sort(sapply(tfs, function(x){
  nrow(filter(finalP,v1==x|v2==x))
}), decreasing=TRUE)
b<-sort(sapply(tfs, function(x){
  nrow(filter(finalM,v1==x|v2==x))
}),decreasing = TRUE)

node.tfs<-names(a)[1:8]
cytoP<-rbind(filter(finalP,v1%in%node.tfs),filter(finalP,v2%in%node.tfs))
write.csv(cytoP, './reports/thesis/cytoP.edge.csv',row.names = FALSE,quote = FALSE)
node.attr1<- data.frame(node=node.tfs,attr='tf')
node.attr2<- data.frame(node=setdiff(unique(union(cytoP$v1,cytoP$v2)),node.tfs),attr='gene')
node.attr<-rbind(node.attr1,node.attr2)
write.csv(node.attr, './reports/thesis/cytoP.node.csv',row.names = FALSE,quote = FALSE)


node.tfs<-names(b)[1:8]
cytoM<-rbind(filter(finalM,v1%in%node.tfs),filter(finalM,v2%in%node.tfs))
write.csv(cytoM, './reports/thesis/cytoM.edge.csv',row.names = FALSE,quote = FALSE)
node.attr1<- data.frame(node=node.tfs,attr='tf')
node.attr2<- data.frame(node=setdiff(unique(union(cytoM$v1,cytoM$v2)),node.tfs),attr='gene')
node.attr<-rbind(node.attr1,node.attr2)
write.csv(node.attr, './reports/thesis/cytoM.node.csv',row.names = FALSE,quote = FALSE)

sum(a!=0)
dim(filter(finalP,v1%in%tfs|v2%in%tfs))
dim(filter(finalP,v1%in%tfs|v2%in%tfs, r>0))
dim(filter(finalP,v1%in%tfs|v2%in%tfs, r<0))
sum(b!=0)
dim(filter(finalM,v1%in%tfs|v2%in%tfs))
dim(filter(finalM,v1%in%tfs|v2%in%tfs, r>0))
dim(filter(finalM,v1%in%tfs|v2%in%tfs, r<0))



sort(table(trs2$tf),decreasing = TRUE)
###############################################################
co<-as.vector(filter(finalM,v1=='IRF1')$v2)
bi <-as.vector(filter(trs2,tf=='IRF1')$target)
intersect(co,bi)


#############################################################3