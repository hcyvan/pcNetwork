source('./R/lib.R')
library(dplyr)


lnc.cis<-pf.get.lnc.cis(up=1000)
lnc.cis<-mutate(lnc.cis,lnc.sym=pf.ensembl2symbol(lncRNA),gene.sym=pf.ensembl2symbol(gene))

lyl<-read.csv('data/lyl/lncRNA.triplex.cis.and.trans.csv')

#####################3
gtrd <- fread('./data/Homo_sapiens_meta_clusters.interval', check.names = TRUE)
gtrd <- gtrd[,.(chrom=str_sub(X.CHROM, 4), start=START, end=END, len=END-START, tf=tfTitle)]

biomart.set <- as.data.table(pf.get.biomart())
biomart.set <- biomart.set[chromosome_name%in%gtrd$chrom,][,.(gene=ensembl_gene_id, transcript=ensembl_transcript_id, tss=transcription_start_site, up=transcription_start_site-1000,chrom=chromosome_name)]

# filter(biomart.set, gene==pf.symbol2emsembl('LMNTD2'))

biomart.map <- biomart.set[gene%in%pf.symbol2emsembl(lyl$gene)]
i<-0
lapply(split(biomart.map,seq(nrow(biomart.map))),function(x){
  i<<-i+1
  print(paste(i/nrow(gtrd)))
})

i <- 0
total <- nrow(gtrd)
genes<-mapply(function(chrom, start, end){
  i <<- i + 1
  if (i%%1000==0){
    print(paste(i%/%1000, total%/%1000))
  }
  tfbs <- biomart.map[[chrom]][up<=end&tss>=start,]
  if (nrow(tfbs)==0){
    gene <- NA
  }else{
    gene <- tfbs$gene[1]
  }
  gene
}, gtrd$chrom, gtrd$start, gtrd$end)


scanTss <- function(gtrd,v3=FALSE) {
  biomart.set <- as.data.table(pf.get.biomart())
  biomart.set <- biomart.set[chromosome_name%in%gtrd$chrom,][,.(gene=ensembl_gene_id, transcript=ensembl_transcript_id, tss=transcription_start_site, up=transcription_start_site-1000,chrom=chromosome_name)]
  biomart.map<-split(biomart.set, by = 'chrom')
  i <- 0
  total <- nrow(gtrd)
  genes<-mapply(function(chrom, start, end){
    i <<- i + 1
    if (i%%1000==0){
      print(paste(i%/%1000, total%/%1000))
    }
    tfbs <- biomart.map[[chrom]][up<=end&tss>=start,]
    if (nrow(tfbs)==0){
      gene <- NA
    }else{
      gene <- tfbs$gene[1]
    }
    gene
  }, gtrd$chrom, gtrd$start, gtrd$end)
  gtrd[,("gene"):=genes]
  if (v3) {
    gtrd[!is.na(gene), .(tf, gene)][,.N, by=c('tf','gene')]
  } else {
    unique(gtrd[!is.na(gene), .(tf, gene)])
  }
}