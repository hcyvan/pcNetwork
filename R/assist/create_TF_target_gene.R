#‘ Create TF-Gene pairs
source('./R/lib.R')
library(data.table)
library(stringr)

tf2gene <- function(pairs) {
  lapply(split(pairs,as.vector(pairs$tf)), function(x){unique(as.vector(x$gene))})
}

### Fimo TSS
fimo.tss.1 <- fread(paste0('./data/fimo.460.tss/fimo.460.tss.1000.1.tsv'), check.names = TRUE)
fimo.tss.2 <- fread(paste0('./data/fimo.460.tss/fimo.460.tss.1000.2.tsv'), check.names = TRUE)
fimo.tss.3 <- fread(paste0('./data/fimo.460.tss/fimo.460.tss.1000.3.tsv'), check.names = TRUE)
fimo.tss.4 <- fread(paste0('./data/fimo.460.tss/fimo.460.tss.1000.4.tsv'), check.names = TRUE)
fimo.tss <- rbind(fimo.tss.1, fimo.tss.2, fimo.tss.3, fimo.tss.4)
name <- str_split(fimo.tss$sequence_name, '\\|', simplify = TRUE)
fimo.tss[,c('gene','transcript','chromosome'):=list(name[,1], name[2], name[,3])]
# V1 TSS + 1e-5
fimo.tss.set <- unique(fimo.tss[p.value<=1e-5,.(tf=motif_alt_id, gene=gene)])
saveRDS(as.data.frame(fimo.tss.set), './cache/fimo.tss.set.rds')
# V2 TSS + 1e-5 + N
fimo.tss.filter <- fimo.tss[p.value<=1e-5]
fimo.tss.filter.unique <- unique(fimo.tss.filter, by=c('motif_alt_id','gene', 'start', 'strand'))
fimo.tss.set.2 <- fimo.tss.filter.unique[,.N, by=c('motif_alt_id','gene')][,.(tf=motif_alt_id,gene=gene,N=N)]
saveRDS(as.data.frame(fimo.tss.set.2), './cache/fimo.tss.set.2.rds')
# V3 TSS + 1e-4 + N
fimo.tss.filter <- fimo.tss[p.value<=1e-4]
fimo.tss.filter.unique <- unique(fimo.tss.filter, by=c('motif_alt_id','gene', 'start', 'strand'))
fimo.tss.set.3 <- fimo.tss.filter.unique[,.N, by=c('motif_alt_id','gene')][,.(tf=motif_alt_id,gene=gene,N=N)]
saveRDS(as.data.frame(fimo.tss.set.3), './cache/fimo.tss.set.3.rds')
### Fimo GSS
fimo.gss <- fread('./data/fimo.460.gss.tsv', stringsAsFactors = F, check.names = TRUE)
fimo.gss.set <- unique(fimo.gss[p.value<1e-4, .(tf=motif_alt_id, gene=sequence_name)])
saveRDS(as.data.frame(fimo.gss.set), './cache/fimo.gss.set.rds')
### TRUST
trrust<-fread('./data/trrust_rawdata.human.tsv')
trrust.mid <- trrust[,.(tf=V1,gene=pf.symbol2emsembl(V2))]
trrust.set <- unique(trrust.mid)[!is.na(gene)]
saveRDS(as.data.frame(trrust.set), './cache/trrust.set.rds')
# V2
trrust.set.2 <- trrust.mid[,.N, by=c('tf','gene')][!is.na(gene)]
saveRDS(as.data.frame(trrust.set.2), './cache/trrust.set.2.rds')
### Gtrd
gtrd <- fread('./data/Homo_sapiens_meta_clusters.interval', check.names = TRUE)
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

###---------------------------- GTRD
#### Gtrd: v1
gtrd.mid <- gtrd[,.(chrom=str_sub(X.CHROM, 4), start=START, end=END, len=END-START, tf=tfTitle)]
#gtrd.set <- scanTss(gtrd.mid)
#saveRDS(as.data.frame(gtrd.set), './cache/gtrd.set.rds')
gtrd.set <- as.data.table(readRDS('./cache/gtrd.set.rds'))
#### Gtrd: v2         peak.caller.count > 2
gtrd.mid.2 <- gtrd[peak.caller.count>2][,.(chrom=str_sub(X.CHROM, 4), start=START, end=END, len=END-START, tf=tfTitle)]
#' gtrd.set.peak.caller.count.2 <- scanTss(gtrd.mid.2)
#' saveRDS(as.data.frame(gtrd.set.peak.caller.count.2), './cache/gtrd.set.peak.call.count.2.rds')
gtrd.set.2 <- as.data.table(readRDS('./cache/gtrd.set.peak.call.count.2.rds'))
#### Gtrd: v3
gtrd.mid.3 <- gtrd[,.(chrom=str_sub(X.CHROM, 4), start=START, end=END, len=END-START, tf=tfTitle)]
gtrd.set.3 <- scanTss(gtrd.mid.3, v3=TRUE)
saveRDS(as.data.frame(gtrd.set.3), './cache/gtrd.set.3.rds')


#### Gtrd Prostate Cancer
cell.type <- read.csv('./data/cell.set.csv', stringsAsFactors = FALSE)
pc <- filter(cell.type,X==1)$x
gtrd.mid.pc <- gtrd[cell.set%in%pc,.(chrom=str_sub(X.CHROM, 4), start=START, end=END, len=END-START, tf=tfTitle)]
gtrd.set.pc <- scanTss(gtrd.mid.pc)
saveRDS(as.data.frame(gtrd.set.pc), './cache/gtrd.set.pc.rds')
gtrd.set.pc <- as.data.table(readRDS('./cache/gtrd.set.pc.rds'))
