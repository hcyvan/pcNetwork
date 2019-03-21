library(stringr)
library(dplyr)
library(pcProfile)


pv.lncRNA <- c('lincRNA',
               'bidirectional_promoter_lncRNA',
               '3prime_overlapping_ncRNA',
               'macro_lncRNA',
               'antisense',
               'sense_overlapping',
               'sense_intronic')
pv.pcg <- c('protein_coding',
            'IG_C_gene',
            'IG_J_gene',
            'TR_D_gene',
            'TR_V_gene',
            'TR_J_gene',
            'IG_D_gene',
            'IG_V_gene',
            'TR_C_gene')



##' Get and cache data
##'
##' Read data from cache if it read before. Otherwise, get data from
##' func and cache it.
##' @title 
##' @return 
##' @author Navy Cheng
##' @param key key to cache data
##' @param func how to get data
##' @param ... 
.pf.cache <- function(key, func, ...){
    env = globalenv()
    if (exists(key, envir = env)) {
        get(key, envir = env)
    } else {
        data <- func(...)
        assign(key, data, envir = env)
        data
    }
}

##' get biomart data
##'
##' get biomart data. Which is 
##' @title 
##' @return biomart
##' @author Navy Cheng
pf.get.biomart <- function() {
    .pf.cache('.biomart', func = function(){
        readRDS('./support/biomart.rds')
    })
}

##' Tranlate gene ensembl id to symbol
##'
##' Tranlate gene ensembl id to symbol
##' @title 
##' @return gene symbol
##' @author Navy Cheng
##' @param ensembl 
pf.ensembl2symbol <- function(ensembl) {
    biomart <- pf.get.biomart()
    biomart[match(ensembl,biomart$ensembl_gene_id),]$external_gene_name
}

##' Tranlate gene ensembl id to biotype
##'
##' Tranlate gene ensembl id to biotype
##' @title 
##' @return gene biotype
##' @author Navy Cheng
##' @param ensembl 
pf.ensembl2biotype <- function(ensembl) {
    biomart <- pf.get.biomart()
    biomart[match(ensembl,biomart$ensembl_gene_id),]$gene_biotype
}

pf.ensembl2chr <- function(ensembl) {
  biomart <- pf.get.biomart()
  biomart[match(ensembl,biomart$ensembl_gene_id),]$chromosome_name
}

##' Covert symbol to gene ensembl id
##'
##' Covert symbol to gene ensembl id
##' @title 
##' @param symbol 
##' @return gene ensembl id
##' @author Navy Cheng
pf.symbol2emsembl <- function(symbol) {
    biomart <- pf.get.biomart()
    biomart[match(symbol,biomart$external_gene_name),]$ensembl_gene_id
}


##' Get differential expression genes
##'
##' @title 
##' @param type
##' @return differential expression genes filtered by GeneType
##' @author Navy Cheng
pf.get.diff <- function(type=c('all', 'pcg', 'lncRNA'),n='5946') {
    type <- unique(match.arg(type, several.ok=TRUE))
    # diff <- readRDS('./cache/diff.3193.rds')
    if (n=='5946') {
      diff <- readRDS('./support/diff.T_N.qlf.log1.5.005.5946.rds')
    } else if (n=='3096') {
      diff <- readRDS('./support/diff.T_N.qlf.1.005.3069.rds')
    } else {
      stop()
    }
    if ('all'%in%type){
        diff
    } else {
        data <- data.frame()
        for(t in type) {
            if (t=='lncRNA') {
                data <- rbind(diff[diff$GeneType%in%pv.lncRNA,], data)
            } else if (t=='pcg') {
                data <- rbind(diff[diff$GeneType%in%pv.pcg,], data)
            }
        }
        data
    }
}

pf.get.diff.ar<-function(){
  readRDS('support/ar.diff.rds')
}

##' Get TCGA-PRAD count matrix
##'
##' @title Get count matrix
##' @return count matrix
##' @author Navy Cheng
pf.get.count <- function(refresh=FALSE) {
    load('./support/prad.rna.count.rda')
    prad.rna.count
}

##' Get log2FPKM matrix
##'
##' @title Get logFPKM matrix
##' @return logFpkm matrix
##' @author Navy Cheng
pf.get.logFpkm <- function() {
  .pf.cache('.logFpkm', func = function(){
    readRDS('./data/prad.rna.logFpkm.rds')
  })
}

##' Filter zFPKM matrix
##'
##' @title Filter logFPKM matrix
##' @param ensembl ensembl ids
##' @return logFPKM matrix
##' @author Navy Cheng
pf.filter.logFpkm <- function(ensembl, rm.na=FALSE) {
  data <- pf.get.logFpkm()
  index <- match(ensembl, rownames(data))
  if (rm.na) {
    index <- na.omit(index)
  }
  data[index,]
}

##' Get z-scored log2FPKM matrix
##'
##' @title Get zFPKM matrix
##' @return zfpkm matrix
##' @author Navy Cheng
pf.get.zfpkm <- function() {
  .pf.cache('.zfpkm', func = function(){
    readRDS('./data/prad.rna.zfpkm.rds')
  })
}

##' Filter zFPKM matrix
##'
##' @title Filter zFPKM matrix
##' @param ensembl ensembl ids
##' @return zFPKM matrix
##' @author Navy Cheng
pf.filter.zfpkm <- function(ensembl, rm.na=FALSE) {
  data <- pf.get.zfpkm()
  index <- match(ensembl, rownames(data))
  if (rm.na) {
    index <- na.omit(index)
  }
  data[index,]
}

##' Get FPKM matrix
##'
##' @title Get FPKM matrix
##' @param refresh reload FPKM
##' @return fpkm matrix
##' @author Navy Cheng
pf.get.fpkm <- function(refresh=FALSE) {
  rds.path <- './cache/data.fpkm.rds'
  csv.path <- './data/data.fpkm.csv'
  if (file.exists(rds.path)&& !refresh) {
    data.fpkm <- readRDS(rds.path)
  } else {
      data.fpkm <- read.delim(csv.path,
                              sep = ',',
                              header = TRUE,
                              stringsAsFactors = FALSE,
                              row.names = 'GeneID')
    rownames(data.fpkm) <- str_split_fixed(rownames(data.fpkm), '\\.',2)[,1]
    saveRDS(data.fpkm, file = rds.path)
  }
  data.fpkm
}

##' Filter FPKM matrix
##'
##' @title Filter FPKM matrix
##' @param ensembl ensembl ids
##' @param refresh reload FPKM
##' @return FPKM matrix
##' @author Navy Cheng
pf.filter.fpkm <- function(ensembl, rm.na=FALSE, refresh=FALSE) {
  data <- pf.get.fpkm(refresh=FALSE)
  index <- match(ensembl, rownames(data))
  if (rm.na) {
    index <- na.omit(index)
  }
  data[index,]
}

##' Filter gene annotation
##' 
##' @title Filter gene annotation
##' @param ensembl ensembl ids
##' @return annotations
##' @author c509
pf.filter.anno <- function(ensembl) {
  biomart <- pf.get.biomart()
  biomart[match(ensembl, biomart$ensembl_gene_id),]
}

#' pcNetwork
getGene2tfMatrix <- function(tf2gene, genes, tfs.filter=NULL){
  gene2tf.split <- split(tf2gene, tf2gene$gene)[genes]
  tfs <- Reduce(union, lapply(gene2tf.split, function(x){x$tf}))
  if (!is.null(tfs.filter)) {
    tfs <- tfs.filter
  }
  gene2tf.l <- lapply(gene2tf.split, function(x){
    tmp <- x$N[match(tfs,x$tf)]
    if(length(tmp)==0){return(NA)}
    tmp
  })
  gene2tf.m <-do.call(rbind, gene2tf.l)
  rownames(gene2tf.m) <- genes
  colnames(gene2tf.m) <- tfs
  gene2tf.m[is.na(gene2tf.m)] <- 0
  gene2tf.m
}

##' calculate multiple cor
##' 
##' @title Filter gene annotation
##' @param m1 matrix,sample is column
##' @param m2 matrix,sample is column
##' @return annotations
##' @author navych
pf.multicor <- function(m1, m2=NULL, method= c('pearson', 'kendall', 'spearman'), rds=NA, rewrite=FALSE, verbose=TRUE) {
  if(!is.na(rds)) {
    if(file.exists(rds) && !rewrite) {
      message(paste('Load data from', rds))
      return(readRDS(rds))
    }
  }
  method <- match.arg(method)
  if (is.null(m2)) {
    data1 <- data2 <- as.matrix(m1)
  } else {
    data1 <- as.matrix(m2)
    data2 <- as.matrix(m1)
  }
  r <- matrix(nrow=ncol(data2),ncol=ncol(data1))
  if (verbose) {
    total <- ncol(data1)
    message(paste('Calculating:', total,'rounds needed!'))
  }
  i <- 0
  p <- apply(data1,2,function(x){
    i<<-i+1
    t0 <- Sys.time()
    if (verbose) {
      cat(paste0(i,'/',total,'\n'))
    }
    j<-0
    apply(data2,2,function(y){
      j<<-j+1
      if (!is.null(data2)|i > j) {
        test <- cor.test(x,y,method = method)
        r[j,i] <<- test$estimate
        test$p.value
      } else {
        NaN
      }
    })
  })
  dimnames(r) <- dimnames(p)
  list(r,p)
  r.melt <- reshape2::melt(r) %>% filter(!is.na(value))
  p.melt <- reshape2::melt(p) %>% filter(!is.na(value))
  ret <- data.frame(v1=r.melt$Var1, v2=r.melt$Var2, r=r.melt$value, p.value=p.melt$value, FDR=p.adjust(p.melt$value, method = 'BH'))
  if (!is.na(rds)) {
    saveRDS(ret, rds)
  }
  ret
}

pf.get.ac<-function(){
  readRDS('data/tf.ac.rds')
}

pf.get.dug<-function(){
  readRDS('./data/coxph.dfs.nui.dug.rds')
}
pf.get.oug<-function(){
  readRDS('./data/coxph.os.nui.oug.rds')
}
pf.get.sugid<-function(){
  oug<-pf.get.oug()
  dug<-pf.get.dug()
  unique(union(oug$gene,dug$gene))
}
pf.get.lnc.cis<-function(up=1000,uniq=TRUE){
  lnc.gene.raw<-readRDS('support/lnc.cis.origin.rds')
  lnc.gene<-data.frame(lncRNA=as.vector(lnc.gene.raw$lncRNA),
                       gene=as.vector(lnc.gene.raw$gene),
                       dist=lnc.gene$dist,
                       stringsAsFactors = FALSE)
  lnc.gene<-filter(lnc.gene,dist<=up)
  lnc.cis<-filter(lnc.gene,lncRNA!=gene)
  if(uniq){
    lnc.cis<-distinct(lnc.cis,lncRNA,gene,.keep_all = TRUE)
  }
  lnc.cis
}

pf.get.lnc.trans<-function(){
  lnc.trans.raw<-read.delim('data/lncRNA.DNA.triplex.txt',
                            stringsAsFactors = FALSE)
  mutate(lnc.trans.raw,lncRNA=Sequence.ID, gene=Duplex.ID)%>%
    distinct(lncRNA,gene,.keep_all = TRUE)%>%
    select(lncRNA,gene)
}
pf.get.lnc.cis.cor<- function(){
  readRDS('support/lnc.cis.cor.rds')
}
pf.get.lnc.trans.cor<- function(){
  readRDS('support/lnc.trans.cor.rds')
}
pf.get.lnctf.cor.indirect<-function(){
  readRDS('support/lnctf.cor.indirect.rds')
}
pf.get.ltg.cis<-function(){
  readRDS('support/lnc.tf.gene.cis.rds')
}
pf.get.ltg.trans<-function(){
  readRDS('support/lnc.tf.gene.trans.rds')
}

pf.get.loc<-function(){
  readRDS('./cache/cell.location.rds')
}



