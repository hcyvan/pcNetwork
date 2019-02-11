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
    biomart[match(ensembl,biomart$ensembl_gene_id),]$hgnc_symbol
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

##' Covert symbol to gene ensembl id
##'
##' Covert symbol to gene ensembl id
##' @title 
##' @param symbol 
##' @return gene ensembl id
##' @author Navy Cheng
pf.symbol2emsembl <- function(symbol) {
    biomart <- pf.get.biomart()
    biomart[match(symbol,biomart$hgnc_symbol),]$ensembl_gene_id
}


##' Get differential expression genes
##'
##' @title 
##' @param type
##' @return differential expression genes filtered by GeneType
##' @author Navy Cheng
pf.get.diff <- function(type=c('all', 'pcg', 'lncRNA')) {
    type <- unique(match.arg(type, several.ok=TRUE))
    # diff <- readRDS('./cache/diff.3193.rds')
    diff <- readRDS('./support/diff.T_N.qlf.1.005.3069.rds')
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
  readRDS('./data/prad.rna.logFpkm.rds')
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
  readRDS('./data/prad.rna.zfpkm.rds')
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
