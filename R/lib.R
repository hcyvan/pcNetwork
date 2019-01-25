library(stringr)
library(dplyr)


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

##' Covert symbol to gene ensembl id
##'
##' Covert symbol to gene ensembl id
##' @title 
##' @param symbol 
##' @return gene ensembl id
##' @author Navy cheng
pf.symbol2emsembl = function(symbol) {
    biomart <- pf.get.biomart()
    biomart[match(symbol,biomart$hgnc_symbol),]$ensembl_gene_id
}


##' Get differential expression genes
##'
##' @title 
##' @param type
##' @return differential expression genes filtered by GeneType
##' @author c509
pc.get.diff <- function(type=c('all', 'pcg', 'lncRNA')) {
    type <- unique(match.arg(type, several.ok=TRUE))
    diff <- readRDS('./cache/diff.3193.rds')
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
