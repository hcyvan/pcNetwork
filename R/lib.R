library(stringr)
library(dplyr)

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
.pc.cache = function(key, func, ...){
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
pc.get.biomart = function() {
    .pc.cache('.biomart', func = function(){
        readRDS('./support/biomart.rds')
    })
}

##' Tranlate gene ensembl id to symbol
##'
##' Tranlate gene ensembl id to symbol
##' @title 
##' @return gene symbol
##' @author Navy Cheng
pc.ensembl2symbol = function(ensembl) {
    biomart <- pc.get.biomart()
    biomart[match(ensembl,biomart$ensembl_gene_id),]$hgnc_symbol
}
##' Covert symbol to gene ensembl id
##'
##' Covert symbol to gene ensembl id
##' @title 
##' @param symbol
##' @return gene ensembl id
##' @author Navy cheng
pc.symbol2emsembl = function(symbol) {
    biomart <- pc.get.biomart()
    biomart[match(symbol,biomart$hgnc_symbol),]$ensembl_gene_id
}
