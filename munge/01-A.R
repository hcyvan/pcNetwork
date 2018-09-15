diffRna <- read.delim(config$diffRnaFile, sep = ',')
diffRna <- diffRna[,-1]
diffLncRna <- filter(diffRna, GeneType %in% config$lncRNA)
diffPcg <- filter(diffRna, GeneType %in% config$PCGs)
diffLncRnaAndPcg <- filter(diffRna, GeneType %in% c(config$lncRNA, config$PCGs))

bioMart <- distinct(mart.export, Gene.stable.ID, .keep_all = TRUE)