library(ProjectTemplate)

load.project()

rnaWithAnno <- left_join(diffLncRnaAndPcg, bioMart, by=c('GeneID'='Gene.stable.ID')) %>%
  select(colnames(diffRna), 
         tss = Transcription.start.site..TSS.,
         geneStart = Transcript.start..bp., 
         geneEnd = Transcript.end..bp., 
         chromosome = Chromosome.scaffold.name)



findGeneAround <- function(genes, up, down) {
  locas <- list()
  for (chromo in unique(genes[,'chromosome'])) {
    if (!is.na(chromo)) {
      print(paste('checking chromosome', chromo))
      chromoGenes <- arrange(filter(genes, chromosome==chromo), tss)
      len <- dim(chromoGenes)[1]
      upNum <- list()
      for (u in up) {
        upNum[[paste0('up', u)]] <- c()
        for(i in seq(len)) {
          tss <- chromoGenes[i, 'tss']
          j <- i
          jNum <- 0
          while (j > 1) {
            j <- j - 1
            if (tss - chromoGenes[j, 'tss'] < u) {
              jNum <- jNum + 1
            } else {
              break
            }
          }
          upNum[[paste0('up', u)]] <- c(upNum[[paste0('up', u)]], jNum)
        }
      }
      downNum <- list()
      for (d in down) {
        downNum[[paste0('down', d)]] <- c()
        for(i in seq(len)) {
          tss <- chromoGenes[i, 'tss']
          k <- i
          kNum <- 0
          while (k < len) {
            k <- k + 1
            if (chromoGenes[k, 'tss'] - tss < d) {
              kNum <- kNum + 1
            } else {
              break
            }
          }
          downNum[[paste0('down', d)]] <- c(downNum[[paste0('down', d)]], kNum)
        }
      }
      locas[[paste0('ch',chromo)]] <- data.frame(chromoGenes, as.data.frame(upNum), as.data.frame(downNum))
    }
  }
  locas
}

geneAround <- findGeneAround(rnaWithAnno,
                               up = c(1000, 5000, 10000, 20000,30000),
                               down = c(1000, 5000, 10000, 20000,30000))

geneAroundFrame = data.frame()
for (l in geneAround) {
  geneAroundFrame <- rbind(geneAroundFrame, l)
}

write.csv(geneAroundFrame, file = './reports/geneAround.csv')
