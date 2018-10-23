library(dplyr)

source('./lib/globals.R')
source('./lib/helpers.R')

data.gene.factor <- read.delim('./data/data.gene.factor.csv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

load('cache/diff.qlf.2877.pairs.rda')


lncNRAvsPCGs <- diff.cor.pairs %>%
                filter((type1%in%config$lncRNA & type2%in%config$PCGs) | (type1%in%config$PCGs & type2%in%config$lncRNA)) %>%
                filter(significant) %>%
                mutate_if(is.factor, as.character)

write.csv(lncNRAvsPCGs[,1:2], 'data/diff.qlf.2877.pairs.ids.csv', row.names = FALSE)

lncNRAvsPCGs %>% group_by(Var1) %>% summarise(count=n()) %>% select(gene=Var1, count=count) -> lncRNAs
lncNRAvsPCGs %>% group_by(Var2) %>% summarise(count=n()) %>% select(gene=Var2, count=count) -> PCGs


intersect2gene <- function(a, b) {
  (data.gene.factor %>% filter(GeneID==a))$Factor -> set.a
  (data.gene.factor %>% filter(GeneID==b))$Factor -> set.b
  list(a=set.a, b=set.b, intersect=intersect(set.a, set.b))
}

intersect2gene.inter.num <- function(a, b) {
  length(intersect2gene(a, b)$intersect)
}


m <- matrix(nrow = length(PCGs$gene), ncol = length(lncRNAs$gene))
colnames(m) <- lncRNAs$gene
rownames(m) <- PCGs$gene

system.time(
  
  for (i in c(1)) {
    for (j in 1:length(PCGs$gene)) {
      print(paste(i,j))
      m[i,j] = intersect2gene.inter.num(lncRNAs$gene[i], PCGs$gene[j])
    }
  }
  
)



outer(lncRNAs$Var1, PCGs$Var2, function(a, b){
  print(length(a))
  print(length(b))
  print(cbind(a,b)[1:10,])
})

a <- 'ENSG00000122548'
b <- 'ENSG00000000005'
c <- intersect2gene(a, b)
intersect2gene.inter.num(a,b)


