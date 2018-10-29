library(dplyr)

source('./lib/globals.R')
source('./lib/helpers.R')

data.gene.factor <- read.delim('./data/data.gene.factor.csv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
modules <- read.delim('./data/diff.qlf.2877.wgcna.color.csv', sep = ',',header = TRUE, stringsAsFactors = FALSE)

genes <- unique(data.gene.factor$GeneID)
tfs <- unique(data.gene.factor$Factor)

tf2genes.cache <- './cache/tf2gene.rda'
gene2tfs.cache <- './cache/gene2tfs.rda'

# ===================== pre-handle data.gene.factor ==================== 
tf2genes <- list()
i <- 1
for (tf in tfs) {
  cat(paste0(i, '\n'))
  i <- i+1
  tf2genes[[tf]] <- (data.gene.factor %>% filter(Factor==tf))$GeneID
}
save(tf2genes, file = tf2gene.cache)

gene2tfs <- list()
i <- 1
for (id in genes) {
  cat(paste0(i, '\n'))
  i <- i+1
  gene2tfs[[id]] <- (data.gene.factor %>% filter(GeneID==id))$Factor
}
save(gene2tfs, file = gene2tfs.cache)

load(tf2genes.cache)
load(gene2tfs.cache)

tf2genes.num <- sapply(tf2genes, length)
barplot(sort(tf2genes.num ))

gene2tfs.num <- sapply(gene2tfs, length)
barplot(sort(gene2tfs.num))
# ---------------------------------------------------------------------



data.gene.factor %>% filter(Factor=='AR') -> gene






summarise(group_by(modules, color),num=n())


black <- (modules %>% filter(color=='black'))$id


# -------------------------------- tfs -------------------------------------
# tfs <- gene.tf(modules$id, cache=FALSE)
tfs.cache <- './cache/diff.qlf.2877.tfs.rda'
# save(tfs, file=tfs.cache)
load(tfs.cache)

gene.tf <- function(ids, cache=TRUE) {
  tfs.tmp <- list()
  if(cache) {
    load(tfs.cache)
    for (id in ids) {
      tfs.tmp[[id]] <- tfs[[id]]
    }
  } else {
    i <- 1
    for (id in ids) {
      cat(paste0(i, '\n'))
      i <- i+1
      tfs.tmp[[id]] <- (data.gene.factor %>% filter(GeneID==id))$Factor
    }
  }
  tfs.tmp
}

tf.gene <- function(tfs) {
  genes.tmp <- list()
  i <- 1
  for (tf in tfs) {
    cat(paste0(i, '\n'))
    i <- i+1
    genes.tmp[[id]] <- (data.gene.factor %>% filter(Factor==tf))$GeneID
  }
  genes.tmp
}

tfs.intersect <- list() 
for (c in unique(modules$color)) {
  color.ids <- (modules %>% filter(color==c))$id
  tfs.intersect[[c]] <- Reduce(intersect, gene.tf(color.ids))
}

tfs.intersect <- list() 
for (c in unique(modules$color)) {
  color.ids <- (modules %>% filter(color==c))$id
  tfs.intersect[[c]] <- Reduce(intersect, gene.tf(color.ids))
}

tfs.num <- sapply(tfs, length)
barplot(sort(tfs.num))























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


