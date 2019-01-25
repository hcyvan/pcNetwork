library(dplyr)
source('./lib/globals.R')
source('./lib/helpers.R')

all <- helper.get.lncRNA.PCG()
pcg <- filter(all, GeneType%in%config$PCGs)
############################################### GET x.2.y
get.lncRNA.2.PCG <- function(s=0, fdr=0.05) {
  env = globalenv()
  key = paste0('.lncRNA.2.PCG.',s)
  if (exists(key, envir = env)) {
    get(key,envir = env)
  } else {
    cor.pairs <- readRDS('./cache/cor.pairs.info.rds')
    data <- filter(cor.pairs, type1%in%config$lncRNA & type2%in%config$PCGs & FDR<fdr & abs(r)>=s)
    lncRNA.2.PCG <- lapply(split(data,as.vector(data$v1)), function(x){as.vector(x$v2)})
    assign(key, lncRNA.2.PCG, envir = env)
    lncRNA.2.PCG
  }
}

get.tf.2.PCG.from.fimo <- function() {
  env = globalenv()
  key = '.tf.2.PCG'
  if (exists(key, envir = env)) {
    get(key, envir = env)
  } else {
    fimo <- read.delim('./data/enricher/fimo.460.2230.tsv', stringsAsFactors = F)
    fimo <- filter(fimo, motif_alt_id!='')
    tf.2.PCG <- lapply(split(fimo,as.vector(fimo$motif_alt_id)), function(x){unique(as.vector(x$sequence_name))})
    assign(key, tf.2.PCG, envir = env)
    tf.2.PCG
  }
}

plot.lncRNA.2.PCG <- function() {
  par(mfrow=c(3,3))
  for (s in seq(0.1, 0.9, 0.1)) {
    sapply(get.lncRNA.2.PCG(s), function(x){length(x)})%>%sort->tmp
    tmp%>%barplot(main=s)
    summary(tmp)
  }
  par(mfrow=c(1,1))
}

plot.tf.2.PCG <- function() {
  sapply(get.tf.2.PCG.from.fimo(), function(x){length(x)})%>%sort->tmp
  tmp%>%barplot()
  summary(tmp)
}
# plot.lncRNA.2.PCG()
# plot.tf.2.PCG()
################################################# tf-pcg
# get.tf.2.PCG.from.enricher <- function() {
#   tf.enricher <- read.delim('./data/enricher.all.bg0.new.csv', sep = ',', stringsAsFactors = F)
#   tf.enricher.list <- setNames(split(tf.enricher, seq(nrow(tf.enricher))), tf.enricher$detail.ID)
#
#   load('./cache/biomart.symbol.biotype.rda')
#
#   lapply(tf.enricher.list, function(tf){
#     symbol <- str_split(tf$detail.geneID, '/')[[1]]
#     filter(biomart.symbol.biotype, hgnc_symbol%in%symbol)$ensembl_gene_id
#   })
# }

#############################################################
devtools::load_all('./package/x2y/')

getFix <- function(s=0, tf='fimo') {
  if (tf=='enricher') {
    tf.2.PCG <- get.tf.2.PCG.from.enricher()
  } else {
    tf.2.PCG <- get.tf.2.PCG.from.fimo()
  }
  lncRNA.2.PCG <- get.lncRNA.2.PCG(s)
  lncRNA.2.PCG.m <- x2yMatrix(lncRNA.2.PCG)
  tf.2.PCG.m <- x2yMatrix(tf.2.PCG)
  fix <- x2yMatrixAdjust(lncRNA.2.PCG.m, tf.2.PCG.m,y.names = pcg$GeneID)
  fix
}


lncTF.all <- function(s=0) {
  fix <- getFix(s)
  message(paste('Calculate Phi of',ncol(fix$a), 'lncRNA and', ncol(fix$b),'TF'))
  lncRNA.tf.list <- xyCor(fix$a, fix$b)
  write.csv(lncRNA.tf.list, file = paste0('./data/lncRNA.tf.list/lncRNA.tf.',tf,'.',s,'.csv'))
  lncRNA.tf.list
}

# system.time(lncRNA.tf.fimo.0.1 <- lncTF.all(0.1))
# system.time(lncRNA.tf.fimo.0.2 <- lncTF.all(0.2))
# system.time(lncRNA.tf.fimo.0.3 <- lncTF.all(0.3))
# system.time(lncRNA.tf.fimo.0.4 <- lncTF.all(0.4))
# system.time(lncRNA.tf.fimo.0.5 <- lncTF.all(0.5))
# system.time(lncRNA.tf.fimo.0.6 <- lncTF.all(0.6))
# system.time(lncRNA.tf.fimo.0.7 <- lncTF.all(0.7))
# system.time(lncRNA.tf.fimo.0.8 <- lncTF.all(0.8))
# system.time(lncRNA.tf.fimo.0.9 <- lncTF.all(0.9))

# save(lncRNA.tf.fimo.0.1,
#      lncRNA.tf.fimo.0.2,
#      lncRNA.tf.fimo.0.3,
#      lncRNA.tf.fimo.0.4,
#      lncRNA.tf.fimo.0.5,
#      lncRNA.tf.fimo.0.6,
#      lncRNA.tf.fimo.0.7,
#      lncRNA.tf.fimo.0.8,
#      lncRNA.tf.fimo.0.9,
#      file = './cache/lncRNA.tf.fimo.x.rda'
#      )

############################################################
load(file = './cache/lncRNA.tf.fimo.x.rda')

fix1 <- getFix(0.1)
fix2 <- getFix(0.2)
fix3 <- getFix(0.3)
fix4 <- getFix(0.4)
fix5 <- getFix(0.5)
fix6 <- getFix(0.6)
fix7 <- getFix(0.7)
fix8 <- getFix(0.8)
fix9 <- getFix(0.9)


lncTP.0.1 <- new('XY2Z', raw=lncRNA.tf.fimo.0.1, x=fix1$a, y=fix1$b)
lncTP.0.2 <- new('XY2Z', raw=lncRNA.tf.fimo.0.2, x=fix2$a, y=fix2$b)
lncTP.0.3 <- new('XY2Z', raw=lncRNA.tf.fimo.0.3, x=fix3$a, y=fix3$b)
lncTP.0.4 <- new('XY2Z', raw=lncRNA.tf.fimo.0.4, x=fix4$a, y=fix4$b)
lncTP.0.5 <- new('XY2Z', raw=lncRNA.tf.fimo.0.5, x=fix5$a, y=fix5$b)
lncTP.0.6 <- new('XY2Z', raw=lncRNA.tf.fimo.0.6, x=fix6$a, y=fix6$b)
lncTP.0.7 <- new('XY2Z', raw=lncRNA.tf.fimo.0.7, x=fix7$a, y=fix7$b)
lncTP.0.8 <- new('XY2Z', raw=lncRNA.tf.fimo.0.8, x=fix8$a, y=fix8$b)
lncTP.0.9 <- new('XY2Z', raw=lncRNA.tf.fimo.0.9, x=fix9$a, y=fix9$b)


save(lncTP.0.1, lncTP.0.2, lncTP.0.3, lncTP.0.4, lncTP.0.5, lncTP.0.6, lncTP.0.7, lncTP.0.8, lncTP.0.9, file = './cache/lncTP.0.x.rda')
load('./cache/lncTP.0.x.rda')


lncTP.0.1
lncTP.0.2
lncTP.0.3
lncTP.0.4
lncTP.0.5
lncTP.0.6
lncTP.0.7
lncTP.0.8
lncTP.0.9

par(mfrow=c(3,3))
for(s in seq(0.1,0.9,0.1)) {
  lnctp <- get(paste0('lncTP.',s))
  cat('---------------------------',s,'----------------------------------------');cat('\n')
  print(lnctp)
  print(sort(lnctp@nodes$y))
  print(ncol(lnctp@detail))
  if (nrow(lnctp@detail) > 0) {
    barplot(sort(lnctp@detail$c11), main=s)
  }
}


write.csv(lncTP.0.3@nodes$y, './data/lnctp.tf.csv')
write.csv(lncTP.0.3@nodes$x, './data/lnctp.lncRNA.csv')
write.csv(lncTP.0.3@detail, './data/lnctp.detail.csv')
