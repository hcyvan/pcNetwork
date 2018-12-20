#' Calculate the relationship of tow x2y struct
#' 
#' @details
#' The tow x2y struct contain the same 'y' set {Y}. The contingency table represent the relationship of 'x' pairs with 'y'
#' 1 represent xa or xb have associate with 'y', and 0 represent no association. So, 11 represent both xa and xb have
#' associate with 'y' and 0 represent neither of them have association.
#' 
#' Phi and p-value of xa and xb is also caculated.
#' 
#' \code{
#'           b\cr
#'           1   0\cr           
#'    a   1  11  10\cr
#'        0  1   0\cr
#'}
#' 
#' @param am xa, input data with a fixed x2y struct
#' @param bm xb, input data with a fixed x2y struct
#' @return \item{name xy.relationship}{Phi, pvalue and contingency table of xa and xb}
#' @export
#' @examples
#' xa2y <- list(x1=c('y1','y2','y3'), x2=c('y2','y3'), x3=c('y1','y2'))
#' xa2b.m <- x2yMatirx(a2b)
#' xb2y <- list(x1=c('y1','y2','y3','y4'), x2=c('y2','y3'), x3=c('y1','y2','y5'))
#' xb2b.m <- x2yMatirx(a2b)
#' fix <- x2yMatrixAdjust(xa2b.m, xb2b.m)
#' xyCor(fix$a, fix$b)

xyCor <- function(am, bm,verbose=F) {
  an <- ncol(am)
  bn <- ncol(bm)
  matrix.phi <- matrix(nrow = an, ncol = bn)
  matrix.p <- matrix(nrow = an, ncol = bn)
  matrix.c.0 <- matrix(nrow = an, ncol = bn)
  matrix.c.1 <- matrix(nrow = an, ncol = bn)
  matrix.c.10 <- matrix(nrow = an, ncol = bn)
  matrix.c.11 <- matrix(nrow = an, ncol = bn)
  
  lapply(seq(ncol(am)), function(i){
    if (verbose) {
      message(paste(i,'/',an))
    }
    lapply(seq(ncol(bm)), function(j){
      phi.test <- cor.test(am[,i],bm[,j])
      con <- phiContingency(am[,i],bm[,j])
      matrix.phi[i,j] <<- phi.test$estimate
      matrix.p[i,j] <<- phi.test$p.value
      matrix.c.0[i,j] <<- con[1]
      matrix.c.1[i,j] <<- con[2]
      matrix.c.10[i,j] <<- con[3]
      matrix.c.11[i,j] <<- con[4]
    })
  })->tmp
  
  dimnames(matrix.phi) <- dimnames(matrix.p) <-
    dimnames(matrix.c.0) <- dimnames(matrix.c.1) <-
    dimnames(matrix.c.10) <- dimnames(matrix.c.11) <- list(colnames(am), colnames(bm))
  
  melt.phi <- reshape2::melt(matrix.phi)
  melt.p <- reshape2::melt(matrix.p)
  melt.c.0 <- reshape2::melt(matrix.c.0)
  melt.c.1 <- reshape2::melt(matrix.c.1)
  melt.c.10 <- reshape2::melt(matrix.c.10)
  melt.c.11 <- reshape2::melt(matrix.c.11)
  fdr <- p.adjust(melt.p$value, method = 'BH')
  
  ret <- data.frame(melt.phi, melt.p$value, fdr, melt.c.11$value, melt.c.10$value, melt.c.1$value, melt.c.0$value)
  colnames(ret) <- c('a', 'b', 'phi', 'p', 'FDR', 'c11', 'c10', 'c1', 'c0')
  ret[order(ret$FDR),]
}

phiContingency <- function(a,b) {
  # a and b looks like c(1,0,1,0,0,1)
  t <- as.data.frame(table(a*10+b))
  ret <- t[match(c(0,1,10,11), t$Var1),2]
  ret[is.na(ret)] <- 0
  ret
}
