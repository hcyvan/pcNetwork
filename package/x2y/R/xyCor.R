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

library(parallel)

xyCor <- function(am, bm, cores=NULL) {
    if (is.null(cores)){
        cores <- detectCores()
    }
    a.list <- mclapply(as.list(as.data.frame(am)), function(x){
        apply(bm,2,function(y){
            phi.test <- cor.test(x,y)
            con <- phiContingency(x,y)
            paste(c(phi.test$estimate, phi.test$p.value, con), collapse='/')
        })    
    }, mc.cores = cores)
    a2b.m <- do.call(rbind, a.list)
    rownames(a2b.m)<-colnames(am)
    a2b <- reshape2::melt(a2b.m)
    m <- str_split(as.vector(a2b$value),'\\/',simplify=TRUE)
    ret <- data.frame(a=a2b$Var1,
                      b=a2b$Var2,
                      phi=as.numeric(m[,1]),
                      p=as.numeric(m[,2]),
                      FDR=p.adjust(as.numeric(m[,2]), method='BH'),
                      c11=as.integer(m[,6]),
                      c10=as.integer(m[,5]),
                      c1=as.integer(m[,4]),
                      c0=as.integer(m[,3]),
                      stringsAsFactors = FALSE)
    ret[order(ret$FDR),]
}

phiContingency <- function(a,b) {
  # a and b looks like c(1,0,1,0,0,1)
  t <- as.data.frame(table(a*10+b))
  ret <- t[match(c(0,1,10,11), t$Var1),2]
  ret[is.na(ret)] <- 0
  ret
}
