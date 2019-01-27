#' Adjust the row of tow xy.struct
#' 
#' Adjust the tow matix so that they share the same row names
#' 
#' @param xa2y Input data with a xy.matrix struct
#' @param xb2y Input data with a xy.matrix struct
#' @return \item{name fix.list}{description a list contain the fixed result}
#' @keywords filename 
#' @export
#' @examples
#' xa2y <- list(x1=c('y1','y2','y3'), x2=c('y2','y3'), x3=c('y1','y2'))
#' xa2b.m <- x2yMatirx(a2b)
#' xb2y <- list(x1=c('y1','y2','y3','y4'), x2=c('y2','y3'), x3=c('y1','y2','y5'))
#' xb2b.m <- x2yMatirx(a2b)
#' x2yMatrixAdjust(xa2b.m, xb2b.m)

x2yMatrixAdjust <- function(xa2y, xb2y, mode=c('intersect', 'union'),y.names=NULL) {
  mode <- match.arg(mode)
  
  xa2y.row <- rownames(xa2y)
  xb2y.row <- rownames(xb2y)
  if (is.null(y.names)) {
    if (mode=='intersect') {
      y.names<-intersect(xa2y.row, xb2y.row)
    } else {
      y.names<-union(xa2y.row, xb2y.row)
    }
  }
  xa2y.ad <- xa2y[match(y.names, xa2y.row),]
  xb2y.ad <- xb2y[match(y.names, xb2y.row),]
  rownames(xa2y.ad) <- rownames(xb2y.ad) <- y.names
  xa2y.ad[is.na(xa2y.ad)] <- 0
  xb2y.ad[is.na(xb2y.ad)] <- 0
  list(a=xa2y.ad, b=xb2y.ad)
}

##' Get xy.struct
##'
##' @title Get xy.struct
##' @param xa x2b paris or list
##' @param xb x2b paris or list
##' @param background a set of `b`
##' @export
##' @return \item{name fix.list}{description a list contain the fixed result}
##' @author c509
getX2yMatrixAdjust <- function(xa, xb, background=NULL) {
  a.m <- x2yMatrix(xa, b=background)
  b.m <- x2yMatrix(xb, b=background)
  fix <- x2yMatrixAdjust(a.m, b.m, y.names = pcg$GeneID)
  fix
}
