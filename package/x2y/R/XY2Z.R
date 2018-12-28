checkXY2Z <- function(object) {
  errors <- character()
  if (!all(unique(object@detail$a)%in%colnames(object@x))) {
    msg <- "object@detail$a should in XY2Z@x"
    errors <- c(errors, msg)
  }
  if (!all(unique(object@detail$b)%in%colnames(object@y))) {
    msg <- "object@detail$b should in XY2Z@y"
    errors <- c(errors, msg)
  }
  if (nrow(object@x)!=nrow(object@y) | !all(rownames(object@x)==rownames(object@y))) {
    msg <- "rownames of object@x and rownames of object@y is not same"
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) TRUE else errors
}

setClass('XY2Z', slots = list(raw='data.frame',
                              FDR='numeric',
                              detail='data.frame',
                              x='matrix',
                              y='matrix',
                              nodes='list'),
         prototype = list(nodes=list(),FDR=0.05),
         validity=checkXY2Z)
##-----------------------------------
setGeneric('getNodes', def=function(object) {
  standardGeneric("getNodes")
})
setMethod("getNodes", 'XY2Z', function(object) {
  x.node<-as.vector(unique(object@detail$a))
  y.node<-as.vector(unique(object@detail$b))
  z.node<-unique(Reduce(union, lapply(split(object@detail,seq(nrow(object@detail))), function(xy){
    contingency <- object@x[,xy$a]*10+object@y[,xy$b]
    rownames(object@x)[contingency==11]
  })))
  object@nodes <- list(x=x.node, y=y.node, z=z.node)
})
##-----------------------------------
setGeneric('getZByX', def=function(object) {
  standardGeneric("getZByX")
})
setMethod("getZByX", 'XY2Z', function(object) {
  lapply(split(object@detail,as.vector(object@detail$a)), function(xys){
    unique(Reduce(union, lapply(split(xys,seq(nrow(xys))), function(xy){
      contingency <- object@x[,xy$a]*10+object@y[,xy$b]
      rownames(object@x)[contingency==11]
    })))
  })
})
##-----------------------------------
setGeneric('getYZByX', def=function(object) {
    standardGeneric("getYZByX")
})
setMethod("getYZByX", 'XY2Z', function(object) {
    lapply(split(object@detail,as.vector(object@detail$a)), function(xys){
        z <- unique(Reduce(union, lapply(split(xys,seq(nrow(xys))), function(xy){
            contingency <- object@x[,xy$a]*10+object@y[,xy$b]
            rownames(object@x)[contingency==11]
        })))
        y <- sapply(split(xys, seq(nrow(xys))), function(xy){
            xy$b
        })
        list(z,y)
    })
})
##-----------------------------------
setGeneric('getZByXY', def=function(object) {
  standardGeneric("getZByXY")
})
setMethod("getZByXY", 'XY2Z', function(object) {
  z <- lapply(split(object@detail,seq(nrow(object@detail))), function(xy){
    contingency <- object@x[,xy$a]*10+object@y[,xy$b]
    rownames(object@x)[contingency==11]
  })
  names(z) <- sapply(split(object@detail, seq(nrow(object@detail))), function(xy){
      paste0(xy$a, '-', xy$b)
  })
  z
})
##-------------------------------------
setGeneric('filterByX', def=function(object, x) {
    standardGeneric("filterByX")
})
setMethod("filterByX", signature(object='XY2Z',x='character'), function(object, x) {
    object@detail <- filter(object@detail, a%in%x)
    object@nodes <- getNodes(object)
    object
})
##-------------------------------------
setMethod("initialize", "XY2Z", function(.Object,...){
  .Object <- callNextMethod()
  detail <- .Object@raw%>%filter(FDR<.Object@FDR, phi>0, c11+c10>(c11+c10+c1+c0)/10)
  detail$a <- as.vector(detail$a)
  detail$b <- as.vector(detail$b)
  .Object@detail <- detail
  .Object@nodes <- getNodes(.Object)
  .Object
})
##-------------------------------------
setMethod('show', 'XY2Z', function(object){
  cat('lncRNA', 'tf', 'pcg','lncRNA-tf', sep = '\t');cat('\n')
  cat(length(object@nodes$x), length(object@nodes$y), length(object@nodes$z),nrow(object@detail), sep = '\t');cat('\n')
})
