##' Enrich two sets
##'
##' @title double enrich
##' @param pair1 pair1
##' @param pair2 pair2
##' @param background background set
##' @param rds rds file path
##' @param cores cpu core numbers
##' @param refresh re caculate
##' @return XY2Z
##' @export
##' @author Navy Cheng
dEnricher <- function(pair1, pair2, background=NULL, rds=NULL, cores=6, refresh=FALSE) {
    if (refresh || is.null(rds) || !file.exists(rds)) {
        message('Calulate XY2Z ...\n')
        m.adj <-getX2yMatrixAdjust(pair1, pair2, background)
        x2y <- xyCor(m.adj$a, m.adj$b, cores=6)
        xy2z <- new('XY2Z', raw=x2y, x=m.adj$a, y=m.adj$b)
        if (!is.null(rds)) {
            saveRDS(xy2z, file = rds)
        }
        xy2z
    } else {
        message(paste0('Load XY2Z from ', rds, '...\n'))
        readRDS(rds)
    }
}
