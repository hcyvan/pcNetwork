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
    m.adj <-getX2yMatrixAdjust(pair1, pair2, background)
    if (refresh || is.null(rds) || !file.exists(rds)) {
        message('Calulate x2y ...\n')
        x2y <- xyCor(m.adj$a, m.adj$b, cores=6)
        if (!is.null(rds)) {
            saveRDS(x2y, file = rds)
        }
    } else {
        message(paste0('Load x2y from ', rds, '...\n'))
        x2y <- readRDS(rds)
    }
    new('XY2Z', raw=x2y, x=m.adj$a, y=m.adj$b)
}
