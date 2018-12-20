#' Change x2y struct to xy.matrix
#' 
#' @details
#' x2y struct is a list like this
#' 
#' \code{
#' List of 3:\cr
#'   $ x1: chr[1:3] 'y1' 'y2' 'y3'\cr
#'   $ x2: chr[1:2] 'y2' 'y3'\cr
#'   $ x3: chr[1:2] 'y1' 'y2'\cr
#'}
#'
#' xy.matrix struct is a matrix like this
#' 
#' \code{
#'     x1  x2  x3\cr
#' y1  1   0   1\cr
#' y2  1   1   1\cr
#' y3  1   1   0\cr
#'}
#' 
#' @param a2b Input data with a x2y struct
#' @return \item{name xy.matrix}{description a value transformed from x2y}
#' @keywords filename 
#' @export
#' @examples
#' a2b <- list(x1=c('y1','y2','y3'), x2=c('y2','y3'), x3=c('y1','y2'))
#' x2yMatirx(a2b)

x2yMatrix <- function(a2b) {
  a <- names(a2b)
  b <- Reduce(union,a2b)
  
  lapply(a2b, function(t){
    as.numeric(b%in%t)
  }) -> a2b.loc
  a2b.m <- matrix(unlist(a2b.loc), ncol = length(a2b.loc))
  dimnames(a2b.m) <- list(b,a)
  a2b.m
}
