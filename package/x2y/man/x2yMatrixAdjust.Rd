% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/x2yMatrixAdjust.R
\name{x2yMatrixAdjust}
\alias{x2yMatrixAdjust}
\title{Adjust the row of tow xy.struct}
\usage{
x2yMatrixAdjust(xa2y, xb2y, mode = c("intersect", "union"),
  y.names = NULL)
}
\arguments{
\item{xa2y}{Input data with a xy.matrix struct}

\item{xb2y}{Input data with a xy.matrix struct}
}
\value{
\item{name fix.list}{description a list contain the fixed result}
}
\description{
Adjust the tow matix so that they share the same row names
}
\examples{
xa2y <- list(x1=c('y1','y2','y3'), x2=c('y2','y3'), x3=c('y1','y2'))
xa2b.m <- x2yMatirx(a2b)
xb2y <- list(x1=c('y1','y2','y3','y4'), x2=c('y2','y3'), x3=c('y1','y2','y5'))
xb2b.m <- x2yMatirx(a2b)
x2yMatrixAdjust(xa2b.m, xb2b.m)
}
\keyword{filename}
