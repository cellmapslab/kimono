% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer_sgl_model.R
\name{calc_alpha}
\alias{calc_alpha}
\title{Calculate alpha based on frobenius norm and group information}
\usage{
calc_alpha(x, group)
}
\value{
Frobenius norm of x
}
\description{
Calculate alpha based on frobenius norm and group information
}
\examples{
x <- matrix(1:20,5,4)
group <- c(1,1,2,2)
calc_alpha(x, group)
}
