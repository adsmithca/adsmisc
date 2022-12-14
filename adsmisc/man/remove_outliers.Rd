% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{remove_outliers}
\alias{remove_outliers}
\title{Remove outliers}
\usage{
remove_outliers(x, na.rm = TRUE, ...)
}
\arguments{
\item{x}{a numeric vector}

\item{na.rm}{TRUE or FALSE}

\item{...}{other arguments to be passed}
}
\value{

}
\description{
Remove outliers
}
\examples{
x <- c(1,2,3,4,5,6,5,4,3,2,1,99); remove_outliers(x)
df <- data.frame(y = seq(1:10), x = c(51,52,53,54,55,55,54,53,52,99), g = c(rep("a",5),rep("b",5)))
library(tidyverse);library(magrittr)
df \%<>\% dplyr::group_by(g) \%>\% dplyr::summarise(x = remove_outliers(x))
}
