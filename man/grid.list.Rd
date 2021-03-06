% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Utils.R
\name{grid.list}
\alias{grid.list}
\title{\code{\link{expand.grid}} as \code{list}}
\usage{
grid.list(pre, post, ..., print = F)
}
\arguments{
\item{\code{...}:}{set of \code{lists} to be combined (as in \code{\link{expand.grid}})}

\item{\code{pre}:}{prepend this \code{list} to each combination}

\item{\code{post}:}{append this \code{list} to each combination}

\item{\code{print}:}{dump all the sets of combined parameter (N.B.: params in \code{pre} and \code{post} are not shown)}
}
\value{
a \code{list} or \code{list}s
}
\description{
performs combination of all the elements in the set of lists provided via \code{...}
and store each combiantion into a \code{list} (\code{\link{expand.grid}} produce a \code{data.frame} instead, check also \code{\link{grid.df}})
if \code{pre} and \code{post} \code{list}s are provided they are prepended or appended to each combination \code{list}
}

