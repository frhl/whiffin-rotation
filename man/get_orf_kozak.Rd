% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_orf_kozak.R
\name{get_orf_kozak}
\alias{get_orf_kozak}
\title{get open reading frames kozak strength}
\usage{
get_orf_kozak(x, f = function(y) get_orf(y, share_stops = T))
}
\arguments{
\item{x}{string. UTR sequence.}

\item{f}{function to be used that takes two arguments, 
Ideally, it should be \code{get_orf}.}
}
\value{
a value of kozak strength (3 strong, 2 moderate and 1 weak)
}
\description{
find kozak strength of overlapping reading frames.
}
