% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VSeqTools_pileup.R
\name{consensus}
\alias{consensus}
\title{Returns a consensus sequence given a VSeqTools_pileup object}
\usage{
consensus(pu)
}
\arguments{
\item{pu}{a VSeqTools_pileup object returned by BAM_pileup}
}
\value{
A character vector with entries A,C,G,T,d,i and length equal to
positions on reference.

Warning:  if the VSeqTools_pileup object is formed from a BAM file
with gaps, the consensus will skip missed positions.
}
\description{
Given a pileup, a consensus sequence is determined at each position
on the reference.  The consensus is returned as a character vector
with possible entries A, C, G, T, d, i.
}
