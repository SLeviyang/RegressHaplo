% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RegressHaplo_workflow.R
\name{haplotypes_to_parameters.pipeline}
\alias{haplotypes_to_parameters.pipeline}
\title{Given a read table and loci, create parameter files for RegressHaplo.  Parameter
files specify y, P, and the haplotypes represented by h}
\usage{
haplotypes_to_parameters.pipeline(out_dir)
}
\arguments{
\item{out_dir}{output directory for parameter files and assumed
directory of read table}
}
\description{
Given a read table and loci, create parameter files for RegressHaplo.  Parameter
files specify y, P, and the haplotypes represented by h
}
\details{
The read table is assumed to be in out_dir with filename read_table.csv.
}
