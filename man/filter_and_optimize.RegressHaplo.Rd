% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RegressHaplo_util.R
\name{filter_and_optimize.RegressHaplo}
\alias{filter_and_optimize.RegressHaplo}
\title{Generate consistent haplotypes for a read table and, if desired, apply
RegressHaplo optimization.}
\usage{
filter_and_optimize.RegressHaplo(
  df,
  global_rho = NULL,
  max_global_dim = 10000,
  max_local_dim = 1200,
  min_cover = 500,
  run_optimization = F
)
}
\arguments{
\item{df}{read table}

\item{global_rho}{If a global fit should be computed, the rho that should
be used.}

\item{max_local_dim}{The maximum number of haplotypes that can be filtered}

\item{min_cover}{The minimum read coverage needed to link across a
read table position}

\item{run_optimization}{Should an optimization be run, or should just
consistent global haplotypes be returned.}

\item{max_gobal_dim}{The maximum number of consistent haplotypes that
should be generated}
}
\value{
A list constaining the elements df, pi, fit, and h.   df
is simply the read_table returned.  h are the global
consistent haplotypes generated after filtering; h is a
character matrix with colnames giving positions. pi and fit
are NA if the optimization is not run, otherwhise pi is a vector
of frequencies with length equal to the number of haplotypes (nrow(h))
and fit is a scalar describing the fit of the solution.
}
\description{
Generates consistent haplotypes by filtering local haplotypes using
the RegressHaplo algorithm to satisfy dimensions requirements, and
if desired then applies the RegressHaplo algorithm globally.
}
\details{
Haplotypes are generated by splitting the read table
positions into loci and then iteratively filtering the local haplotypes
using the RegressHaplo algorithm until all combinations of local haplotypes
have dimension less than max_global_dim.  At the outset, loci are defined
as positions linked by reads, but if a locus has too many consistent
haplotypes (> max_local_dim), then a locus is split in half until the
dimension is reduced.  This allows application of the RegressHaplo
algorithm locally.

To run an optimization, run_optimization must be TRUE and a global_rho
must be provided.
}
