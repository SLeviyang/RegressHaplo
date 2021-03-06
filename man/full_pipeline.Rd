% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RegressHaplo_workflow.R
\name{full_pipeline}
\alias{full_pipeline}
\title{Execute full RegressHaplo pipeline}
\usage{
full_pipeline(bam_file, out_dir, max_num_haplotypes = 800,
  min_cover = 500, rho_vals = NULL, start_pos = NULL,
  end_pos = NULL, sig = 0.01, num_trials = 700, heavy_tail = T)
}
\arguments{
\item{bam_file}{bam file}

\item{out_dir}{output directory}

\item{rho_vals}{The values for rho, the penalty parameter, that will be used in
the regression.  If NULL then RegressHaplo will choose values.}

\item{start_pos}{Position on the reference at which the reconstruction begins.}

\item{end_pos}{Position on the reference at which the reconstruction ends}

\item{sig}{The significance level at which variants should be called.}

\item{num_trials}{number of trials to run in order to find optimal solution.
This number of trials is run for each rho value, unless rho=NULL, in which
case this is the total number of trials.}

\item{heavy_tail}{See comments in help of bam_to_variant_calls.pipeline}

\item{max_num_haplotype}{The maximum number of haplotype over which the regression
will be performed.  Lower values mean faster run times but poorer inference.
Any number above 1200 will lead to very slow run times.}
}
\description{
Execute full RegressHaplo pipeline
}
