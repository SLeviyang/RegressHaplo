
#' Given a bam file, create a read table
#' @param bam_file path to bam file
#' @param out_dir output directory, will be created if needed
#' @param min_freq minimum frequency to consider as variation
#'
#' @return path to read table file
#' @export
bam_to_read_table.RegressHaplo <- function(bam_file, out_dir, min_freq=.01)
{
  if (!dir.exists(out_dir))
    dir.create(out_dir)

  last_out_dir_char <- substr(out_dir, nchar(out_dir), nchar(out_dir))
  if (last_out_dir_char != "/")
    stop("output directory should be terminated with a //")


  pu <- BAM_pileup(bam_file)
  vc <- variant_calls(pu, min_freq=min_freq)[,c("A", "C", "G", "T", "d")]

  nt_pos <- which(apply(vc, 1, function(cr) sum(cr)>1))

  df <- read_table(bam_file, nt_pos=nt_pos)
  read_table_file <- paste(out_dir, "read_table.csv", sep="")
  write.table(df, read_table_file, sep=",", row.names=F)

  return (read_table_file)
}

#' Given a read table, create parameter files for RegressHaplo
#'
#' @param out_dir output directory for parameter files and assumed
#' directory of read table
#'
#' @details The read table is assumed to be in out_dir with filename read_table.csv
#' @return NULL
#' @export
read_table_to_parameters.RegressHaplo <- function(out_dir, min_cover=500,
                                                  create_loci_file=T)
{
  read_table_file <- paste(out_dir, "read_table.csv", sep="")

  df <- read.table(read_table_file, header=T, sep=",", check.names = F,
                   stringsAsFactors = F)
  df <- clean.read_table(df, min_count=20)

  par <- parameters.RegressHaplo(df,
                                 max_global_dim=1200,
                                 max_local_dim=1200,
                                 min_cover=min_cover)

  y <- matrix(par$y, ncol=1)
  P <- par$P
  h <- par$h

  y_file <- paste(out_dir, "y.csv", sep="")
  P_file <- paste(out_dir, "P.csv", sep="")
  h_file <- paste(out_dir, "h.csv", sep="")

  write.table(y, y_file, sep=",", row.names=F, col.names=F)
  write.table(P, P_file, sep=",", row.names=F, col.names=F)
  write.table(h, h_file, sep=",", row.names=F,
              col.names=T)

  if (create_loci_file) {
    loci_file <- paste(out_dir, "loci.csv", sep="")
    con <- file(loci_file, open="w")
    loci_pos <- sapply(par$loci, function(locus) {
      paste(locus, collapse=", ")
    })
    writeLines(loci_pos, con)
    close(con)
  }

  return (NULL)
}

#' Given a parameter files, run RegressHaplo to produce multiple solutions
#'
#' @param out_dir output directory for parameter files and assumed
#' directory of parameter files
#'
#' @details The parameter files assumed to exist are y.csv, P.csv, h.csv
#' @return path to haplotype file
#' @export
parameters_to_solutions.RegressHaplo <- function(out_dir, num_trials=100)
{
  y <- get_y_workflow.RegressHaplo(out_dir)
  P <- get_P_workflow.RegressHaplo(out_dir)

  solutions <- solutions.RegressHaplo(y, P,
                                      num_trials=num_trials,
                                      rho_vals=c(.1,1,2,5),
                                      kk_vals=2)

  sol_file <- file <- paste(out_dir, "solutions.csv", sep="")
  write.table(solutions, sol_file, sep=",",
              row.names = F, col.names = F)

  return (sol_file)
}

#' Given a RegressHaplo solutions file, determine optimal haplotypes
#'
#' @param out_dir output directory for haplotype file and assumed
#' directory of solutions file
#'
#' @return path to solutions file
#' @export
solutions_to_haplotypes.RegressHaplo <- function(out_dir)
{
  s <- get_solutions_workflow.RegressHaplo(out_dir)
  h <- get_h_workflow.RegressHaplo(out_dir)
  rhs <- RegressHaploSolutions(s, h)

  H <- best_fit.RegressHaploSolutions(rhs, K_val = NULL)$H

  H_file <- paste(out_dir, "HAPLO.csv", sep="")
  write.Haplo(H, H_file)

  return (H_file)
}

#' Execute full RegressHaplo workflow a
#' @param bam_file bam file
#' @param out_dir output directory
#' @param min_freq minimum frequency below which variation is ignored
#' @param num_trials number of trials to run in order to find optimal solution
workflow.RegressHaplo <- function(bam_file, out_dir, min_freq=.01,
                                  num_trials=100)
{
  bam_to_read_table.RegressHaplo(bam_file, out_dir, min_freq=min_freq)
  read_table_to_parameters.RegressHaplo(out_dir)
  parameters_to_solutions.RegressHaplo(out_dir, num_trials=num_trials)
  H_file <- solutions_to_haplotypes.RegressHaplo(out_dir)

  return (H_file)
}

#####
# helper functions
get_y_workflow.RegressHaplo <- function(out_dir)
{
  file <- paste(out_dir, "y.csv", sep="")

  y <- read.table(file)
  y <- as.numeric(y[,1])
  return (y)
}

get_P_workflow.RegressHaplo <- function(out_dir)
{
  file <- paste(out_dir, "P.csv", sep="")

  P <- read.table(file, sep=",")
  P <- as.matrix(P)
  colnames(P) <- NULL
  rownames(P) <- NULL

  return (P)
}

get_h_workflow.RegressHaplo <- function(out_dir)
{
  file <- paste(out_dir, "h.csv", sep="")

  h <- read.table(file, sep=",", header=T,
                  colClasses = "character",
                  check.names = F)
  h <- as.matrix(h)

  return (h)
}

get_solutions_workflow.RegressHaplo <- function(out_dir)
{
  sol <- paste(out_dir, "solutions.csv", sep="")
  sol_m <- read.table(sol, sep=",")
  sol_m <- as.matrix(sol_m)
  colnames(sol_m) <- NULL

  return (sol_m)
}

get_HAPLO_workflow.RegressHaplo <- function(out_dir)
{
  hfile <- paste(out_dir, "HAPLO.csv", sep="")

  H <- read.Haplo(hfile)
  return (H)
}
