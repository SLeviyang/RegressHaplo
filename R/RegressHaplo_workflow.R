
#' Given a bam file, create variant call file
#'
#' @param bam_file path to bam file
#' @param out_dir output directory, will be created if needed.  variant call file
#' will be placed in the output directory
#'
#' @return path to variant call file
#'
#' @details variant calls are given by a data.frame with 6 columns
#' named c("pos", "A", "C", "G", "T", "-") pos gives the position at which a
#' true variant exists.  The other columns have logical entries,
#' with TRUE meaning that a true variant exists with the corresponding nucleotide.
#' @export
bam_to_variant_calls.pipeline <- function(bam_file, out_dir,
                                          sig=.01, heavy_tail=T)
{
  out_dir <- fix_out_dir(out_dir)

  pu <- BAM_pileup(bam_file, max_depth=5000, min_base_quality=0,
                   min_mapq=0)

  #bonferroni correction
  bf_sig <- sig/nrow(pu)

  vc <- variant_calls(pu, sig=bf_sig, heavy_tail=heavy_tail)
  vc_pos <- get_variant_call_pos(vc)

  vc <- vc[vc_pos,]
  names(vc) <- c("pos", "A", "C", "G", "T", "-", "i")

  variant_call_file <- paste(out_dir, "variant_calls.csv", sep="")
  write.csv(vc, variant_call_file, row.names = F)

  return (variant_call_file)
}


#' Given a bam file, create a read table
#' @param bam_file path to bam file
#' @param out_dir output directory, will be created if needed
#'
#' @details Converts bam_file to a raw read table, written to
#' $outdir/raw_read_table.csv and then applies error correction
#' to reads and creates a filtered read table, written to
#' $outdir/read_table.csv
#'
#' @return path to read table file
#' @export
bam_to_read_table.pipeline <- function(bam_file,
                                       out_dir)
{
  out_dir <- fix_out_dir(out_dir)
  variant_calls <- get_variant_calls.pipeline(out_dir)

  # construct raw read table
  df <- read_table(bam_file, variant_calls=variant_calls)
  read_table_file <- paste(out_dir, "raw_read_table.csv", sep="")
  write.table(df, read_table_file, sep=",", row.names=F)

  # clean up read table to eliminate empty rows
  df <- clean.read_table(df, min_count=0,
                         remove_outside_reads=T,
                         remove_empty_cols=F,
                         remove_non_variant_pos = F,
                         remove_deletions = F,
                         remove_partial_cover_reads = F)

  # filter reads based on Poisson error model to construct read table
  df_filter <- error_filter.read_table(df, .007, .01/600)
  read_table_file <- paste(out_dir, "read_table.csv", sep="")
  write.table(df_filter, read_table_file, sep=",", row.names=F)

  return (read_table_file)
}

#' Given a read table, split variable positions into loci
#'
#' @param out_dir output directory for parameter files and assumed
#' directory of read table
#'
#' @details The read table is assumed to be in out_dir with filename read_table.csv.
#' @return NULL
#' @export
read_table_to_loci.pipeline <- function(out_dir, min_cover=500,
                                        max_num_haplotypes=1200)
{
  df <- get_read_table.pipeline(out_dir)

  # split into loci based on gaps
  sdf_initial <- split_unlinked_loci.read_table(df, sort_loci=F,
                                                min_cover=min_cover)

  # check each locus to see if the number of haplotypes is less
  # than max_num_haplotypes.  if less, then remove from sdf and
  # add to sdf_final, otherwise split
  sdf_final <- list()
  sdf <- sdf_initial
  while (length(sdf) > 0) {
    sdf_next <- list()
    for (i in 1:length(sdf)) {
      cdf <- sdf[[i]]

      hp <- consistent_haplotypes.read_table(cdf, rm.na=T,
                                             max_num_haplotypes=max_num_haplotypes)
      # hp will be null if number of haplotypes exceeds max_num_haplotypes
      if (!is.null(hp)) {
        cat("accepted locus ", i, "with num haps:", length(hp), "\n")
        sdf_final <- append(sdf_final, list(cdf))
      }
      else {
        cat("splitting locus ", i, "\n")
        split_cdf <- split.read_table(cdf)
        sdf_next <- append(sdf_next, split_cdf)
      }
    } # end for loop
    sdf <- sdf_next
    cat("processed:", length(sdf_final), "tobe:", length(sdf), "\n")
  } # end while loop

  loci <- lapply(sdf_final, function(cdf) {
    as.numeric(names(dplyr::select(cdf, -count)))
  })

  save_loci.pipeline(loci, out_dir)

  return (NULL)
}

#' Given loci, create global haplotypes
#'
#' @param out_dir output directory for parameter files and assumed
#' directory of read table
#'
#' @return NULL
#' @export
loci_to_haplotypes.pipeline <- function(out_dir,
                                        max_num_haplotypes=1200)
{
  df <- get_read_table.pipeline(out_dir)
  loci <- get_loci.pipeline(out_dir)

  # create loci based read tables
  ldf <- lapply(loci, function(cpos) {
    cdf <- df[,c("count", cpos)]
    cdf <- clean.read_table(cdf, min_count=0, remove_empty_cols=F)
    cdf <- regroup.read_table(cdf)
  })

  # first try no filtering
  local_rho <- 0
  nhaps <- number_global_haplotypes.RegressHaplo(local_rho, ldf, max_num_haplotypes)

  # if filtering is required then try to hit target within .8-1 of max_num_haplotypes
  if (nhaps > max_num_haplotypes) {
    # setup bisection
    left_rho <- 0
    left_n <- nhaps

    # make right less than .8*max_num_haplotypes
    right_rho <- 1
    repeat {
      right_n <- number_global_haplotypes.RegressHaplo(right_rho, ldf, max_num_haplotypes)
      if (right_n > .8*max_num_haplotypes)
        right_rho <- right_rho + 2
      else
        break
    }

    # just to be sure let's not allow more than 20 iterations
    counter <- 0
    repeat {
      center_rho <- (left_rho + right_rho)/2
      center_n <- number_global_haplotypes.RegressHaplo(center_rho, ldf, max_num_haplotypes)
      ratio <- center_n/max_num_haplotypes
      if (ratio >= .8 & ratio <= 1) {
        local_rho <- center_rho
        break
      } else if (center_n > max_num_haplotypes) {
        left_rho <- center_rho
        left_n <- center_n
      }
      else {
        right_rho <- center_rho
        right_n <- center_n
      }
      cat("left_n right_n", left_n, right_n, "\n")

      counter <- counter + 1
      if (counter == 20)
        break
    }
  }

  # use local rho to find haplotypes
  h_local <- Map(function(cdf, i) {

      haps <- consistent_haplotypes.read_table(cdf, rm.na=T,
                                               max_num_haplotypes = max_num_haplotypes)
      if (is.null(haps)) {
        cat("locus ", i, "has read table with too many paths\n")
        cat("read_table_to_loci.pipeline was not run.  If it was, please open issue in github.")
        stop("bad locus")
      }
      par <- penalized_regression_parameters.RegressHaplo(cdf, haps)
      rh <- penalized_regression.RegressHaplo(par$y, par$P, rho=local_rho, kk=2, verbose=F)
      pi <- get_pi.RegressHaplo(rh)

      filtered_haps <- haps[pi>0,,drop=F]

      return (filtered_haps)
  }, ldf, 1:length(ldf))

  # form all permutations
  h_consistent <- haplotype_permute.RegressHaplo(h_local)
  cat("number of global haplotypes", nrow(h_consistent), "\n")

  colnames(h_consistent) <- pos_names.read_table(df)

  out_dir <- fix_out_dir(out_dir)
  h_file <- paste(out_dir, "h.csv", sep="")
  write.table(h_consistent, h_file, sep=",", row.names=F, col.names=T)

  return (NULL)
}

#' Given a read table and loci, create parameter files for RegressHaplo.  Parameter
#' files specify y, P, and the haplotypes represented by h
#'
#' @param out_dir output directory for parameter files and assumed
#' directory of read table
#'
#' @details The read table is assumed to be in out_dir with filename read_table.csv.
#' @return NULL
#' @export
haplotypes_to_parameters.pipeline <- function(out_dir)
{
  df <- get_read_table.pipeline(out_dir)
  h <- get_h.pipeline(out_dir)

  par <- penalized_regression_parameters.RegressHaplo(df, h)

  y <- matrix(par$y, ncol=1)
  P <- par$P

  out_dir <- fix_out_dir(out_dir)
  y_file <- paste(out_dir, "y.csv", sep="")
  P_file <- paste(out_dir, "P.csv", sep="")

  write.table(y, y_file, sep=",", row.names=F, col.names=F)
  write.table(P, P_file, sep=",", row.names=F, col.names=F)

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
parameters_to_solutions.pipeline <- function(out_dir, num_trials=100)
{
  y <- get_y.pipeline(out_dir)
  P <- get_P.pipeline(out_dir)

  solutions <- solutions.RegressHaplo(y, P,
                                      num_trials=num_trials,
                                      rho_vals=c(.1,1,5,10,20),
                                      kk_vals=2)

  out_dir <- fix_out_dir(out_dir)
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
solutions_to_haplotypes.pipeline <- function(out_dir)
{
  s <- get_solutions.pipeline(out_dir)
  h <- get_h.pipeline(out_dir)

  rhs <- RegressHaploSolutions(s, h)

  H <- best_fit.RegressHaploSolutions(rhs, K_val = NULL)$H
  browser()

  out_dir <- fix_out_dir(out_dir)
  H_file <- paste(out_dir, "final_haplo.csv", sep="")
  write.Haplo(H, H_file)

  return (H_file)
}

#' Execute full RegressHaplo pipeline
#' @param bam_file bam file
#' @param out_dir output directory
#' @param variant_calls.  A data.frame specifying positions and nucleotides on
#' reference to be taken as true variants, see details.
#' @param num_trials number of trials to run in order to find optimal solution
#'
#' @export
full_pipeline <- function(bam_file, out_dir, variant_calls,
                                  num_trials=100)
{
  bam_to_read_table.pipeline(bam_file, variant_calls=variant_calls, out_dir)
  read_table_to_parameters.pipeline(out_dir)
  parameters_to_solutions.pipeline(out_dir, num_trials=num_trials)
  H_file <- solutions_to_haplotypes.pipeline(out_dir)

  return (H_file)
}

#####
# functions to access pipeline data

#' @export
get_variant_calls.pipeline <- function(out_dir)
{
  out_dir <- fix_out_dir(out_dir)
  file <- paste(out_dir, "variant_calls.csv", sep="")

  df <- read.table(file, header=T, sep=",", check.names = F,
                   stringsAsFactors = F)
  return (df)
}

#' @export
get_read_table.pipeline <- function(out_dir, raw=F)
{
  out_dir <- fix_out_dir(out_dir)
  if (raw)
    read_table_file <- paste(out_dir, "raw_read_table.csv", sep="")
  else
    read_table_file <- paste(out_dir, "read_table.csv", sep="")

  df <- read.table(read_table_file, header=T, sep=",", check.names = F,
                   stringsAsFactors = F)
  return (df)
}

save_loci.pipeline <- function(loci, out_dir)
{
  out_dir <- fix_out_dir(out_dir)

  loci_file <- paste(out_dir, "loci.csv", sep="")
  loci_pos <- sapply(loci, function(locus) {
    paste(locus, collapse="+")
  })
  loci_df <- data.frame(locus=1:length(loci_pos), pos=loci_pos)
  write.csv(loci_df, loci_file, row.names = F)

  return (NULL)
}

#' @export
get_loci.pipeline <- function(out_dir)
{
  out_dir <- fix_out_dir(out_dir)

  loci_file <- paste(out_dir, "loci.csv", sep="")
  loci_df <- read.csv(loci_file, header=T, stringsAsFactors = F)

  loci <- strsplit(loci_df$pos, split="\\+")

  return (loci)
}

#' @export
get_y.pipeline <- function(out_dir)
{
  out_dir <- fix_out_dir(out_dir)
  file <- paste(out_dir, "y.csv", sep="")

  y <- read.table(file)
  y <- as.numeric(y[,1])
  return (y)
}

#' @export
get_P.pipeline <- function(out_dir)
{
  out_dir <- fix_out_dir(out_dir)
  file <- paste(out_dir, "P.csv", sep="")

  P <- read.table(file, sep=",")
  P <- as.matrix(P)
  colnames(P) <- NULL
  rownames(P) <- NULL

  return (P)
}

#' @export
get_h.pipeline <- function(out_dir)
{
  out_dir <- fix_out_dir(out_dir)
  file <- paste(out_dir, "h.csv", sep="")

  h <- read.table(file, sep=",", header=T,
                  colClasses = "character",
                  check.names = F)
  h <- as.matrix(h)

  return (h)
}

#' @export
get_solutions.pipeline <- function(out_dir)
{
  out_dir <- fix_out_dir(out_dir)
  sol <- paste(out_dir, "solutions.csv", sep="")
  sol_m <- read.table(sol, sep=",")
  sol_m <- as.matrix(sol_m)
  colnames(sol_m) <- NULL

  return (sol_m)
}

#' @export
get_haplo.pipeline <- function(out_dir)
{
  out_dir <- fix_out_dir(out_dir)
  hfile <- paste(out_dir, "final_haplo.csv", sep="")

  H <- read.Haplo(hfile)
  return (H)
}


#' Add a / to out_dir if it is not there
fix_out_dir <- function(out_dir)
{
  last_out_dir_char <- substr(out_dir, nchar(out_dir), nchar(out_dir))
  if (last_out_dir_char != "/")
    out_dir <- paste(out_dir, "/", sep="")

  return (out_dir)
}
