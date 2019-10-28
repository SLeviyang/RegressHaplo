
#' Given a bam file, create variant call file
#'
#' @param bam_file path to bam file
#' @param out_dir output directory, will be created if needed.  variant call file
#' will be placed in the output directory
#' @param start_pos position on reference at which to start looking for variant
#' calls
#' @param end_pos position on reference at which to end looking for variant calls
#' @param sig significance level at which variants will be called.  An automatic
#' bonferroni correction is applied to account for the number of positions in the
#' reference.
#' @param heavy_tail If T, then a betabinomial is used to call errors following Gerstung et al
#' 2011. If F, then a poisson is used to call variants following Wang et al 2007
#'
#' @return path to variant call file
#'
#' @references Wang,C. et al. (2007) Characterization of mutation
#' spectra with ultra-deep pyrosequencing: application to HIV-1 drug
#' resistance. Genome Res., 17, 1195–1201,
#' @references Gerstung et al. (2011) Reliable detection of subclonal
#' single-nucleotide variants in tumour cell populations.  Nature
#' Communications.
#'
#' @details variant calls are given by a data.frame with 6 columns
#' named c("pos", "A", "C", "G", "T", "-") pos gives the position at which a
#' true variant exists.  The other columns have logical entries,
#' with TRUE meaning that a true variant exists with the corresponding nucleotide.
#' @export
bam_to_variant_calls.pipeline <- function(bam_file, out_dir,
                                          start_pos=NULL, end_pos=NULL,
                                          sig=.01, heavy_tail=T)
{
  out_dir <- fix_out_dir(out_dir)

  pu <- BAM_pileup(bam_file, max_depth=5000, min_base_quality=0,
                   min_mapq=0)
  pu_file <- paste(out_dir, "pile_up.csv", sep="")
  write.csv(pu, pu_file, row.names = F)

  #bonferroni correction
  bf_sig <- sig/nrow(pu)

  # if sig==0, no error correction
  if (sig==0)
    vc <- variant_calls(pu, min_freq=1E-4)
  else
    vc <- variant_calls(pu, sig=bf_sig, heavy_tail=heavy_tail)

  vc_pos <- get_variant_call_pos(vc)
  if (!is.null(start_pos))
    vc_pos <- vc_pos[vc_pos >= start_pos]
  if (!is.null(end_pos))
    vc_pos <- vc_pos[vc_pos <= end_pos]

  vc <- vc[is.element(vc$pos, vc_pos),]
  names(vc) <- c("pos", "A", "C", "G", "T", "-", "i")

  variant_call_file <- paste(out_dir, "variant_calls.csv", sep="")
  write.csv(vc, variant_call_file, row.names = F)

  return (variant_call_file)
}


#' Given variant calls and a bam file, create a read table
#' @param bam_file path to bam file
#' @param out_dir output directory, will be created if needed
#' @param use_raw_read_table If there is an existing raw_read_table.csv file, should
#' it be used as a raw read table?  See details.
#' @param sig significance levels at which to filter reads.
#'
#' @details Converts bam_file to a raw read table, written to
#' $outdir/raw_read_table.csv and then applies error correction
#' to reads and creates a filtered read table, written to
#' $outdir/read_table.csv.   The parameter use_raw_read_table can
#' be used to skip the step of creating the raw read table, which
#' is very time consuming.  This is particularly useful, if the
#' user wants to experiment with different error levels to form
#' an acceptable read table.
#'
#' @return path to read table file
#' @export
variant_calls_to_read_table.pipeline <- function(bam_file,
                                       out_dir,
                                       use_raw_read_table=F,
                                       sig=.01,
                                       debug=F)
{
  out_dir <- fix_out_dir(out_dir)
  variant_calls <- get_variant_calls.pipeline(out_dir)
  pu <- get_pile_up.pipeline(out_dir)

  # construct raw read table
  if (!use_raw_read_table) {
    df <- read_table(bam_file, variant_calls=variant_calls, pu=pu,
                   debug=debug)
    read_table_file <- paste(out_dir, "raw_read_table.csv", sep="")
    write.table(df, read_table_file, sep=",", row.names=F)
  } else
    df <- get_read_table.pipeline(out_dir, raw=T)

  # clean up read table to eliminate empty rows
  print("cleaning read table")
  df <- clean.read_table(df, min_count=0,
                         remove_outside_reads=T,
                         remove_empty_cols=T,
                         remove_non_variant_pos = F,
                         remove_deletions = F,
                         remove_partial_cover_reads = F)

  # filter reads based on Poisson error model to construct read table
  # and minimum threshold for each read partition
  print("setting up for read table error correction")
  pu <- get_pile_up.pipeline(out_dir)
  error_rate <- get_error_rate(pu, split_by_nuc=F)
  
  # if there's an error, filter with it
  if (error_rate > 1E-4)
    df_filter <- error_filter.read_table(df, error_freq=error_rate, sig=sig)
  else
    df_filter <- df

  print("cleaning read table again")
  df_filter <- clean.read_table(df_filter, min_count=0,
                         remove_outside_reads=T,
                         remove_empty_cols=T,
                         remove_non_variant_pos = F,
                         remove_deletions = F,
                         remove_partial_cover_reads = F)

  read_table_file <- paste(out_dir, "read_table.csv", sep="")
  write.table(df_filter, read_table_file, sep=",", row.names=F)

  return (read_table_file)
}

#' Given a read table, split variable positions into loci
#'
#' @param out_dir output directory for parameter files and assumed
#' directory of read table
#' @param max_num_haplotypes the maximum number of haplotypes over the full reference
#' that will be considered.
#'
#' @details The read table is assumed to be in out_dir with filename read_table.csv.  The
#' parameter min_cover can (and should) be used to force loci to break across low coverage
#' regions of the reference.
#' @return NULL
#' @export
read_table_to_loci.pipeline <- function(out_dir,
                                        max_num_haplotypes=1200)
{
  df <- get_read_table.pipeline(out_dir)

  # split into loci based on gaps
  print("making initial locus choices")
  sdf_initial <- split_unlinked_loci.read_table(df, sort_loci=F,
                                                min_cover=0)

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
      # if there are no consistent haplotypes,
      # hp can come back as NULL or nrow==0, this needs to be corrected
      # in consistent haplotypes computations, so that NULL is only return,
      # but for now we catch it here
      if (!is.null(hp))
        if (nrow(hp)==0)
          hp <- NULL

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
    cat("processed:", length(sdf_final), "unprocessed:", length(sdf), "\n")
  } # end while loop

  # we need to sort sdf_final since the loci might not be in order
  first_pos <- sapply(sdf_final, function(cdf) as.numeric(pos_names.read_table(cdf)[1]))
  sdf_final <- sdf_final[order(first_pos)]

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
                                        max_num_haplotypes=1000)
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

  cat("initial number of global haplotypes:", nhaps, "\n")
  # if filtering is required then try to hit target within .8-1 of max_num_haplotypes
  if (nhaps > max_num_haplotypes) {
    print("too many global haplotypes, reducing number of local haplotypes!")
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
      cat("left_n right_n center_rho", left_n, right_n, center_rho, "\n")

      if (ratio >= .8 & ratio <= 1) {
        break
      } else if (center_n > max_num_haplotypes) {
        left_rho <- center_rho
        left_n <- center_n
      }
      else {
        right_rho <- center_rho
        right_n <- center_n
      }


      counter <- counter + 1
      if (counter == 20)
        break
    }
    # if we exited due to counter, choose rho as right_rho
    # to fall below the max_num_haplotypes target
    if (counter==20)
      local_rho <- right_rho
    else
      local_rho <- center_rho
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

  par <- penalized_regression_parameters.RegressHaplo(df, h, position_fit=F)

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
#' @param num_trials number of starting points to try in attempting to find
#' a global min for each rho values
#' @param rho_values rho values to use in finding global min.
#'
#' @details The parameter files assumed to exist are y.csv, P.csv, h.csv
#' @return path to haplotype file
#' @export
parameters_to_solutions.pipeline <- function(out_dir, num_trials=700,
                                             rho_vals=c(.1,1,5,10,20))
{
  # if rho vals are null, try to find rho values that capture
  # large ran and produce different fit values.
  if (is.null(rho_vals)) {
    cat("Since RHO values have not been provided, running 50 regressions to determine rho values.\n")
    rho_vals <- exp(seq(log(.001), log(10), length.out=50))
    parameters_to_solutions.pipeline(out_dir, num_trials=1,
                                     rho_vals = rho_vals)
    m <- get_solutions.pipeline(out_dir)
    fits <- round(10^5*m[1,])/10^5

    df_rho <- data.frame(rho=rho_vals, fit=fits)
    rho_final <- plyr::daply(df_rho, .(fits), function(cdf) {
      mean(cdf$rho)
    })

  } else
    rho_final <- rho_vals


  num_trials_per_rho <- ceiling(num_trials/length(rho_final))

  y <- get_y.pipeline(out_dir)
  P <- get_P.pipeline(out_dir)

  solutions <- solutions.RegressHaplo(y, P,
                                      num_trials=num_trials_per_rho,
                                      rho_vals=rho_final,
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
#' @param K_val pick the best solution out of solutions with K_val
#' haplotypes.  If NULL, then optimal K_val is used.
#'
#' @return path to solutions file
#' @export
solutions_to_haplotypes.pipeline <- function(out_dir, K_val=NULL)
{
  s <- get_solutions.pipeline(out_dir)
  h <- get_h.pipeline(out_dir)

  rhs <- RegressHaploSolutions(s, h)

  H <- best_fit.RegressHaploSolutions(rhs, K_val = K_val)$H

  out_dir <- fix_out_dir(out_dir)
  H_file <- paste(out_dir, "final_haplo.csv", sep="")
  write.Haplo(H, H_file)

  return (H_file)
}

#' Given RegressHaplo haplotypes file, produce fasta file showing
#' haplotypes over the full reference with frequencies as part of
#' sequence names
#'
#' @param bam_file bam file containing reads for haplotype reconstruction
#' @param out_dir output directory for haplotype file and assumed
#' directory of solutions file
#'
#' @return path to fasta file
#' @export
haplotypes_to_fasta.pipeline <- function(bam_file, out_dir)
{
  haplo <- get_haplo.pipeline(out_dir)
  haplo_m <- get_hap.Haplo(haplo)
  pi <- get_freq.Haplo(haplo)
  nhaps <- length(pi)

  variable_pos <- get_variable_positions.pipeline(out_dir)
  colnames(haplo_m) <- variable_pos

  pu <- get_pile_up.pipeline(out_dir)

  consensus <- consensus(pu)
  all_pos <- 1:nrow(pu)

  out_haps <- sapply(all_pos, function(i) {
    if (is.element(i, variable_pos))
      return (haplo_m[,as.character(i)])
    else
      return (rep(consensus[i], nhaps))
  })
  hap_s <- apply(out_haps, 1, paste, collapse="")
  dna <- BStringSet(hap_s)
  names(dna) <- paste("haplotype", 1:nhaps, "_", round(10^4*pi)/10^4, sep="")

  out_dir <- fix_out_dir(out_dir)
  fasta_file <- paste(out_dir, "final_haplo.fasta", sep="")
  writeXStringSet(dna, fasta_file)

  return (fasta_file)
}

#' Execute full RegressHaplo pipeline
#' @param bam_file bam file
#' @param out_dir output directory
#' @param max_num_haplotype The maximum number of haplotype over which the regression
#' will be performed.  Lower values mean faster run times but poorer inference.
#' Any number above 1200 will lead to very slow run times.
#' @param rho_vals The values for rho, the penalty parameter, that will be used in
#' the regression.  If NULL then RegressHaplo will choose values.
#' @param start_pos Position on the reference at which the reconstruction begins.
#' @param end_pos Position on the reference at which the reconstruction ends
#' @param sig The significance level at which variants should be called.
#' @param num_trials number of trials to run in order to find optimal solution.
#' This number of trials is run for each rho value, unless rho=NULL, in which
#' case this is the total number of trials.
#' @param heavy_tail See comments in help of bam_to_variant_calls.pipeline
#'
#' @export
full_pipeline <- function(bam_file, out_dir,
                          max_num_haplotypes=800,
                          rho_vals=NULL,
                          start_pos=NULL, end_pos=NULL,
                          sig=.01, num_trials=700, heavy_tail=T)
{
  if (!is.null(start_pos) & !is.null(end_pos))
    if (end_pos < start_pos) {
      cat("REGRESSHAPLO WARNING:\n",
          "End position of reconstruction must be greater than start position!\n",
          "Pipeline has been stopped.\n")
      return (NULL)
    }

  if (!dir.exists(out_dir)) {
    cat("REGRESSHAPLO WARNING:\n")
    cat(out_dir, "does not exist\n")
    cat("Please create the output directory and then call full_pipeline again\n")
    return (NULL)
  }
  cat("Making variant calls...\n")
  bam_to_variant_calls.pipeline(bam_file, out_dir,
                                start_pos=start_pos, end_pos=end_pos,
                                sig=sig, heavy_tail=heavy_tail)

  # check if there are variant calls
  variant_calls <- get_variant_calls.pipeline(out_dir)
  if (nrow(variant_calls)==0) {
    cat("REGRESSHAPLO WARNING:", "There are no variable positions called!\n",
        "The region is homogeneous up to NGS calls.\n",
        "Pipeline has been stopped!\n")
    return (NULL)
  }

  cat("Constructing read table (this may take a while)...\n")
  variant_calls_to_read_table.pipeline(bam_file, out_dir, sig=sig)

  cat("Constructing regions/loci...\n")
  read_table_to_loci.pipeline(out_dir, max_num_haplotypes=max_num_haplotypes)

  cat("Constructing haplotypes...\n")
  loci_to_haplotypes.pipeline(out_dir, max_num_haplotypes=max_num_haplotypes)

  cat("Preparing matrices and vectors for regression...\n")
  haplotypes_to_parameters.pipeline(out_dir)

  cat("Solving regressions..\n")
  parameters_to_solutions.pipeline(out_dir, num_trials=num_trials,
                                               rho_vals=rho_vals)


  cat("Choosing the best regression solution for the reconstruction..\n")
  solutions_to_haplotypes.pipeline(out_dir)
  haplotypes_to_fasta.pipeline(bam_file, out_dir)

  cat("REGRESSHAPLO reconstruction complete!\n")
  cat("See final_haplo.fasta in", out_dir, "directory for the reconstruction \n")

  return (NULL)
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
get_pile_up.pipeline <- function(out_dir)
{
  out_dir <- fix_out_dir(out_dir)
  file <- paste(out_dir, "pile_up.csv", sep="")

  df <- read_BAM_pileup(file)
  
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

  # read once to find out how many columns.  This needs to be done because
  # read.table is interpreting "T" and "TRUE" in cases when the column is all
  # "T"
  df_h <- read.table(read_table_file, header=T, sep=",", check.names = F,
                     nrow=1,
                   stringsAsFactors = F)
  colcl <- c("integer", rep("character", ncol(df_h)-1))

  df <- read.table(read_table_file, header=T, sep=",", check.names = F,
                   colClasses = colcl,
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

  # BUG FIX HERE BY Stephen Shank and Dave Bouvier
  #loci <- strsplit(as.character(loci_df$pos), split="\\+")
  loci <- strsplit(as.character(loci_df$pos), split="\\+")

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

#' @export
get_fasta.pipeline <- function(out_dir)
{
  out_dir <- fix_out_dir(out_dir)
  hfile <- paste(out_dir, "final_haplo.fasta", sep="")

  dna <- readBStringSet(hfile)
  freq_names <- names(dna)
  freq <- sapply(freq_names, function(s) strsplit(s, split="_")[[1]][2])
  freq <- as.numeric(freq)

  return (list(haplotypes=as.character(dna), freq=freq))
}

#' @export
get_variable_positions.pipeline <- function(out_dir)
{
  h <- get_h.pipeline(out_dir)
  pos <- as.numeric(colnames(h))

  return (pos)
}


#' Add a / to out_dir if it is not there
fix_out_dir <- function(out_dir)
{
  last_out_dir_char <- substr(out_dir, nchar(out_dir), nchar(out_dir))
  if (last_out_dir_char != "/")
    out_dir <- paste(out_dir, "/", sep="")

  return (out_dir)
}


############################################################
# methods to manipulate solutions file

#' Summary of regression fits over all starting points
#'
#' @param out_dir RegressHaplo pipeline directory
#'
#' @return A data.frame with columns (rho) rho value for fitting,
#' (K) number of haplotypes reconstruted,
#' (fit) fit of haplotype reconstruction, (solution_number).  Each
#' row of the data.frame corresponds to a solution of the optimization
#' for a given starting point and rho value.
#' @export
get_solutions_summary.pipeline <- function(out_dir)
{
  s <- get_solutions.pipeline(out_dir)
  h <- get_h.pipeline(out_dir)

  rhs <- RegressHaploSolutions(s, h)
  df <- dplyr::select(rhs$df_stats, -kk)

  return (df)
}

#' Haplotype reconstruction of given solution
#'
#' @param out_dir, RegressHaplo pipeline directory
#' @param i solution number
#'
#' @return a data.frame with columns (haplotype) haplotypes reconstructed, (freq)
#' frequencies of haplotype, and the haplotype number from the list of global
#' haplotypes in h.csv (haplotype_number)
#' @export
get_solutions_haplotype_reconstruction.pipeline <- function(out_dir, i)
{
  s <- get_solutions.pipeline(out_dir)
  h <- get_h.pipeline(out_dir)

  if (ncol(s) < i) {
    cat("you are asking for solution", i, "but only", ncol(s), "solutions exist", "\n")
    stop("stopping!")
  }

  freqs <- s[5:nrow(s),i]
  hap_ind <- which(freqs > 0)
  hap <- apply(h[hap_ind,,drop=F], 1, paste, collapse="")

  df <- data.frame(haplotype=hap, freq=freqs[hap_ind],
                   haplotype_number=hap_ind, stringsAsFactors = F)

  return (df)
}
