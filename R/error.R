###########################################################
# functions to analyze single position frequency errors

#' Determine the frequency error at variable positions
#'
#' @param bam_file path to bam file
#' @param haplotypes character.vector of haplotypes
#' @param freq numeric.vector of haplotype frequencies
#' @param min_coverage minimum coverage needed to compute errors, otherwise
#' an NA is returned
#'
#' @return a data.frame containing columns pos, error that give
#' reference position and frequency error at the position.  Errors are the
#' sum of absolute value of frequency error over the nucleotides.
single_position.error <- function(bam_file, haplotypes, freq,
                                  variable_pos,
                                  minimum_coverage=1000)
{
  pu <- BAM_pileup(bam_file, max_depth=2000)
  hap_m <- plyr::aaply(haplotypes, 1, function(s) strsplit(s, split="")[[1]])

  if (nrow(pu) != ncol(hap_m)) {
    warning("haplotypes and BAM pileup do not have the same positions")
  }

  nhaps <- length(freq)
  npos_hap <- ncol(hap_m)

  pu_nuc_names <- c("A", "C", "G", "T", "d")

  hap_nuc_names <- c("A", "C", "G", "T", "-")
  template <- rep(0, length(pu_nuc_names))
  names(template) <- hap_nuc_names

  pu_pos <- match(pu$pos, variable_pos)
  pu <- pu[!is.na(pu_pos),]

  errors <- plyr::adply(pu, 1, function(cpu) {
    cpos <- cpu$pos
    if (cpos > npos_hap)
      return (data.frame(pos=cpos, error=NA))

    if (cpu$coverage < minimum_coverage)
      return (data.frame(pos=cpos, error=NA))

    hap_freqs <- template
    for (j in 1:nhaps) {
      cnuc <- hap_m[j,cpos]
      hap_freqs[cnuc] <- hap_freqs[cnuc] + freq[j]
    }
    pu_freqs <- as.numeric(cpu[1,pu_nuc_names])
    e <- sum(abs(hap_freqs-pu_freqs))

    return (data.frame(pos=cpos, error=e))
  }, .id=NULL, .expand = F)

  return (errors)
}


#' Plot the position errors
#'
#' @param df data.frame with pos and error columns, both numeric
#' @param boxplot If T, show a boxplot that considers all errors collectivly,
#' otherwise show a barplot giving error at each position
#'
plot_single_position.error <- function(df, boxplot=T)
{
  # remove NA positions
  if (!is.element("label", names(df)))
    df$label <- factor("")

  ind <- apply(df, 1, function(x) !any(is.na(x)))
  df <- df[ind,]

  p <- ggplot()
  if (boxplot)
    p <- p + geom_boxplot(mapping=aes(x=label, y=error),
                          data=df)# position="dodge", stat="identity")
  else
    p <- p + geom_bar(mapping=aes(x=pos, y=error, fill=label),
                      data=df, position="dodge", stat="identity")

  return (p)
}


#################################################################################
# functions for position pair error analysis
create_position_pair_comparison_read_table.error <- function(animal, week)
{
  variable_pos <- variable_positions.error(animal, week, minimum_coverage = 1000)
  bf <- get_BAM_file(animal, week)
  pu <- BAM_pileup(bf)
  vc <- variant_calls(pu, min_freq=0.007)
  vc_ind <- match(variable_pos, vc$pos)
  if (any(is.na(vc_ind)))
    stop("problem creating comparison table", animal, week, "\n")
  vc <- vc[vc_ind,]

  rt <- read_table(bf, variant_calls=vc)
  rt <- clean.read_table(rt, min_count=0, remove_outside_reads = T,
                         remove_empty_cols = T)

  out_dir <- get_output_dir(animal, week)
  rt_file <- paste(out_dir, "comparison_read_table_", animal, "_w", week, ".csv", sep="")
  write.table(rt, rt_file, sep=",", row.names=F)

  return (rt_file)
}

get_position_pair_comparison_read_table.error <- function(animal, week)
{
  out_dir <- get_output_dir(animal, week)
  read_table_file <- paste(out_dir, "comparison_read_table_", animal, "_w", week, ".csv", sep="")

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


#' Get position pairs with threshold coverage
#'
#' @return a data.frame with columns pos1, pos2 giving the positions of the pairs
position_pairs.error <- function(animal, week, minimum_pair_coverage=500)
{
  variable_pos <- variable_positions.error(animal, week, minimum_coverage = 1000)
  rt <- get_position_pair_comparison_read_table.error(animal, week)

  # consider all variable pos pairs with 200 base pair
  var_pos_pairs <- adply(variable_pos, 1, function(pos) {
    near_pos_ind <- which(abs(pos-variable_pos) < 250 & variable_pos > pos)
    pos_pairs <- variable_pos[near_pos_ind]

    if (length(pos_pairs)==0)
      return (NULL)

    return (data.frame(pos1=pos, pos2=pos_pairs))
  }, .id = NULL)

  var_pos_pairs_valid <- ddply(var_pos_pairs, .(pos1, pos2), function(cdf) {
    pos1 <- cdf$pos1; pos2 <- cdf$pos2
    pos1_s <- as.character(pos1); pos2_s <- as.character(pos2)
    rt_sub <- subset.read_table(rt, pos=c(pos1, pos2))

    # remove all rows with an NA, since that only covers one position
    rt_sub <- filter(rt_sub, !is.na(rt_sub[,2]) & !is.na(rt_sub[,3]))

    total_count <- sum(rt_sub$count)
    return (data.frame(count=total_count))

    #  freq <- rt_sub$count/total_count
    #  df <- data.frame(nuc1=rt_sub[,2], nuc2=rt_sub[,3], freq=freq)

    #  return (df)
  }, .drop=T)

  var_pos_pairs_valid <- filter(var_pos_pairs_valid, count >= minimum_pair_coverage)
  return (dplyr::select(var_pos_pairs_valid, -count))
}


#' Extract the position pair frequencies from raw read table
#'
#' @param pair_df data.frame produed by position_pairs.error
true_position_pair_frequencies.error <- function(animal, week, pair_df)
{
  rt <- get_position_pair_comparison_read_table.error(animal, week)

  var_pos_freq <- ddply(pair_df, .(pos1, pos2), function(cdf) {
    pos1 <- cdf$pos1; pos2 <- cdf$pos2
    pos1_s <- as.character(pos1); pos2_s <- as.character(pos2)
    rt_sub <- subset.read_table(rt, pos=c(pos1, pos2))

    # remove all rows with an NA, since that only covers one position
    rt_sub <- filter(rt_sub, !is.na(rt_sub[,2]) & !is.na(rt_sub[,3]))

    total_count <- sum(rt_sub$count)
    freq <- rt_sub$count/total_count
    df <- data.frame(nuc1=rt_sub[,2], nuc2=rt_sub[,3],
                     nuc=paste(rt_sub[,2], rt_sub[,3], sep=""),
                     freq=freq, stringsAsFactors = F)

    return (df)
  }, .drop=T)

  return (var_pos_freq)
}

#' Extract the position pair frequencies from haplotypes
#'
#' @param pair_df data.frame produed by position_pairs.error
haplo_position_pair_frequencies.error <- function(haplo, pair_df)
{
  h <- get_hap.Haplo(haplo)
  pi <- get_freq.Haplo(haplo)

  var_pos_freq <- ddply(pair_df, .(pos1, pos2), function(cdf) {
    pos1 <- cdf$pos1; pos2 <- cdf$pos2
    nuc1 <- h[,pos1]; nuc2 <- h[,pos2]

    df <- data.frame(nuc1=nuc1, nuc2=nuc2, freq=pi)
    df <- ddply(df, .(nuc1, nuc2), function(adf)
      data.frame(nuc=paste(adf$nuc1[1], adf$nuc2[1], sep=""),
                 freq = sum(adf$freq), stringsAsFactors = F),
      .drop=T)

    return (df)
  }, .drop=T)

  return (var_pos_freq)
}

#' Calculate error for each position pair
#'
#' @param df1 A data.frame containing pos1, pos2, nuc1, nuc2, nuc (pasted nuc1,nuc2),
#' freq (produced by true_position_pair_frequencies.error or
#' haplo_position_pair_frequencies.error)
#'
#' @return a data.frame containing pos1, pos2, error where error is sum of absolute
#' frequency differences at each position over true/predicted nucleotide pairs.
calculate_position_pair_error_df.error <- function(df1, df2)
{
  pos_pair_df <- unique(dplyr::select(df1, pos1, pos2))
  pos_pair_df2 <- unique(dplyr::select(df2, pos1, pos2))

  # first check that pairs agree
  if(!all(pos_pair_df, pos_pair_df2))
    stop("position pairs not equal")

  error_df <- ddply(pos_pair_df, .(pos1, pos2), function(cdf) {

    cpos1 <- cdf$pos1[1]; cpos2 <- cdf$pos2[1]

    df1_filt <- dplyr::filter(df1, pos1==cpos1 & pos2==cpos2)
    df2_filt <- dplyr::filter(df2, pos1==cpos1 & pos2==cpos2)

    all_nuc_pairs <- unique(c(df1_filt$nuc, df2_filt$nuc))
    template1 <- rep(0, length(all_nuc_pairs))
    names(template1) <- all_nuc_pairs
    template2 <- template1

    template1[df1_filt$nuc] <- df1_filt$freq
    template2[df2_filt$nuc] <- df2_filt$freq

    cerror <- sum(abs(template1 - template2))
    return (data.frame(error=cerror))
  }, .drop=T)

  return (error_df)
}

pair_position_error_df.error <- function(animal, week)
{
  pair_position_df <- position_pairs.error(animal, week)
  true_pair_df <- true_position_pair_frequencies.error(animal, week,
                                                       pair_position_df)

  H_RH <- get_RegressHaplo_Haplo(animal, week)
  H_PH <- get_PredictHaplo_Haplo(animal, week)

  H_list <- list(H_RH, H_PH)
  names(H_list) <- c("RegressHaplo", "PredictHaplo")

  H_pair_df_list <- lapply(H_list, function(haplo) {
    df <- haplo_position_pair_frequencies.error(haplo,
                                                pair_position_df)
    return (df)
  })

  errors_list <- lapply(H_pair_df_list, function(pair_df) {
    v <- calculate_position_pair_error_df.error(pair_df,
                                                true_pair_df)$error
    return (v)
  })
  names(errors_list) <- names(H_list)

  errors_df <- pair_position_df
  for (i in 1:length(H_list)) {
    cn <- names(H_list)[i]
    errors_df[[cn]] <- errors_list[[cn]]
  }

  return (errors_df)
}

#' Plot the position errors
#'
#' @param boxplot If T, show a boxplot, otherwise barplot
plot_pair_position_error.error <- function(animal, week)
{
  df <- pair_position_error_df.error(animal, week)

  df_m <- melt(df, id.vars =c("pos1", "pos2"))

  p <- ggplot()
  p <- p + geom_boxplot(mapping=aes(x=variable, y=value),
                        data=df_m)

  return (p)
}



#############################################################
get_single_position_error.error_pipeline <- function(out_dir)
{
  file <- paste(out_dir, "single_position_errors.csv", sep="")
  df <- read.csv(file, header=T, colClasses = c("factor", "numeric", "factor"))

  return (df)
}

