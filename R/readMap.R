#' readMap constructor
#'
#' A readMap object maps reads to haplotypes, and then
#' determines the sampled and predicted frequencies of the
#' haplotypes at each nucleotide position
#'
#' @param df a read table
#' @param h a haplotype matrix
#' @param pi haplotype frequencies
#'
#' @return A readMap object is a list, with an element of each
#' nucleotide position in the haplotypes/read table.  Within each
#' element is a list containing predicted_freq (pi),
#' sampled_freq (frequency of haplotypes based on read mapping),
#' coverage (number of reads covering the nucleotide position),
#' alleles (nucleotides of haplotypes).  Sampled_freq may not sum
#' to 1 if some reads do not map to a haplotype.  Reads map only
#' if they match haplotype at every read position.
readMap <- function(df, h, pi)
{
  npos <- ncol(h)

  # K x R, assignment frequency matrix
  M_freq <- read_assignment_prob_matrix.readMap(df, pi, h)

  # R x P, does read cover position?
  M_coverage <- reads_covering_positions.read_table(df)

  counts <- df$count
  ncounts <- length(counts)
  M_counts <- sapply(1:ncounts, function(i) {
    as.numeric(M_coverage[i,])*counts[i]
  })
  if (class(M_counts)=="matrix")
    M_counts <- t(M_counts)
  else
    M_counts <- matrix(M_counts, nrow=ncounts, ncol=1)

  info <- lapply(1:npos, function(i) {
    pos_counts <- M_counts[,i]
    coverage <- sum(pos_counts)
    freq <- (M_freq %*% pos_counts)/coverage

    return (list(sampled_freq=freq,
                 predicted_freq=pi,
                 coverage=coverage,
                 alleles=h[,i]))
  })

  names(info) <- colnames(h)
  return (info)
}

get_mapping.readMap <- function(rm, pos)
{
  pos <- as.character(pos)
  ind <- which(names(rm)==pos)

  if (length(ind)==0)
    return (NA)
  return (rm[[ind]])
}

###------ methods to assign reads to haplotypes in naive way
#' Calculates probability of read assignments to haplotypes
#'
#' @param df A read table
#' @param pi Numeric vector of haplotype frequencies
#' @param h haplotype matrix
#'
#' @return A (K x R) matrix with K equaling the number of haplotypes and
#' R equaling the number of read groups (rows in df).
#'
#' @details Entry k,r of the returned matrix is the probability that
#' a read r is assigned to haplotype k.  Some reads
#' may not be assignable, specifically when a read contains a nucleotide
#' at a position for which all phi have probability 0.  In this case, 0
#' is the assigned probability.
#' @export
read_assignment_prob_matrix.readMap <- function(df, pi, h)
{
  df_nucs <- select(df, -count)
  if (ncol(df_nucs) != ncol(h)) {
    stop("number positions in read table doesn't match haplotype lengths")
  }
  K <- length(pi)

  # for each read type, generate assignment prob given by
  #  pi[k] * I(read matches haplotype k)
  # where s is an index over the positions in the read that vary
  # and r_s is the nucleotide at position s of the read
  # if no reads match, return NA
  marginal_prob <- apply(df_nucs, 1, function(nucs) {

    ind <- which(!is.na(nucs))

    if (length(ind)==0)
      return (rep(NA, K))

    nuc_vals <- nucs[ind]

    hap_prob <- apply(h, 1, function(c_haplo) {
      as.numeric(all(nuc_vals == c_haplo[ind]))
    })

    probs <- pi*hap_prob

    # read doesn't match any haplotype
    if (sum(probs)==0)
      return (rep(0, K))

    return (probs/sum(probs))
  })

  if (class(marginal_prob) != "matrix")
    marginal_prob <- matrix(marginal_prob, nrow=1)

  return (marginal_prob)
}

#####################################################################
accuracy.readMap <- function(rm)
{
  acc <- sapply(rm, function(crm) sum(crm$sampled_freq))
  return(acc)
}

plot.readMap <- function(rm, df=NULL, animal=NULL, legend=T)
{
  pos <- names(rm)
  pi <- rm[[1]]$predicted_freq
  npos <- length(pos)
  nhaps <- length(pi)

  # get information for labeling nucs as T/F or D (divergent)
  if (!is.null(animal)) {
    founder_df <- get_founder_table(animal)
    active_pos <- as.numeric(pos)
    founder_df <- filter(founder_df, is.element(pos, active_pos))
  }
  else
    founder_df <- data.frame(founder=rep("X", nhaps))


  plot_df <- Map(function(cmf, pind) {
    alleles <- cmf$alleles
    founder_nuc <- founder_df$founder[pind]
    div <- ifelse(founder_nuc==alleles, "T/F", "D")

    data.frame(xmin=(pind-1), xmax=pind,
               freq=cmf$sampled_freq,
               nuc=alleles,
               div=div,
               hap=1:nhaps,
               pos=pind,
               ntpos=pos[pind],
               stringsAsFactors = F)
  }, rm, 1:npos)

  plot_df <- do.call(rbind, plot_df)

  # split on position, then produce ymin, ymax for rectangle plots
  hap_df <- ddply(plot_df, .(pos), function(pos_df) {
    ymax <-  cumsum(pos_df$freq)
    ymin <- c(0,ymax[-length(ymax)])
    mutate(pos_df, ymax=ymax, ymin=ymin)
  })

  # add rows to represent pi for comparison
  pi_df <- data.frame(xmin=-1, xmax=0, freq=pi,
                      nuc="", div="H",
                      hap=1:nhaps,
                      pos=0,
                      ntpos=NA,
                      ymax=cumsum(pi),
                      ymin=c(0, cumsum(pi)[-length(pi)]),
                      stringsAsFactors = F)

  data_df <- rbind(pi_df, hap_df)

  # create data.frame to split unlinked positions
  if (!is.null(df)) {
    active_pos <- pos
    linkage_ind <- unlinked_pos.read_table(df, min_cover=500)[active_pos]
    linkage_pos <- c(0, (0:(npos-1))[linkage_ind], npos)

    linkage_df <- data.frame(x=linkage_pos,
                             xend=linkage_pos,
                             y=0, yend=1)
  }


  p <- ggplot()

  p <- p + geom_rect(mapping=aes(xmin=xmin, xmax=xmax,
                                 ymin=ymin, ymax=ymax,
                                 fill=div),
                                 data=data_df,
                                 col="black", size=.5)

  p <- p + geom_text(mapping=aes(x=(xmin+xmax)/2, y=(ymax+ymin)/2,
                                 label=nuc),
                     data=data_df)

  if (!is.null(df))
    p <- p + geom_segment(mapping=aes(x=x, y=y,
                                      xend=xend, yend=yend),
                          data=linkage_df, col="white", size=4)


  yave <- (cumsum(pi) + c(0, cumsum(pi)[-length(pi)]))/2
  p <- p + scale_y_continuous(expand = c(0,0),
                              breaks=yave,
                              labels=round(100*pi)/100,
                              limits=c(0, 1))
  p <- p + scale_x_continuous(expand = c(0,0),
                              breaks=seq(from=1/2,
                                         by=1,
                                         length.out=npos),
                              labels=pos)
  p <- p + theme(axis.text=element_text(size=15))
  p <- p + xlab("position") + ylab(NULL)
  if (!legend)
    p <- p + theme(legend.position="none")

  return (p)
}
