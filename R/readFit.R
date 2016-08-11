#' readFit constructor
#'
#' @param df read table
#' @param h haplotype matrix
#' @param pi haplotype frequencies, may be NULL
#' @param position If T, assess fits based on nucleotide
#' positions. If F, assess fits based on reads.
#'
#' @details A readFit object is a list.  Each entry in the list
#' contains the entries sampled_freq, predicted_freq, coverage,
#' alleles, P.   sampled_freq and predicted_freq are numeric vectors
#' with an entry for each allele.  Alleles is a character vector.
#' sampled_freq always sums to 1, since alleles are those seen in the
#' data, while predicted_freq may sum to less than 1 when a haplotype
#' does not match any allele.
#'
#' readFit objects are built based on nucleotide positions or reads.
#' For position based readFits, the alleles are some subset of A,C,G,T
#' which are sampled at a given position.  Coverage is the number of
#' reads that cover the position.  For read based fits, the alleles
#' are nucleotide sequences corresponding to a given sequence template
#' (see tempate in read_table objects).  Coverage is the number of reads
#' with a given template.
#'
#' @return a readFit object
#' @export
readFit <- function(df, h, pi=NULL, position=F)
{
  if (position)
    rf <- position_fit.readFit(df, h, pi=pi)
  else
    rf <- read_fit.readFit(df, h, pi=pi)

  return (rf)
}


#' Return a readFit object base on nucleotide positions
position_fit.readFit <- function(df, h, pi=NULL)
{
  # create a template for each position
  npos <- ncol(h)
  if (npos==1)
    template_m <- matrix(T, nrow=1)
  else {
    template_m <- diag(npos)
    template_m <- apply(template_m, 2, as.logical)
  }

  out <- template_fit.readFit(template_m, df, h, pi=pi)

  names(out) <- colnames(h)

  return (out)
}

#' Build a readFit object base on reads
read_fit.readFit <- function(df, h, pi=NULL)
{
  # get all templates from read table
  template_m <- templates.read_table(df)
  out <- template_fit.readFit(template_m, df, h, pi=pi)

  return (out)
}

#' Build a readFit object base on templates.  This is a
#' generic function used to build readFit objects for
#' both position and read based fits
#' @export
template_fit.readFit <- function(template_m, df, h, pi=NULL)
{
  npos <- ncol(h)
  K <- nrow(h)
  ntemplates <- nrow(template_m)

  # for each template, generate all sampled_freqs, alleles,
  # P matrix, predicted_freqs
  info <- lapply(1:ntemplates, function(i) {
    ctemplate <- template_m[i,]
    allele_counts <- template_alleles.read_table(ctemplate, df, match=F)
    sampled_freqs <- allele_counts/sum(allele_counts)
    coverage <- sum(allele_counts)
    alleles <- names(sampled_freqs)
    nalleles <- length(alleles)

    # prepare nucs to be passed to haplotype_match
    nucs <- sapply(1:nalleles, function(i) {
      cnucs <- rep("+", npos)
      cnucs[ctemplate] <- strsplit(alleles[i], split="")[[1]]
      return (cnucs)
    })
    if (class(nucs) != "matrix")
      nucs <- matrix(nucs, nrow=nalleles)
    else
      nucs <- t(nucs)

    P <- haplotype_match.readFit(nucs, h)
    if (is.null(pi))
      predicted_freqs <- rep(NA, K)
    else
      predicted_freqs <- as.numeric(P %*% pi)

    return (list(sampled_freq=sampled_freqs,
                 predicted_freq=predicted_freqs,
                 coverage=coverage,
                 alleles=alleles,
                 P=P))
  })

  return (info)
}

#' Determine haplotypes that match a nucleotide pattern
#'
#' Given a nucleotide pattern covering some portion of
#' read table positions, return a vector of 1's and 0's
#' representing haplotypes that match the patterh
#'
#' @param nucs A character vector equal to the number of positions
#' in the read table.  Positions at which nuceotides are not specified
#' have value "+".
#' @param df read table
#' @param haplotype matrix
#'
#' @return A numeric vector of 1's and 0's
haplotype_match.readFit <- function(nucs, h)
{
  if (class(nucs) != "matrix")
    nucs <- matrix(nucs, nrow=1)

  if (ncol(nucs) != ncol(h))
    stop("nucs vector has wrong length")

  h_match <- apply(nucs, 1, function(cnucs) {
    active_pos <- which(cnucs != "+")
    x <- apply(h, 1, function(ch) all(ch[active_pos]==cnucs[active_pos]))
    as.numeric(x)
  })

  if (class(h_match) != "matrix")
    h_match <- matrix(h_match, ncol=nrow(h))
  else
    h_match <- t(h_match)

  return (h_match)
}

#' Visualize a readfit object
#'
#' @param rf a readFit object
#' @param outfile If not null, figure is written to outfile.
plot.readFit <- function(rf, outfile=NULL)
{
  ng <- length(rf)
  if (!is.null(names(rf)))
    labs <- names(rf)
  else
    labs <- paste("g", 1:ng, sep="")

  p_list <- lapply(1:ng, function(i) {
    crf <- rf[[i]]
    cp <- crf$predicted_freq
    cs <- crf$sampled_freq
    coverage <- crf$coverage
    alleles <- crf$alleles
    locus <- labs[i]
    df_p <- data.frame(freq=cp,
               coverage=coverage,
               type="predicted",
               alleles=alleles)
    df_s <- data.frame(freq=cs,
                       coverage=coverage,
                       type="sampled",
                       alleles=alleles)
    df <- rbind(df_p, df_s)

    p <- ggplot()
    p <- p + geom_bar(mapping=aes(x=alleles, y=freq, fill=type),
                      position="dodge", data=df,
                      stat="identity")
    if (i != 1)
      p <- p + theme(legend.position="none")
    p <- p + theme(axis.title.x = element_blank())
    p <- p + ggtitle(paste(locus, coverage, sep="-"))

    return (p)
  })

  m <- grid.arrange(grobs=p_list)

  if (!is.null(outfile)) {
    jpeg(outfile, width=2500, height = 1000)
    print(m)
    dev.off()
  }

  return (m)
}

