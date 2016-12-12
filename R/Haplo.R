#' Haplo object construct
#'
#' A Haplo object contains a collection of haplotypes,
#' their frequency, and total count
#'
#' @param h A haplotype matrix (matrix of characters with rows giving haplotypes)
#' @param pi A numeric vector giving frequencies of haplotypes
#' @param count The total coverage, may be NA if this is not well-defined
#'
#' @details A Haplo object is a list with elements h, pi, count, h_seq.
#' h is a character matrix with rows giving haplotypes, pi is a numeric
#' vector of frequencies, count is not currently used and h_seq is
#' a characer vector of haplotypes (rows of h pasted).
#'
#' @return A Haplo object
#' @export
Haplo <- function(h, pi, count=NA)
{
  if (class(h) != "matrix" | class(pi) != "numeric")
    stop("h or pi are not of right class")
  if (nrow(h) != length(pi))
    stop("h and pi do not imply the same number of haplotypes")

  h_seq <- aaply(h, 1, paste, collapse="")
  # if h is empty, make this a empty character vector
  if (length(h_seq)==0)
    h_seq <- as.character(c())

  ind <- order(pi, decreasing = T)

  h <- h[ind,,drop=F]
  pi <- pi[ind]
  h_seq <- h_seq[ind]

  H <- list(h=h, pi=pi, count=count, h_seq=h_seq)

  class(H) <- c("Haplo", "list")
  return (H)
}

#' @export
get_hap.Haplo <- function(H)
{
  return (H$h)
}

#' @export
get_nhap.Haplo <- function(H)
{
  return (nrow(get_hap.Haplo(H)))
}

#' @export
get_hap_seq.Haplo <- function(H)
{
  return (H$h_seq)
}

#' @export
get_freq.Haplo <- function(H)
{
  return (H$pi)
}

#' @export
get_count.Haplo <- function(H)
{
  return (H$count)
}

#' @export
get_nhap.Haplo <- function(H)
{
  return (length(H$pi))
}

#' Set the position names for the characters (nucleotides) forming the haplotypes.
#'
#' Alters the colnames of the haplotype matrix
#' to reflect pos names.   This is meant for visual purposes.   The Haplo object never accesses
#' these names.
#'
#' @param H A Haplo object
#' @param pos_names A character vector giving position names
set_pos_names.Haplo <- function(H, pos_names)
{
  h <- get_hap.Haplo(H)
  if (ncol(h) != length(pos_names))
    stop("h and pos_names imply different number of characters in haplotypes")

  colnames(H$h) <- pos_names

  return (H)
}

###########################################################################
plot.Haplo <- function(H, p=NULL, facet_label=NULL)
{
  pos <- rep(NA, get_nh.Haplo(H))
  HL <- HaploLocus(H, pos)

  return (plot.HaploLocus(HL, p, facet_label))
}

#' Returns Haplo object with identical haplotype frequencies merged.
unique.Haplo <- function(H)
{
  hseq <- get_hap_seq.Haplo(H)
  h <- get_hap.Haplo(H)
  pi <- get_freq.Haplo(H)

  hseq_u <- unique(hseq)
  inds <- lapply(hseq_u, function(hh) which(hh==hseq))

  pi_u <- sapply(inds, function(cinds) sum(pi[cinds]))
  h_u <- lapply(inds, function(cinds) h[cinds[1],,drop=F])
  h_u <- do.call(rbind, h_u)

  H_u <- Haplo(h_u, pi_u)
  return (H_u)
}

write.Haplo <- function(H, outfile)
{
  h <- get_hap.Haplo(H)
  f <- get_freq.Haplo(H)

  h_f <- cbind(f, h)

  write.table(h_f, outfile, row.names=F, col.names=F,
              sep=",")

  return (NULL)
}

read.Haplo <- function(infile)
{
  h_f <- read.table(infile, sep=",", colClasses="character")
  f <- as.numeric(h_f[,1])
  h <- as.matrix(h_f[,2:ncol(h_f)])

  return (Haplo(h, f))
}

#' Compare frequencies and haplotypes of H1 using
#' the frequencies and haplotype of H2 as a base
#'
#' Calculate the distance from the H1 haplotype to the
#' closest H2 haplotype.  Calculate the difference between
#' the frequencies of H1 and the frequencies of H2
compare.Haplo <- function(H1, H2)
{
  h1 <- get_hap.Haplo(H1)
  h2 <- get_hap.Haplo(H2)

  npos_h1 <- ncol(h1)
  npos_h2 <- ncol(h2)

  if (npos_h1 != npos_h2)
    stop("haplotypes have different number of positions")

  pi1 <- get_freq.Haplo(H1)
  pi2 <- get_freq.Haplo(H2)

  nh1 <- length(pi1)
  nh2 <- length(pi2)

  dist_m <- matrix(0, nrow=nh1, ncol=nh2)
  for (i in 1:nh1)
    for (j in 1:nh2)
      dist_m[i,j] <- sum(h1[i,] != h2[j,])

  dist_df <- adply(1:nh1, 1, function(hap_num) {
    ind <- which.min(dist_m[hap_num,])
    dist <- dist_m[hap_num, ind]
    data.frame(hap_num=hap_num,
               hap_match_num=ind,
               dist=dist)
  }, .id=NULL)



  nh_min <- min(nh1, nh2)
  pi_compare <- data.frame(hap_num=1:nh_min,
                           pi1=pi1[1:nh_min],
                           pi2=pi2[1:nh_min])

  return (list(dist=dist_df, pi_diff=pi_compare))
}
