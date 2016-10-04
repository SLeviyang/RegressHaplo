#' Construct a HaploLocus object
#'
#' A HaploLocus object combines a Haplo object with positions and a locus name
#'
#' @param H Haplo object
#' @param pos positions corresponding to each column of the haplotype matrix in the Haplo object
#' @param locus_name A string
#'
#' @return A HaploLocus object
#' @export
HaploLocus <- function(H, pos, locus_name=NA)
{
 # if (class(H)[1] != "Haplo" | class(pos) != "numeric")
#    stop("H or pos of wrong class")

  h <- get_hap.Haplo(H)
  if (ncol(h) != length(pos))
    stop("H and pos do not imply the same number of haplotype characeters (nucleotides)")

  locus_names <- as.character(locus_name)
  H <- set_pos_names.Haplo(H, as.character(pos))
  HL <- list(H=H, pos=pos, locus_name=locus_name)

  class(HL) <- c("HaploLocus", "list")
  return (HL)
}

get_Haplo.HaploLocus <- function(HL)
{
  return (HL$H)
}

get_pos.HaploLocus <- function(HL)
{
  return (HL$pos)
}

get_name.HaploLocus <- function(HL)
{
  return (HL$locus_name)
}

#######################################################################
plot.HaploLocus <- function(HL, p=NULL, facet_label=NULL)
{
  HLoci <- HaploLoci(list(HL))
  return (plot.HaploLoci(HLoci, p, facet_label))
}

#' Split a HaploLocus object into multiple loci to form a HaploLoci object
split.HaploLocus <- function(HL, pos_list, locus_names=NA)
{
  pos <- get_pos.HaploLocus(HL)

  # check that each vector of positions is non-empty
  pos_list_length <- sapply(pos_list, length)
  if (!all(pos_list_length > 0))
    stop("each entry of pos_list must be of positive length")

  # check the each vector of positions is a subset of Haplo pos
  sub <- sapply(pos_list, function(cpos) length(setdiff(cpos, pos))==0)
  if (!all(sub))
    stop("each entry of pos_list must be a subset of the haplotype positions")

  H <- get_Haplo.HaploLocus(HL)
  freq <- get_freq.Haplo(H)
  h <- get_hap.Haplo(H)

  HL_list <- Map(function(cpos, clocus_name) {
    pos_match <- match(pos, cpos)
    active_ind <- which(!is.na(pos_match))

    new_H <- Haplo(h[,active_ind,drop=F], freq)
    new_HL <- HaploLocus(new_H, cpos, clocus_name)

    return (new_HL)
  }, pos_list, locus_names)

  return (HaploLoci(HL_list))
}
