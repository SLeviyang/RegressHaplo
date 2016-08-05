#' Construct for HaploLoci object
#'
#' A HaploLoci object is a collection of HaploLocus objects with utilities
#' for accessing the individual loci using their locus names
#'
#' @param HL_list A list of HaploLocus objects
#'
#' @return A HaploLoci object
#' @export
HaploLoci <- function(HL_list)
{
  if (class(HL_list) != "list")
    stop("HL_list must be a list")

  if (length(HL_list)==0)
    return (NULL)

  class_vec <- sapply(HL_list, function(HL) class(HL)[1])
  if (!all(class_vec=="HaploLocus"))
    stop("HL_list must contain only HaploLocus objects")

  locus_names <- sapply(HL_list, get_name.HaploLocus)

  #check that names are unique
  not_na_ind <- !is.na(locus_names)
  if (length(unique(locus_names[not_na_ind])) != sum(not_na_ind))
    stop("HL_list contains loci with the same name.  All names must be unique, if not NA")

  names(HL_list) <- locus_names

  class(HL_list) <- c("HaploLoci", "list")
  return (HL_list)
}

get_locus_names.HaploLoci <- function(HL_list)
{
  return (names(HL_list))
}

subset.HaploLoci <- function(HL_list, ind=NULL, locus_names=NULL)
{
  HL_list <- get_HaploLocus.HaploLoci(HL_list, ind, locus_names)
  class(HL_list) <- c("HaploLoci", "list")

  return (HL_list)
}

get_HaploLocus.HaploLoci <- function(HL_list, ind=NULL, locus_names=NULL)
{
  # if both ind and locus_names are NULL return a list of the HaploLocus objects
  if (is.null(ind) & is.null(locus_names))
    return (HL_list)

  if (!is.null(ind)) {
    return (HL_list[ind])
  }

  if (class(locus_names) != "character")
    stop("locus_names must be a character vector")

  class(HL_list) <- c("list")

  return (HL_list[locus_names])
}
#############################################################
#  methods to manipulate/analyze HaploLoci object
loci_df.HaploLoci <- function(HL)
{
  HL_list <- get_HaploLocus.HaploLoci(HL)
  locus_names <- get_locus_names.HaploLoci(HL)
  names(HL_list) <- locus_names

  df <- ldply(HL_list, function(Hlocus) {
    pos <- get_pos.HaploLocus(Hlocus)
    data.frame(start=min(pos), end=max(pos))
  }, .id="locus")

  return (df)
}
nloci.HaploLoci <- function(HL_list)
{
  return (length(get_locus_names.HaploLoci(HL)))
}

# Get positions over all loci
pos.HaploLoci <- function(HLoci)
{
  HLocus_list <- get_HaploLocus.HaploLoci(HLoci)
  pos_list <- lapply(HLocus_list, get_pos.HaploLocus)
  pos_list <- unlist(pos_list)
  names(pos_list) <- NULL

  return (pos_list)
}

# Get haplotypes across all loci
haplotypes.HaploLoci <- function(HLoci)
{
  HLocus_list <- get_HaploLocus.HaploLoci(HLoci)
  H_list <- lapply(HLocus_list, get_Haplo.HaploLocus)
  h_list <- lapply(H_list, get_hap.Haplo)

  locus_names <- get_locus_names.HaploLoci(HLoci)
  names(h_list) <- locus_names

  return (h_list)
}

nhaplotypes.HaploLoci <- function(HLoci)
{
  h_list <- haplotypes.HaploLoci(HLoci)
  nh <- sapply(h_list, nrow)

  return (nh)
}

# get global haplotypes across all loci by permuting locus haplotypes
global_haplotypes.HaploLoci <- function(HLoci)
{
  h_list <- haplotypes.HaploLoci(HLoci)
  all_pos <- pos.HaploLoci(HLoci)

  nh <- sapply(h_list, nrow)
  if (prod(nh) > 10^4) {
    cat("number of haplotypes", prod(nh), "\n")
    return (NA)
  }

  if (length(h_list)==1)
    return (h_list[[1]])

  con_haps_seqs <- lapply(h_list, function(h)
    apply(h, 1, paste, collapse=""))
  merged_con_haps_seqs <- as.matrix(expand.grid(con_haps_seqs))
  merged_haps <- aaply(merged_con_haps_seqs, 1, function(seqs) {
    allseqs <- paste(seqs, collapse="")
    split_seqs <- strsplit(allseqs, split="")[[1]]
  }, .drop=F)

  colnames(merged_haps) <- all_pos

  return (merged_haps)
}

#' Plot a HaploLoci object
#'
#' @param HL A HaploLoci object
#' @param p A ggplot object.  If null then initialize a ggplot object, otherwise
#' just add a layer
#' @param facet_label Either NULL or a string.  If not NULL, then the string is
#' add to all layer data.frame's under the column facet_label.  Used to
#' facet multiple HaploLoci objects
plot.HaploLoci <- function(HL, p=NULL, facet_label=NULL)
{
  HLocus_list <- get_HaploLocus.HaploLoci(HL)
  all_pos <- pos.HaploLoci(HL)
  npos <- length(all_pos)

  plot_df <- ldply(HLocus_list, function(HLocus) {
    pos <- get_pos.HaploLocus(HLocus)
    H <- get_Haplo.HaploLocus(HLocus)

    h <- get_hap.Haplo(H)
    h_seq <- get_hap_seq.Haplo(H)
    freq <- get_freq.Haplo(H)

    pos_ind <- match(pos, all_pos)
    ymax <- cumsum(freq)
    ymin <- c(0, ymax[-length(ymax)])

    full_df <- data.frame(posmin=min(pos), posmax=max(pos),
                 allele=h_seq, freq=freq,
                 xmin=min(pos_ind), xmax=max(pos_ind)+1,
                 ymin=ymin, ymax=ymax)

#     npos <- length(pos)
#     locus_df <- lapply(1:npos, function(i) {
#       nucs <- factor(h[,i], levels = c("A", "C", "G", "T", "-"))
#       data.frame(pos=pos[i], allele=seqs, freq=freq,
#                  xmin=pos_ind[i], xmax=pos_ind[i]+1,
#                  ymin=ymin, ymax=ymax)
#     })
#     full_df <- do.call(rbind, locus_df)
    return (full_df)
  }, .id=NULL)

  loci_df <- loci_df.HaploLoci(HL)
  loci_plot_df <- ddply(loci_df, .(locus), function(locus_df) {
    start_ind <- which(all_pos == locus_df$start[1])
    end_ind <- which(all_pos == locus_df$end[1])+1

    df1 <- mutate(locus_df, x=start_ind, xend=start_ind,
                  y=0, yend=1)
    df2 <- mutate(locus_df, x=end_ind, xend=end_ind,
                  y=0, yend=1)
    rbind(df1, df2)
  })

  if (!is.null(facet_label)) {
    plot_df <- mutate(plot_df, facet_label=facet_label)
    loci_df <- mutate(loci_df, facet_label=facet_label)
  }

  new_plot <- is.null(p)
  if (new_plot)
    p <- ggplot()


  p <- p + geom_rect(mapping=aes(xmin=xmin, xmax=xmax,
                                 ymin=ymin, ymax=ymax,
                                 fill=allele),
                     data=plot_df, col="black", size=.5)

  p <- p + geom_text(mapping=aes(x=(xmin+xmax)/2, y=(ymax+ymin)/2,
                                   label=allele),
                       data=plot_df)

  p <- p + geom_segment(mapping=aes(x=x, y=y, xend=xend, yend=yend),
                        data=loci_plot_df, col="yellow", size=4)

  if (new_plot) {
    p <- p + scale_y_continuous(expand = c(0,0),
                                breaks=NULL,
                                limits=c(0, 1))

    p <- p + scale_x_continuous(expand = c(0,0),
                                breaks=(1:npos)+.5,
                                labels=all_pos)


    #p <- p + scale_fill_manual(values=c("red", "blue", "brown", "green",
    #                                    "grey"), drop = FALSE)
    p <- p + theme(legend.position="none")

    p <- p + theme(axis.text=element_text(size=12))
    p <- p + xlab("position") + ylab("frequency")
  }

  return (p)
}

#############################################################
## methods connecting HaploLoci object to other objects

#' Restrict a read table to positions covered by loci
#'
#' @param HLoci a HaploLoci object
#' @param df a read table
#' @param min_count minimum value for keeping a read
#'
#' @return A read table
filter_read_table.HaploLoci <- function(HLoci, df, min_count=0)
{
  all_pos <- pos.HaploLoci(HLoci)
  df <- subset.read_table(df, pos=all_pos)
  df <- clean.read_table(df, min_count=min_count)

  return (df)
}

#' Apply RegressHaplo to haplotypes given by a HaploLoci object and
#' data given by a read table
apply_RegressHaplo.HaploLoci <- function(HLoci, df, rho, plot=T,
                                         join_loci=F)
{
  df_filtered <- filter_read_table.HaploLoci(HLoci, df, min_count=20)
  h <- global_haplotypes.HaploLoci(HLoci)

  rh <- compute_solution.RegressHaplo(df_filtered, h, rho)
  h_out <- get_h.RegressHaploLocal(rh)
  pi_out <- get_pi.RegressHaploLocal(rh)
  pos <- as.numeric(colnames(h_out))

  H <- Haplo(h_out, pi_out)
  H_locus <- HaploLocus(H, pos)

  if (plot) {
    HLocus_list <- get_HaploLocus.HaploLoci(HLoci)
    pos_list <- lapply(HLocus_list, get_pos.HaploLocus)
    locus_names <- get_locus_names.HaploLoci(HLoci)

    H_locus_split <- split.HaploLocus(H_locus, pos_list,
                                      locus_names)

    p <- plot(HLoci, facet_label="local")
    p <- plot(H_locus_split, p=p, facet_label="global")
    p <- p + facet_grid(facet_label ~ .)
    print(p)
  }

  if (join_loci)
    return (H_locus)
  else {
    HLocus_list <- get_HaploLocus.HaploLoci(HLoci)
    pos_list <- lapply(HLocus_list, get_pos.HaploLocus)
    locus_names <- get_locus_names.HaploLoci(HLoci)

    H_locus_split <- split.HaploLocus(H_locus, pos_list,
                                      locus_names)
    return (H_locus_split)
  }

  return (H_locus)
}



