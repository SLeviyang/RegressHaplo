#' Read table constructor
#'
#' Create a read table from a BAM file and variant calls
#'
#' @param bam_file bam file
#' @param bai_file index file.  If missing then .bai appendix on bam_file is assumed.
#' @param variant_calls.  A data.frame specifying positions and nucleotides on
#' reference to be taken as true variants, see details.
#' @param pu A BAM pileup object returned by BAM_pileup method.  If NULL, a pile up
#' will be created
#'
#' @return
#' A data.frame with columns:  count pos1 pos2 .. posn.
#' Each row of the data.frame
#' corresponds to a collection of reads that
#' are identical at the variable positions.
#' The count column gives the number of reads
#' with the values specified in the posi columns.  Every read is listed in the table, so
#' that the sum of the count column gives the total number of reads.
#' The posi are the positions specified by variant_calls, i.e. posi==as.character(variant_calls$pos[i]).
#' Entries in the posi columns are either A/C/G/T/-/NA.  NA means the read did not cover that variable
#' position.
#'
#' @details variant_calls has 6 columns which must be names as c("pos", "A", "C", "G", "T", "-").
#' pos gives the position at which a true variant exists.  The other columns have logical entries,
#' with TRUE meaning that a true variant exists with the corresponding nucleotide.
#'
#' @export
read_table <- function(bam_file,
                       bai_file=NULL,
                       variant_calls,
                       pu=NULL,
                       debug=F)
{
  if (is.null(bai_file)) {
    bai_file <- paste(bam_file, ".bai", sep="")
  }
  cat("parsing", bam_file, "\n")

  # on first pass retrieve all variants at called positions
  # then filter for correct variants, see bottom of function
  nt_pos <- variant_calls$pos
  start_pos <- min(nt_pos)
  end_pos <- max(nt_pos)

  bam_header <- scanBamHeader(bam_file)[[1]]
  ref_name <- names(bam_header$targets)

  if (debug) {
    print("BAM header and reference")
    print(bam_header)
    print(ref_name)
    print("Variable positions")
    print(nt_pos)
  }

  which <- GRanges(seqnames = ref_name,
                   ranges=IRanges(start_pos, end_pos))
  what <- c("seq", "pos", "qname")
  param <- ScanBamParam(which=which, what=what)

  if (debug) {
    print("BAM input parameters")
    print(param)
  }


  # this throws warnings if pairings are not matched
  # let's not freak out users, so turn off warnings
  oldw <- getOption("warn")
  options(warn = -1)

  if (debug) {
    print("Preparing to read alignments.  GenomicAlignments ver:")
    print(package.version("GenomicAlignments"))
  }
  ga_pair <- readGAlignmentPairs(bam_file, bai_file,
                                 param=param)

  if (debug) {
    print("Full Alignment")
    print(ga_pair)
  }

  options(warn = oldw)


  # check if this file has paired-end reads
  if (length(ga_pair)==0) {
    print("PROCESSING FILE AS SINGLE END READS")
    paired_end <- F
    ga <- readGAlignments(bam_file, bai_file, param=param)
    read_strings <- single_end_read_table(ga, nt_pos)
  } else {
    print("PROCESSING FILE AS PAIRED END READS")
    paired_end <- T
    read_strings <- paired_end_read_table(ga_pair, nt_pos,
                                          debug=debug)
  }

  read_strings_t <- table(read_strings)

  # put the reads in order of pos
  read_first_pos <- sapply(names(read_strings_t), function(s) {
    # if read doesn't cover variable positions put it last
    if (s=="outside")
      return (10^6)

    as.numeric(strsplit(s, split=":")[[1]][1])
  })
  read_strings_t <- read_strings_t[order(read_first_pos)]

  template <- matrix(as.character(NA), nrow=1, ncol=length(nt_pos))
  colnames(template) <- as.character(nt_pos)

  # construct rows of read table
  ig <- 1
  df <- adply(names(read_strings_t), 1, function(s) {

    m <- template
    if (s=="outside") {
      return (data.frame(m, stringsAsFactors = F))
    }

    ssplit <- strsplit(s, split="/")[[1]]
    spos <- sapply(ssplit, function(ss) strsplit(ss, split=":")[[1]][1])
    snuc <- sapply(ssplit, function(ss) strsplit(ss, split=":")[[1]][2])

    m[,spos] <- snuc

    return (data.frame(m, stringsAsFactors = F))
  }, .expand=T, .id=NULL)

  # add counts to front of data.frame (how to use mutate to do this?)
  df <- cbind(as.numeric(read_strings_t), df)
  names(df) <- c("count", as.character(nt_pos))

  class(df) <- c("read_table", "data.frame")

  # now remove variants which are not called
  cat("filtering for variant calls", "\n")
  if (is.null(pu))
    pu <- BAM_pileup(bam_file, max_depth=5000, min_base_quality=0, min_mapq=0)
  consensus_values <- consensus(pu)[variant_calls$pos]
  df <- filter_true_variants.read_table(df, variant_calls, consensus_values)

  return (df)
}


  ########################## non-paired end
#' A helper function that creates read tables from single end reads
#'
#' Should not be called directly by user, but instead is accessed
#' by calling read_table
#'
single_end_read_table <- function(ga, nt_pos)
{
  seq <- as.character(mcols(ga)$seq)
  seq_on_ref <- create_refspace_seq(seq, cigar(ga))
  seq_on_ref_split <- strsplit(seq_on_ref, split="")

  nreads <- length(seq)

  read_strings <- sapply(1:nreads, function(i) {
    if (i %% 1000 == 0)
      cat("processing read", i, "of", nreads, "\n")

    pos <- mcols(ga)$pos[i]#  bam_info$pos[i]
    cseq <- seq_on_ref_split[[i]]
    all_pos <- pos:(pos+length(cseq)-1)

    variant_pos <- is.element(all_pos, nt_pos)

    # if the read doesn't cover variable pos, return "outside"
    if (!any(variant_pos))
      return ("outside")

    paste(all_pos[variant_pos], cseq[variant_pos], sep=":", collapse="/")
  })

  return (read_strings)
}

#' A helper function that creates read tables from paired end reads
#'
#' Should not be called directly by user, but instead is accessed
#' by calling read_table
#'
paired_end_read_table <- function(ga_pair,
                                  nt_pos,
                                  debug=F)
{
  ga1 <- GenomicAlignments::first(ga_pair)
  ga2 <- GenomicAlignments::last(ga_pair)

  seq1 <- as.character(mcols(ga1)$seq)
  seq2 <- as.character(mcols(ga2)$seq)

  if (debug) {
    print("Full pairing information")
    print(ga_pair)
    print("First end information")
    print(ga1)
    print("Second end information")
    print(ga2)

    nseq1 <- min(length(seq1), 2)
    nseq2 <- min(length(seq2), 2)

    print("First end sequences")
    if (nseq1 > 0)
      print(seq1[1:nseq1])
    print("Second end sequences")
    if (nseq2 > 0)
      print(seq2[1:nseq2])
  }

  seq_on_ref1 <- create_refspace_seq(seq1, cigar(ga1))
  seq_on_ref2 <- create_refspace_seq(seq2, cigar(ga2))

  seq_on_ref_split1 <- strsplit(seq_on_ref1, split="")
  seq_on_ref_split2 <- strsplit(seq_on_ref2, split="")

  nreads <- length(seq1)

  read_strings <- sapply(1:nreads, function(i) {
    if (i %% 1000 == 0)
      cat("processing read", i, "of", nreads, "\n")

    f_pos1 <- mcols(ga1)$pos[i]
    f_pos2 <- mcols(ga2)$pos[i]

    # make sure the first seq is downstream
    if (f_pos1 < f_pos2) {
      cseq1 <- seq_on_ref_split1[[i]]
      cseq2 <- seq_on_ref_split2[[i]]
      pos1 <- mcols(ga1)$pos[i]
      pos2 <- mcols(ga2)$pos[i]
    } else {
      cseq1 <- seq_on_ref_split2[[i]]
      cseq2 <- seq_on_ref_split1[[i]]
      pos1 <- mcols(ga2)$pos[i]
      pos2 <- mcols(ga1)$pos[i]
    }

    all_pos1 <- pos1:(pos1+length(cseq1)-1)
    all_pos2 <- pos2:(pos2+length(cseq2)-1)

    # if there is overlap, clip second sequence
    overlap_all_pos <- base::intersect(all_pos1, all_pos2)
    if (length(overlap_all_pos) > 0) {
      overlap_ind <- is.element(all_pos2, overlap_all_pos)
      cseq2 <- cseq2[!overlap_ind]
      all_pos2 <- all_pos2[!overlap_ind]
    }

    cseq <- c(cseq1, cseq2)
    all_pos <- c(all_pos1, all_pos2)

    variant_pos <- is.element(all_pos, nt_pos)

    # if the read doesn't cover variable pos, return "outside"
    if (!any(variant_pos))
      return ("outside")

    sout <- paste(all_pos[variant_pos],
                  cseq[variant_pos], sep=":", collapse="/")

    return (sout)
  })

  return (read_strings)
}

########################################################
#  methods that analyze read_table objects

#' Returns positions of read table as character vector
#'
#' @return A character vector
pos_names.read_table <- function(df)
{
  return (names(dplyr::select(df, -count)))
}

#' Return the consensus sequence at each position of the read table
consensus.read_table <- function(df)
{
  df_nucs <- dplyr::select(df, -count)
  nucs <- apply(df_nucs, 2, function(cc) {
    nuc_t <- table(cc)
    if (length(nuc_t)==0)
      return (NA)
    max_ind <- which.max(nuc_t)
    return (names(nuc_t)[max_ind])
  })

  names(nucs) <- names(df_nucs)
  return (nucs)
}

#' Return coverage at each position in read table
#'
#' @return A numeric vector
coverage.read_table <- function(df)
{
  counts <- df$count
  df_nucs <- dplyr::select(df, -count)
  coverage <- apply(df_nucs, 2, function(cc) {
    active <- which(!is.na(cc))
    sum(counts[active])
  })

  names(coverage) <- names(df_nucs)

  return (coverage)
}


#' Return all alleles matching a read template
#'
#' @param template A logical vector
#' @param df a read table
#' @param match if T then only return alleles that comprise
#' an entire read.  if F then return all alleles that comprise
#' any subsequence of a read.
#'
#' @details a template is a logical vector of length
#' equal to the number of positions in the read table
#' and with T corresponding to nucleotide positions that are
#' part of the template.
#'
#' @return A named numeric vector with entries giving allele counts
#' and names giving alleles
template_alleles.read_table <- function(template, df, match=F)
{
  df_nucs <- dplyr::select(df, -count)
  counts <- df$count
  if (length(template) != ncol(df_nucs))
    stop("template length does not match number of positions in read table")

  active_pos <- which(template)

  # remove all reads that do not have nucleotides at all active positions
  ind <- apply(df_nucs, 1, function(cread) {
    covered_pos <- which(!is.na(cread))
    length(setdiff(active_pos, covered_pos)) == 0
  })

  df_nucs <- df_nucs[ind,,drop=F]
  counts <- counts[ind]
  if (nrow(df_nucs)==0)
    return (NULL)

  # if we want match, remove all reads that have nucleotides outside active
  # positions
  if (match) {
    ind <- apply(df_nucs, 1, function(cread) {
      covered_pos <- which(!is.na(cread))
      length(setdiff(covered_pos, active_pos)) == 0
    })
    df_nucs <- df_nucs[ind,,drop=F]
    counts <- counts[ind]
  }

  if (nrow(df_nucs)==0)
    return (NULL)

  alleles <- apply(df_nucs, 1, function(cread) paste(cread[active_pos], collapse=""))
  unique_alleles <- unique(alleles)

  unique_allele_counts <- sapply(unique_alleles, function(ca)
    sum(counts[which(alleles==ca)]))

  names(unique_allele_counts) <- unique_alleles

  return (unique_allele_counts)
}

#' Return read indices corresponding to the template
#'
#' @param template A logical vector or matrix. If matrix, rows
#' contain templates.  Templates must be unique
#' @param df a read table
#'
#' @return A list containing indices in df or read groups (rows)
#' that have the given templates
template_indices.read_table <- function(template, df)
{
  if (!is.matrix(template))
    template <- matrix(template, nrow=1)

  template_s <- apply(template, 1, function(s)
    paste(as.numeric(s), collapse=""))

  if (length(template_s) != length(unique(template_s)))
    stop("templates must be unique")

  df_reads <- as.matrix(dplyr::select(df, -count))
  df_s <- apply(df_reads, 1, function(s) {
    z <- as.numeric(!is.na(s))
    paste(z, collapse="")
  })

  ind <- lapply(template_s, function(s) which(df_s==s))

  return (ind)
}

#' Retrieve all templates matching at least one read in the
#' read table
#'
#' Each read determines a template, with T at the positions
#' covered by the read, but multiple reads can have the same
#' template.  This function returns all unique templates.
#'
#' @param df read table
#'
#' @return a logical matrix with each row giving a template
templates.read_table <- function(df)
{
  df_nucs <- dplyr::select(df, -count)
  template_m <- matrix(!is.na(df_nucs),
                       nrow=nrow(df_nucs),
                       ncol=ncol(df_nucs))

  return (unique(template_m))
}

#' For each position, return reads that cover the position

#' @return A R by P logical matrix where R is number of reads and P is
#' number of positions.  A TRUE entry means the read covers the position
reads_covering_positions.read_table <- function(df)
{
  df_nucs <- dplyr::select(df, -count)

  m <- sapply(1:ncol(df_nucs), function(i) {
    ind <- !is.na(df_nucs[,i])
    return (ind)
  })
  if (class(m) != "matrix")
    m <- matrix(m, ncol=ncol(df_nucs))

  colnames(m) <- names(df_nucs)

  return (m)
}

#' For each haplotype, return reads that match haplotype at covered positions
#'
#' @param df read table
#' @param h haplotype matrix
#'
#' @return A R by K logical matrix where R is number of reads and K is
#' number of haplotypes.  A TRUE entry means the read covers the haplotype.
reads_covering_haplotypes.read_table <- function(df, h)
{
  df_nucs <- dplyr::select(df, -count)
  pos_m <- reads_covering_positions.read_table(df)
  if (ncol(pos_m) != ncol(h))
    stop("haplotype and read table do not have same number of positions")

  nread <- nrow(df)
  m <- sapply(1:nread, function(i) {
    covered_pos <- which(pos_m[i,])
    cread <- df_nucs[i,covered_pos]
    apply(h, 1, function(ch) all(ch[covered_pos]==cread))
  })
  m <- t(m)

  return (m)
}

#' Return nucleotide counts at each position
#'
#' @param df read_table
#'
#' @return A matrix with 5 rows and a column for each position in df
nucs_at_pos.read_table <- function(df)
{
  df_nucs <- dplyr::select(df, -count)

  M <- sapply(1:ncol(df_nucs), function(i) {
    all_nuc_vec <- rep(0, 5)
    names(all_nuc_vec) <- c("A", "C", "G", "T", "-")

    ind <- which(!is.na(df_nucs[,i]))
    if (length(ind)==0)
      return (all_nuc_vec)

    d <- data.frame(count=df$count[ind],
                    nuc=df_nucs[ind,i], stringsAsFactors = F)
    total_count <- sum(d$count)
    nuc_df <- ddply(d, .(nuc), function(nuc_df)
              data.frame(freq=sum(nuc_df$count)/total_count,
                         nuc=nuc_df$nuc[1], stringsAsFactors = F))

    all_nuc_vec[nuc_df$nu] <- nuc_df$freq

    return (all_nuc_vec)
  })

  rownames(M) <- c("A", "C", "G", "T", "-")
  colnames(M) <- names(df_nucs)

  return (M)
}

#' Return the start position of each read group in a read table
#'
#' @param df read_table object
#'
#' @return a numeric vector of start positions
start_pos.read_table <- function(df)
{
  all_pos <- all_pos.read_table(df)

  return(sapply(all_pos, min))
}

#' Return the end position of each read group in a read table
#'
#' @param df read_table object
#'
#' @return a numeric vector of start positions
end_pos.read_table <- function(df)
{
  all_pos <- all_pos.read_table(df)

  return(sapply(all_pos, max))
}

#' Return the positions covered by each read group in a read table
#'
#' @param df read_table object
#'
#' @return a list of numeric vectors containing read positions
all_pos.read_table <- function(df)
{
  df_nucs <- dplyr::select(df, -count)
  nrg <- nrow(df_nucs)
  pos <- as.numeric(names(df_nucs))

  all_pos <- lapply(1:nrg, function(i) {
    ind <- which(!is.na(df_nucs[i,]))
    pos[ind]
  })

  return (all_pos)
}

#' Return the nucleotides in each read group of a read table
#'
#' @param df read_table object
#'
#' @return a list of character vectors containing read nucleotides
seq.read_table <- function(df)
{
  df_nucs <- dplyr::select(df, -count)
  nrg <- nrow(df_nucs)

  seq <- lapply(1:nrg, function(i) {
    ind <- which(!is.na(df_nucs[i,]))
    as.character(df_nucs[i,ind])
  })

  return (seq)
}

#' Create edge matrix for read graph
#'
#' @param df read_table object
#' @param sparse If TRUE, remove redundant edges
#'
#' @return a square matrix with dimension equal to the number
#' of read groups in the read graph
adjacency_matrix.read_table <- function(df, sparse=T)
{
  # if there are too many rows, computation should not be done because
  # it will be slow and cause a memory problem
  if (nrow(df) > 200)
    return (NULL)

  dfs <- split_paired_ends.read_table(df)
  if (nrow(dfs$pairs) != 0) {
    s <- paste("adjacency matrix cannot be constructed for paired reads with gap.",
                "Call split_paired_ends.read_table", sep=" ")
    stop(s)
  }
  start_pos <- start_pos.read_table(df)
  end_pos <- end_pos.read_table(df)
  pos <- all_pos.read_table(df)
  all_pos <- setdiff(names(df), "count")
  seq <- seq.read_table(df)
  nseq <- length(seq)
  nreads <- nrow(df)

  sp_order <- order(start_pos, end_pos)

  m <- matrix(0, nrow=nseq, ncol=nseq)

  # if there is only 1 read, then m=0 and there is
  # no need to check for sparsity
  if (nseq==1)
    return (m)

  for (i_sp in 1:(nseq-1))
    for (j_sp in (i_sp+1):nseq) {
      # use start pos order, so we know
      # read i start pos <= read j start pos

      i <- sp_order[i_sp]; j <- sp_order[j_sp];
      shared_pos <- intersect(pos[[i]], pos[[j]])

      if (length(shared_pos)==0)
        next

      i_shared_ind <- which(is.element(pos[[i]], shared_pos))
      j_shared_ind <- which(is.element(pos[[j]], shared_pos))

      # check that overlap nucleotides match
      if(!all(seq[[i]][i_shared_ind] == seq[[j]][j_shared_ind]))
        next

      # if the ith read doesn't go beyond the jth read
      i_end_pos <- max(pos[[i]])
      j_end_pos <- max(pos[[j]])
      if (i_end_pos <= j_end_pos)
        m[i,j] <- 1

    }

  # if there are no edges, then just return, no need to check for sparse
  if (all(m==0))
    return (m)

  npos <- ncol(df) - 1
  nm_list <- max(2, nreads-1)
  m_list <- list(rep(NA, nm_list))
  m_list[[1]] <- m %*% m
  for (i in 2:nm_list) {
    m_list[[i]] <- m_list[[i-1]] %*% m
  }
  m_pow <- sapply(m_list, as.numeric)

  # a redundant edge exists if m[i,j]=1 and m^k[i,j]=1 for k > 1
  if (sparse) {
    m_vec <- as.numeric(m)
    redundant <- sapply(1:length(m_vec), function(i)
             m_vec[i]==1 & any(m_pow[i,]>=1))
    m_sparse <- ifelse(redundant, 0, m_vec)
    m_sparse <- matrix(m_sparse, nrow=nrow(m), ncol=ncol(m))
    return (m_sparse)
  }
  else
    return (m)
}

#' Determines loci in read table that have no overlapping reads
#'
#' @param df a read table object
#'
#' @return A named logical vector with an entry for every column, other than
#' count, in the read table data.frame.  Names of entries are the corresponding
#' column names (positions) in the read tabel data.frame.
#'
#' @details If an entry in the returned logical vector is true then all
#' positions prior to the entry are unlinked to all positions corresponding
#' to the entry and greater.
unlinked_pos.read_table <- function(df, min_cover)
{
  df_nucs <- dplyr::select(df, -count)
  npos <- ncol(df_nucs)
  pos_ch <- pos_names.read_table(df)
  pos <- as.numeric(pos_ch)

  if (npos==1) {
    unlinked <- F
    names(unlinked) <- pos_ch
    return (unlinked)
  }

  #pos <- as.numeric(names(df_nucs))


  linkage_break <- sapply(2:npos, function(i) {
    # go through read groups and check for read that spans <i, >=i
    spans <- apply(df_nucs, 1, function(rr) {
      any(!is.na(rr[1:(i-1)])) & any(!is.na(rr[i:npos]))
    })

    # we have a linkage break if there are no reads that span
    if (all(!spans))
      return (T)
    else
      return (sum(df$count[spans]) < min_cover)
  })

  unlinked <- c(F, linkage_break)
  #names(unlinked) <- setdiff(names(df), "count")
  names(unlinked) <- pos_ch

  return (unlinked)
}

#' Split reads which have a gap between the paired ends
#'
#' @details reads
#' with no gap are not split, even though they might be from
#' different ends
#'
#' @return a data.frame with split reads and a matrix containing
#' rows that are paired in the returned data.frame
split_paired_ends.read_table <- function(df)
{
  df_nucs <- dplyr::select(df, -count)
  # if not paired_read, return 0, otherwise return index in df_nucs
  #  that splits the two reads

  paired_reads <- apply(df_nucs, 1, function(cread) {
    active_ind <- which(!is.na(cread))
    if (length(active_ind)==0)
      return (0)
    split_ind <- setdiff(min(active_ind):max(active_ind), active_ind)

    if (length(split_ind)==0)
      return (0)
    else
      return (min(split_ind))
  })

  df_np <- filter(df, paired_reads==0)
  df_p <- filter(df, paired_reads>0)
  split_ind <- paired_reads[paired_reads > 0]

  if (length(split_ind)==0)
    return (list(df=df_np, pairs=matrix(nrow=0, ncol=2)))

  # adjust for count column
  split_ind <- split_ind + 1

  df_p_split_list <- Map(function(row_ind, i) {
    cread <- df_p[row_ind,]


    first_read <- cread
    first_read[i:ncol(cread)] <- NA
    second_read <- cread
    second_read[2:i] <- NA

    return (rbind(first_read, second_read))
  }, 1:nrow(df_p), split_ind)
  df_p_split <- do.call(rbind, df_p_split_list)

  m_split <- matrix(1:nrow(df_p_split), byrow=T,
                    nrow=nrow(df_p_split)/2, ncol=2)
  colnames(m_split) <- c("first", "last")

  df_out <- rbind(df_p_split, df_np)

  return (list(df=df_out, pairs=m_split))
}

#' A helper function that create haplotype permutations from a list
#' of haplotype matrices
#'
#' @param con_haps A list of haplotype matrices
merge_haplotypes.read_table <- function(con_haps)
{
  if (length(con_haps)==1)
    return (con_haps[[1]])

  con_haps_seqs <- lapply(con_haps, function(h)
    apply(h, 1, paste, collapse=""))
  merged_con_haps_seqs <- expand.grid(con_haps_seqs)
  merged_haps <- t(apply(merged_con_haps_seqs, 1, function(seqs) {
    unlist(lapply(seqs, function(cs) strsplit(cs, split="")[[1]]))
  }))
  merged_haps_pos <- unlist(lapply(con_haps, colnames))

  colnames(merged_haps) <- merged_haps_pos
  ind <- order(as.numeric(merged_haps_pos))
  merged_haps <- merged_haps[,ind,drop=F]

  return (merged_haps)
}

#' A helper function that converts read graph paths to haplotype matrices
#'
#' @param paths a list of paths generated by igraph function all_simple_paths
#' @param df the read table from which paths was built
#'
#' @return A haplotype matrix
paths_to_haplotypes.read_table <- function(paths, df)
{
  seq <- seq.read_table(df)
  pos <- all_pos.read_table(df)

  all_pos <- as.numeric(pos_names.read_table(df))
  pos_ind <- lapply(pos, function(cpos) match(cpos, all_pos))

  npaths <- length(paths)
  npos <- length(all_pos)

  m <- matrix(NA, nrow=npaths, ncol=npos)
  colnames(m) <- as.character(all_pos)

  for (i in 1:npaths) {
    cpath <- paths[[i]]
    seq_vals <- unlist(seq[cpath])
    pos_ind_vals <- unlist(pos_ind[cpath])

    m[i,pos_ind_vals] <- seq_vals
  }

  return (m)
}

#' Split read table into loci with no spanning reads, create consistent
#' haplotypes for each locus, and then combine to create full haplotypes
consistent_haplotypes_across_loci.read_table <- function(df_in, min_cover=500,
                                                         max_num_haplotypes=20000)
{
  sdf <- split_unlinked_loci.read_table(df_in, min_cover=min_cover, sort_loci = F)
  con_haps <- lapply(sdf, consistent_haplotypes.read_table, rm.na=T,
                     max_num_haplotypes=max_num_haplotypes)

  # if any of the haplotypes are NULL, meaning too many paths, return NULL
  if (any(sapply(con_haps, is.null)))
    return (NULL)

  merged_haps <- merge_haplotypes.read_table(con_haps)

  return (merged_haps)
}

#' Split read table into loci with no spanning single end reads, create
#' consistent haplotypes for each locus, and then combine to create full
#' haploptypes
consistent_haplotypes.read_table <- function(df_in, rm.na=F,
                                             max_num_haplotypes=20000)
{
  paired_split_out <- split_paired_ends.read_table(df_in)
  df <- paired_split_out$df
  pairs <- paired_split_out$pairs

  # if there are no paired reads, then do not split or filter
  if (nrow(pairs)==0) {
    haps <- consistent_haplotypes_single_end.read_table(df_in, rm.na=rm.na,
                                                        max_num_haplotypes=max_num_haplotypes)
    return (haps)
  }

  # split to create consistent haplotypes ignoring pairing of reads
  sdf <- split_unlinked_loci.read_table(df, min_cover=0, sort_loci = F)
  con_haps <- lapply(sdf, consistent_haplotypes_single_end.read_table, rm.na=rm.na,
                     max_num_haplotypes=max_num_haplotypes)

  # if any of the haplotypes are NULL, meaning too many paths, return NULL
  if (any(sapply(con_haps, is.null))) {
    return (NULL)
  }

  merged_haps <- merge_haplotypes.read_table(con_haps)

  # filter with paired reads
  consistent_haps <- filter_haplotypes_with_paired_ends.read_table(df_in,
                                                                    merged_haps)

  return (consistent_haps)
}

#' Quickly determine if the number of paths in an adjacency matrix exceed
#' a lower limit
#'
#' @param m adjacency matrix
#' @param paths_cutoff threshold for number of paths
#'
#' @details number of paths are determined by dynamic programming algorith
#' that is O(V+E).  The quick part is that once the number of paths
#' exceeds num_paths, the algorithm exits.   This is needed when the number
#' of edges is very large and we would waste time computing the total
#' number of paths.
#'
#' @return NA if number of paths exceeds num_paths, otherwise
#' returns the number of paths.
paths_exceed_limit.read_table <- function(m, paths_cutoff)
{
  nv <- nrow(m)

  # find current roots (vertices with no parents)
  current_roots <- which(sapply(1:nv, function(i) all(m[,i]==0)))
  # find leaves
  leaves <- which(sapply(1:nv, function(i) all(m[i,]==0)))

  # make global root and sink
  mplus <- matrix(0, nrow=nv+2, ncol=nv+2)
  mplus[2:(nv+1),2:(nv+1)] <- m
  mplus[1,current_roots+1] <- 1
  mplus[leaves+1,nv+2] <- 1

  g <- graph.adjacency(mplus, mode="directed")
  topo_order <- topo_sort(g, mode="out")
  nv <- length(topo_order)

  num_paths <- rep(NA, nv)
  num_paths[nv] <- 1

  for (i in (nv-1):1) {
    cv <- topo_order[i]

    child_ind <- which(mplus[cv,]==1)
    np <- num_paths[child_ind]

    # debug!
    if (any(is.na(np)))
      stop("bug in path counts")

    num_paths[cv] <- sum(np)

    # for now let's allow the number of paths to be returned
    if (num_paths[cv] > paths_cutoff)
      return (NA)
  }

  return (num_paths[1])
}

consistent_haplotypes_single_end.read_table <- function(df, rm.na=F,
                                                        max_num_haplotypes=20000)
{
  m <- adjacency_matrix.read_table(df, sparse=T)
  # if df is too large m will not be computed
  if (is.null(m))
    return (NULL)

  z <- paths_exceed_limit.read_table(m, max_num_haplotypes)
  cat("number of paths in locus", z, "\n")

  # quick check to see if the number of paths is too large
  if (is.na(z))
    return (NULL)

  nv <- nrow(m)

  # find current roots (vertices with no parents)
  current_roots <- which(sapply(1:nv, function(i) all(m[,i]==0)))
  # find leaves
  leaves <- which(sapply(1:nv, function(i) all(m[i,]==0)))

  # make global root
  mplus <- matrix(0, nrow=nv+1, ncol=nv+1)
  mplus[2:(nv+1),2:(nv+1)] <- m
  mplus[1,current_roots+1] <- 1
  leaves <- leaves + 1

  g <- graph.adjacency(mplus, mode="directed")

  # get paths between root and leaves
  print("computing all paths in read graph")
  paths <- all_simple_paths(g, from=1, to=leaves,
                            mode="out")
  cat("found", length(paths), "paths\n")

  if (length(paths) > max_num_haplotypes) {
    return (NULL)
  }

  paths_shifted <- lapply(paths, function(p) p[-1] - 1)
  haplotypes_all_paths <- paths_to_haplotypes.read_table(paths_shifted, df)

  # haplotypes may not be unique because different paths might
  # correspond to same haplotype
  haplotypes_s <- apply(haplotypes_all_paths, 1, paste, collapse="")
  haplotypes_s_unique <- names(table(haplotypes_s))
  # reconstructing haplotypes from string is complicated by NA, we
  # can't just split
  haplotypes_ind <- match(haplotypes_s_unique, haplotypes_s)
  haplotypes <- haplotypes_all_paths[haplotypes_ind,,drop=F]

  colnames(haplotypes) <- pos_names.read_table(df)

  if (rm.na) {
    ind <- apply(haplotypes, 1, function(h) all(!is.na(h)))
    haplotypes <- haplotypes[ind,,drop=F]
  }

  return (haplotypes)
}


filter_haplotypes_with_paired_ends.read_table <- function(df, haps)
{
  if ((ncol(df)-1) != ncol(haps))
    stop("inconsistent positions between read table and haplotype matrix")

  seq <- seq.read_table(df)
  pos <- all_pos.read_table(df)

  df_pos <- as.numeric(pos_names.read_table(df))
  pos_ind <- lapply(pos, function(cpos) match(cpos, df_pos))

  nreads <- nrow(df)
  nhaps <- nrow(haps)
  npos <- ncol(haps)
  covered <- matrix(F, nrow=nhaps, ncol=npos)

  for (readi in 1:nreads) {
    read_seq <- seq[[readi]]
    read_pos_ind <- pos_ind[[readi]]
    for (hapi in 1:nhaps) {
      matched <- haps[hapi,read_pos_ind] == read_seq
      if (all(matched))
        covered[hapi,read_pos_ind] <- T
    }
  }

  consistent <- apply(covered, 1, all)
  haps_out <- haps[consistent,,drop=F]

  return (haps_out)
}

############################
# methods to edit read_table objects

#' Filter read table for variant calls
#'
#' @param df read_table
#' @param variant_calls variant calls data.frame, see read_table() help
#' @param consensus_values consensus nucleotides at variable positions
#'
#' @details nrow of variant calls must equal length of conseunsus values
#'
#' @return a data.frame in which errors (i.e. not called varaints) are replaced
#' by the consensus value.
filter_true_variants.read_table <- function(df, variant_calls, consensus_values)
{
   df_vp <- pos_names.read_table(df)
   vc_vp <- as.character(variant_calls$pos)
   if (length(df_vp) != length(vc_vp))
     stop("read table variable positions do not match variant calls")
   else if (any(df_vp != vc_vp))
     stop("read table variable positions must mach variant call positions")

   if (length(consensus_values) != nrow(variant_calls))
     stop("consensus values do not match variant calls")

   nvp <- length(variant_calls$pos)
   vc <- as.matrix(dplyr::select(variant_calls, -pos))
   nucs <- colnames(vc)

   for (i in 1:nvp) {
     error_nucs <- nucs[!vc[i,]]
     ind <- which(is.element(df[,i+1], error_nucs))
     df[ind,i+1] <- consensus_values[i]
   }

   df <- regroup.read_table(df)
   return (df)
}

#' Clean read table
#'
#' @param df read table
#' @param min_count remove rows with count below
#' @param remove_outside_reads remove rows with all NA
#' @param remove_empty_cols remove cols with all NA
#' @param remove_non_variant_pos remove cols with only 1 value other than NA
#' @param remove_deletions replace "-" with NA
#' @param remove_partial_cover_reads remove rows with one or more NA
clean.read_table <- function(df, min_count=10,
                           remove_outside_reads=T,
                           remove_empty_cols=T,
                           remove_non_variant_pos=F,
                           remove_deletions=F,
                           remove_partial_cover_reads=F)
{
  df <- filter(df, count >= min_count)

  if (remove_empty_cols) {
    df_var <- dplyr::select(df, -count)
    empty_col <- apply(df_var, 2, function(x) {
      all(is.na(x))
    })
    df <- dplyr::select(df, which(c(T,!empty_col)))
  }

  if (remove_deletions) {
    for (i in 2:ncol(df)) {
      ind <- which(df[,i]=="-")
      if (length(ind)>0)
        df[ind,i] <- NA
    }
  }

  if (remove_non_variant_pos) {
    df_var <- dplyr::select(df, -count)
    no_variation <- apply(df_var, 2, function(x) {
      length(table(x))==1
    })
    df <- dplyr::select(df, which(c(T,!no_variation)))
  }

  if (remove_outside_reads) {
    df_var <- dplyr::select(df, -count)
    outside_ind <- apply(df_var, 1, function(x) all(is.na(x)))
    df <- filter(df, !outside_ind)
  }

  if (nrow(df)==0)
    return (df)

  if (remove_partial_cover_reads) {
    df_var <- dplyr::select(df, -count)
    partial_ind <- apply(df_var, 1, function(x) any(is.na(x)))
    df <- filter(df, !partial_ind)
  }

  if (nrow(df)==0)
    return (df)

  df <- regroup.read_table(df)

  return(df)
}

#' Remove reads that are below error noise.
#'
#' @param df read table
#' @param error_freq per position error rate
#' @param sig significance level at which to filter out errors.
#'
#' @return a read table with errors that are below noise threshold
#' removed
error_filter.read_table <- function(df, error_freq, sig)
{
   temp <- templates.read_table(df)
   ntemp <- length(temp)

   # row in df corresponding to each template
   temp_read_ind <- template_indices.read_table(temp, df)
   # counts for each row in df corresponding to each template
   temp_counts <- lapply(temp_read_ind, function(ind) df$count[ind])
   # number of variable positions in each template
   temp_num_pos <- apply(temp, 1, sum)


   keep_ind <- Map(function(ind, counts, npos, i) {
     total <- sum(counts)

     #if (total < 100)
    #   return (NULL)
     # if sig==0, then no errors and include all reads
     if (sig==0)
       return (ind)

     mu <- total*error_freq
     # do we reject a single read group
     if (ppois(1, mu) > (1 - sig)) {
       return (NULL)
     }

     above_noise <- ppois(counts, mu) > (1 - sig)
     above_noise_ind <- ind[above_noise]
     above_noise_counts <- counts[above_noise]

     #if (sum(above_noise_counts) < 100)
      # return (NULL)

     return (above_noise_ind)
   }, temp_read_ind, temp_counts, temp_num_pos,
   1:length(temp_read_ind))

   keep_ind_vec <- unlist(keep_ind)
   df_filtered <- df[keep_ind_vec,]
  # bad_ind <- setdiff(1:nrow(df), keep_ind_vec)

   return (df_filtered)
}

#' Merges identical read over an edited read table
#'
#' @param df a read_table object
#'
#' @return a read_table object
#'
#' @details After removing columns for positions that are
#' not of interest, a read table can contain multiple rows
#' corresponding to read groups that are identical at the remaining positions.
#' This functions joins those reads and returns a read_table
#' with unique read groups.
regroup.read_table <- function(df)
{
  nreads <- nrow(df)
  df_nucs <- dplyr::select(df, -count)

  # read_strings uniquely identify each read
  read_strings <- apply(df_nucs, 1, paste, collapse="/")
  df <- mutate(df, read_strings=read_strings)

  # split on read_string
  regrouped_df <- ddply(df, .(read_strings), function(rdf) {
    count <- sum(rdf$count)
    rdf$count[1] <- count
    rdf[1,]
  })
  regrouped_df <- dplyr::select(regrouped_df, -read_strings)

  return (regrouped_df)
}

#' Split a read_table into multiple read_tables representing unlinked loci
split_unlinked_loci.read_table <- function(df, sort_loci=T, min_cover=1000)
{
  linkage_ind <- unlinked_pos.read_table(df, min_cover=min_cover)
  # create a vector of positions with a , separating linked
  # positions and a / separting unlinked
  linkage_string <- mapply(function(pos, linkage_break, entry_num) {
    if (entry_num==1)
      return (pos)

    if (linkage_break)
      return (paste("/", pos, sep=""))
    else
      return (paste(",", pos, sep=""))
  }, names(linkage_ind), linkage_ind, 1:length(linkage_ind))
  linkage_string <- paste(linkage_string, collapse="")

  # unpack linkage string
  linked_pos <- strsplit(linkage_string, split="/")[[1]]
  linked_pos <- lapply(linked_pos, function(cpos)
    strsplit(cpos, split=",")[[1]])

  df_unlinked <- lapply(linked_pos, function(lp) {
    cdf_in <- df[,c("count", lp)]
    cdf <- clean.read_table(cdf_in, min_count=0,
                            remove_non_variant_pos = F,
                            remove_deletions = F)
    return (cdf)
  })

  if (sort_loci) {
    npos <- sapply(df_unlinked, ncol)
    df_unlinked <- df_unlinked[order(npos, decreasing=T)]
  }

  return (df_unlinked)
}

#' Split a read table into two read tables.
#'
#' @return A list of two read tables labeled df1 and df2.  If the read table
#' has only 1 position, then NULL is returned.
split.read_table <- function(df)
{
  pos <- as.numeric(pos_names.read_table(df))
  npos <- length(pos)
  if (npos==1)
    return (NULL)

  pos1 <- pos[1:round(npos/2-.01)]
  pos2 <- setdiff(pos, pos1)
  df1 <- subset.read_table(df, pos=pos1)
  df2 <- subset.read_table(df, pos=pos2)

  return (list(df1=df1, df2=df2))

}

join_unlinked_loci.read_table <- function(df_list)
{
  if (class(df_list) != "list")
    stop("df_list must be a list")

  nlist <- length(df_list)
  if (nlist==0)
    return (NULL)

  if (nlist==1)
    return (df_list[[1]])

  df_out <- df_list[[1]]
  for (i in 2:nlist)
    df_out <- join_unlinked_loci_pair.read_table(df_out, df_list[[i]])

  return (df_out)
}

join_unlinked_loci_pair.read_table <- function(df1, df2)
{
  df1_nrow <- nrow(df1)
  df2_nrow <- nrow(df2)

  all_count <- c(df1$count, df2$count)

  df1_nucs <- dplyr::select(df1, -count)
  df2_nucs <- dplyr::select(df2, -count)

  df1_addon <- data.frame(matrix(as.character(NA), nrow=nrow(df2), ncol=ncol(df1_nucs)),
                          stringsAsFactors = F)
  names(df1_addon) <- names(df1_nucs)
  df2_addon <- data.frame(matrix(as.character(NA), nrow=nrow(df1), ncol=ncol(df2_nucs)),
                          stringsAsFactors = F)
  names(df2_addon) <- names(df2_nucs)

  df1_full_nucs <- rbind(df1_nucs, df1_addon)
  df2_full_nucs <- rbind(df2_addon, df2_nucs)

  df_full <- cbind(all_count, df1_full_nucs, df2_full_nucs)
  names(df_full) <- c("count", names(df1_nucs), names(df2_nucs))

  return (df_full)
}

#' Limit read table to a subset of positions
subset.read_table <- function(df, pos=NULL, start_pos=NULL,
                              end_pos=NULL)
{
  all_pos <- as.numeric(setdiff(names(df), "count"))

  if (is.null(pos)) {
    if (is.null(start_pos) | is.null(end_pos))
      stop("if pos is NULL then start_pos and end_pos must be given")

    pos <- all_pos[all_pos >= start_pos & all_pos <= end_pos]
    pos <- c("count", as.character(pos))
  }
  else
    pos <- unique(c("count", as.character(pos)))

  pos <- intersect(pos, names(df))

  df_new <- df[,pos,drop=F]
  df_new <- clean.read_table(df_new, remove_non_variant_pos=F,
                             min_count = 0)

  return (df_new)

}

################
plot_cover.read_table <- function(df, min_cover=1000)
{
  # for viewability, remove some reads
  #df <- clean.read_table(df, min_count=50, remove_non_variant_pos = T)

  # sort the reads so that we can visually compare reads of equal length
  n_nucs <- apply(df, 1, function(rr) sum(!is.na(rr)))
  start_pos <- start_pos.read_table(df)
  ind_order <- order(start_pos, n_nucs, decreasing = F)
  df <- df[ind_order,]

  df_nucs <- dplyr::select(df, -count)
  npos <- ncol(df_nucs)
  pos <- names(df_nucs)
  nreads <- nrow(df_nucs)

  df_rect <- adply(df_nucs, 1, function(cdf) {
    active <- which(!is.na(cdf))
   # seq <- paste(cdf[active], collapse="")
    start <- min(active)
    end <- max(active)

  #  gaps <- setdiff(start:end, active)
  #  if (length(gaps)>0) browser()

    seq_split <- ifelse(is.na(cdf[start:end]), "=", cdf[start:end])
    seq <- paste(seq_split, collapse="")

    return (data.frame(seq=seq, start=start, end=end,
                       stringsAsFactors = F))
  }, .expand=F)

  covered_m <- matrix(F, nrow=nrow(df), ncol=npos)
  height <- rep(NA, nreads)
  for (r in 1:nreads) {
    place <- 1
    read_pos <- df_rect$start[r]:df_rect$end[r]
    while(any(covered_m[place,read_pos]))
      place <- place + 1
    height[r] <- place

    covered_m[place,read_pos] <- T
  }

  df_rect <- mutate(df_rect, height=height)
  df_rect <- mutate(df_rect, count=(df$count))


  df_text <- adply(df_rect, 1, function(cdf) {
    pos <- cdf$start:cdf$end + 1/2
    y <- cdf$height - 1/2
    seq <- strsplit(cdf$seq, split="")[[1]]
    if (length(seq) != length(pos))  {
      stop("non-paired read!")
    }
    data.frame(pos=pos, y=y, seq=seq,
               stringsAsFactors = F)
  })

  linkage_v <- unlinked_pos.read_table(df, min_cover = min_cover)
  split_pos <- which(linkage_v)
  if (length(split_pos) > 0)
    df_linkage <- data.frame(x=split_pos, xend=split_pos,
                           y=0, yend=max(df_rect$height))
  else
    df_linkage <- NULL

  p <- ggplot()

  p <- p + geom_rect(mapping=aes(xmin=start, xmax=end+1,
                             ymin=height-1, ymax=height,
                             fill=count),
                     data=df_rect, col="white", size=1)

  p <- p + geom_text(mapping=aes(x=pos, y=y, label=seq),
                     data=df_text, col="white", size=5)

  if (!is.null(df_linkage))
    p <- p + geom_segment(mapping=aes(x=x, xend=xend, y=y, yend=yend),
                        data=df_linkage, col="yellow", size=2)



  p <- p + scale_y_continuous(expand = c(0,0),
                              breaks=1,
                              labels="")


  p <- p + scale_x_continuous(expand = c(0,0),
                              breaks=(1:npos)+.5,
                              labels=as.numeric(pos))
  p <- p + theme(axis.text=element_text(size=11))
  p <- p + xlab("position")


  return (p)
}

plot_graph.read_table <- function(df)
{
  paired_split_out <- split_paired_ends.read_table(df)
  df <- paired_split_out$df

  counts <- df$count
  start_pos <- start_pos.read_table(df)
  unique_start_pos <- sort(unique(start_pos))
  spread <- 5
  start_ind <- spread*sapply(start_pos, function(sp) which(unique_start_pos==sp))
  unique_start_ind <- sort(unique(start_ind))

  seq <- seq.read_table(df)
  m <- adjacency_matrix.read_table(df)

  nrg <- length(seq)

  # vertex data.frame
  v_df <- data.frame(group=1:nrg,
                     start_ind=start_ind,
                     start_pos=start_pos,
                     count=counts,
                     seq=sapply(seq, paste, collapse=""))

  v_df <- ddply(v_df, .(start_pos), function(sp_df) {
    xmin <- sp_df$start_ind[1]
    ymin <- seq(from=1, to=10, length.out=nrow(sp_df))
    info_s <- paste(sp_df$group, sp_df$seq, sp_df$count, sep=" ")
    mutate(sp_df, xmin=xmin, xmax=xmin+(.5*spread),
           ymin=ymin, ymax=ymin+.5, info_s=info_s)
  })

  #edge data.frame
  ne <- sum(as.numeric(m))
  e_m <- matrix(nrow=ne, ncol=4)
  colnames(e_m) <- c("x", "xend", "y", "yend")

  counter <- 1
  for (i in 1:nrow(m))
    for (j in 1:ncol(m)) {
      if (m[i,j]==1) {
        i_ind <- min(which(v_df$group == i))
        j_ind <- min(which(v_df$group == j))
        e_m[counter,] <- c(v_df$xmax[i_ind],
                           v_df$xmin[j_ind],
                           (v_df$ymin[i_ind]+v_df$ymax[i_ind])/2,
                           (v_df$ymin[j_ind]+v_df$ymax[j_ind])/2)
        counter <- counter+1
      }
    }

  e_df <- as.data.frame(e_m)

  p <- ggplot()
  p <- p + geom_rect(mapping=aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax),
                     data=v_df, fill="white", size=1, col="black")
  p <- p + geom_text(mapping=aes(x=(xmin+xmax)/2, y=(ymin+ymax)/2, label=info_s),
                     data=v_df, size=5, vjust=-2.3)
  if (nrow(e_df) > 0)
    p <- p + geom_segment(mapping=aes(x=x, y=y, xend=xend, yend=yend),
                        data=e_df,
                        arrow = arrow(length = unit(0.03, "npc")),
                        size=1.2)

  p <- p + scale_y_continuous(expand = c(0,0),
                              breaks=NULL,
                              labels=NULL,
                              limit=c(0,12))
  p <- p + scale_x_continuous(expand = c(0,0),
                              breaks=unique_start_ind,
                              labels=unique_start_pos,
                              limit=c(min(unique_start_ind)-1,
                                      max(unique_start_ind)+.6*spread))

  p <- p + xlab("start position") + ylab("")
  p <- p + theme(axis.text=element_text(size=20))
  return (p)
}

