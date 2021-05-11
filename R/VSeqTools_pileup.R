#' Read a VSeqTools_pileup object
#' 
#' @param filename A csv file 
#' 
#' @return a VSeqTools_pileup object
read_BAM_pileup <- function(filename)
{
  df <- read.table(filename, header=T, sep=",", check.names = F,
                   stringsAsFactors = F)
  
  class(df) <- c("VSeqTools_pileup", "data.frame")
  
  return (df)
}


#' Calculate pileup for each position of reference from BAM file
#'
#' Using the Rsamtools pileup function, the counts of A,C,G,T,i,d
#' are determined at each position on reference
#'
#' @param bam_file bam file
#' @param bai_file If not passed, then assumed to be same as bam_file
#' with .bai appendix.
#' @param start_pos Start position of pileup relative to reference.
#' @param end_pos End position of pileup relative to reference.
#' @param max_depth maximum number of reads from which to form the pileup.
#' This parameter is passed to RSamTools PileupParam function.
#' @param min_base_quality Passed to PileupParam
#' @param min_mapq Passed to PileupParam
#' @param strand If BAM file is paired end, which strand to use: "plus" or "minus".  Should
#' be set to NULL if all reads to be used or if BAM file is not paired.
#'
#' @return a VSeqTools_pileup object which inherits data.frame
#' and has columns: pos, A, C, G, T, i, d, Coverage.
#' pos gives the position on the reference.  Coverage is total reads
#' at that position.  The rest of the columns give read freq at that position
#' with the particular nt value.
#'
BAM_pileup <- function(bam_file,
                  bai_file=NULL,
                  max_depth=2000,
                  min_base_quality=30,
                  min_mapq=30,
                  strand=NULL)
{

  if (is.null(bai_file)) {
    bai_file <- paste(bam_file, ".bai", sep="")
  }

  if (!is.null(strand)) {
    p_param <- PileupParam(max_depth=max_depth,
                           min_base_quality=min_base_quality,
                           min_mapq=min_mapq,
                           include_insertions = T,
                           distinguish_strands = T)

    p <- pileup(bam_file, index=bai_file, pileupParam = p_param)
    if (strand=="plus")
      p <- filter(p, strand=="+")
    else if (strand=="minus")
      p <- filter(p, strand=="-")
    else
      stop("strand argument must be plus, minus, or NULL")

  }
  else {
   p_param_no_strand <- PileupParam(max_depth=max_depth,
                         min_base_quality=min_base_quality,
                         min_mapq=min_mapq,
                         include_insertions = T,
                         distinguish_strands = F)
   p <- pileup(bam_file, index=bai_file,
                  pileupParam = p_param_no_strand)
  }

  min_pos <- min(p$pos)
  max_pos <- max(p$pos)
  

  get_variant_df <- function(pileup_df) {

    template <- data.frame(pos=NA, A=0, C=0, G=0, T=0, d=0, i=0,
                           coverage=0)
    variant_df <- ddply(pileup_df, .(pos), function(df) {
      coverage <- sum(df$count)
      freq <- df$count/coverage
      nucs <- as.character(df$nucleotide)
      nuc_type <- sapply(nucs, function(cn) {
            switch(cn, A="A", C="C", G="G", T="T", "+"="i", "-"="d")
      })

      sum_df <- template
      sum_df[,nuc_type] <- freq
      sum_df$pos <- df$pos[1]
      sum_df$coverage  <- coverage

      return (sum_df)
    })
  }

  df <- get_variant_df(p)

  class(df) <- c("VSeqTools_pileup", "data.frame")

  return (df)
}

get_pileup_nt_names <- function()
{
  c("A", "C", "G", "T", "d", "i")
}

##########################################
# - methods that use pileup

#' Returns a consensus sequence given a VSeqTools_pileup object
#'
#' Given a pileup, a consensus sequence is determined at each position
#' on the reference.  The consensus is returned as a character vector
#' with possible entries A, C, G, T, d, i.
#'
#' @param pu a VSeqTools_pileup object returned by BAM_pileup
#'
#' @return
#' A character vector with entries A,C,G,T,d,i and length equal to
#' positions on reference.
#'
#' Warning:  if the VSeqTools_pileup object is formed from a BAM file
#' with gaps, the consensus will skip missed positions.
#'
consensus <- function(pu)
{
  if (!is.element("VSeqTools_pileup", class(pu))) {
    stop("consensus expects VSeqTools_pileup object")
  }

  # treat pileup as data.frame
  df <- pu

  nucs <- get_pileup_nt_names()

  consensus <- sapply(1:nrow(df), function(i) {
    ind <- which.max(df[i,nucs])
    return (nucs[ind])
  })

  return (consensus)
}

#' Returns variant calls given a VSeqTools_pileup object
#'
#' Given a pileup, variants are called accorging to a frequency
#' cutoff (min_freq) or a significance level(sig)
#'
#' @param pu a VSeqTools_pileup object returned by BAM_pileup
#' @param min_freq Minimum frequency above which a variant is called.
#' @param sig Significance levels above which a variant is called
#' @param heavy_tail If T, use betabinomial error distribution, otherwise
#' use poisson error distribution. Only relavent when !is.null(sig)
#'
#' @return
#' A VSeqTools_variant_call object which inherits the data.frame class.
#' VSeqTools_variant_call objects are
#' data.frames with
#' columns pos, A, C, G, T, d, i.  pos gives the position on
#' the reference.  Entries in other columns are logical depending
#' on whether associated variant/deletion/insertion is called
#'
#' @details Exactly one of min_freq and sig must be non-NULL.  This
#' determines whether cutoff or poisson algorithm is used.
#'
variant_calls <- function(pu, min_freq=NULL,
                          sig=NULL, heavy_tail=NULL)
{
  if (!is.element("VSeqTools_pileup", class(pu))) {
    stop("variant_calls expected VSeqTools_pileup object")
  }

  if (is.null(min_freq) & is.null(sig))
    stop("either min_freq or sig must be specified")
  if (!is.null(min_freq) & !is.null(sig))
    stop("min_freq and sig cannot both be non-NULL")

  if (is.null(sig))
    return (variant_calls_cutoff(pu, min_freq=min_freq))
  else
    return (variant_calls_sig(pu, sig=sig, heavy_tail=heavy_tail))


#   nucs <- get_pileup_nt_names()
#
#   df <- pu
#   df_nucs <- as.matrix(df[,nucs])
#
#   df_vc <- matrix(as.numeric(df_nucs) > min_freq, nrow=nrow(df_nucs),
#                   ncol=ncol(df_nucs))
#   df <- cbind(df$pos, df_vc)
#   df <- data.frame(df)
#   names(df) <- c("pos", nucs)
#
#   class(df) <- c("VSeqTools_variant_call", "data.frame")
#
#   return (df)
}

#' Get variant call positions from a VSeqTools_variant_call
#' object
#'
#' @param vc VSeqTools_variant_call object
#'
#' @return a numeric vector giving positions relative to reference at which
#' variants are called.
get_variant_call_pos <- function(vc)
{
  nucs <- setdiff(names(vc), c("pos"))
  vc_nucs <- vc[,nucs]

  ind <- which(apply(vc_nucs, 1, function(cr) sum(cr)>1))
  vals <- sapply(ind, function(cind) {
    cvals <- nucs[which(vc_nucs[cind,] != 0)]
    paste(cvals, collapse="")
  })

  out <- vc$pos[ind]
  names(out) <- vals

  return (out)
}

#' Estimate error rate of pileup
#'
#' @param pu a VSeqTools_pileup object returned by BAM_pileup
#' @param split_by_nuc If T, an error rate is returned for each nucleotide,
#' otherwise a single error rate is returned representing sum of all possible errors
#'
#' @return A vector of error rates.
get_error_rate <- function(pu, split_by_nuc=F)
{
  nucs <- get_pileup_nt_names()

  # construct an error rate for each nucleotide
  errors <- lapply(nucs, function(cn) {
    freqs <- as.numeric(pu[,cn])
    # throw out anything > 3% as cutoff
    ind <- freqs < .03

    data.frame(error_freq=freqs[ind], coverage=pu$coverage[ind])
  })
  names(errors) <- nucs

  mean_error_rate <- sapply(errors, function(cdf) {
    p <- sum(cdf$error_freq*cdf$coverage)/sum(cdf$coverage)
    return (p)
  })
  names(mean_error_rate) <- nucs

  if (!split_by_nuc)
    mean_error_rate <- sum(mean_error_rate)

  return (mean_error_rate)
}

#' Calls nt positions at which a BAM file differs from a base BAM file
#'
#' The pileup of the first BAM file is compared to the consensus of the
#' second BAM file and positions at which there are disagreements are called
#'
#' @param BAMfile BAM file
#' @param b_BAMfile base BAM file to which BAMfile is to be compared
#' @param min_divergenece_freq minimum frequency for which divergence is called
#'
#' @details insertions are ignored in calling divergent positions
#'
#' @return
#' A numeric vector of called nt positions.
divergence_variant_calls <- function(BAMfile, b_BAMfile,
                                min_divergence_freq=0.02)
{
  # get consensus for base
  bpu <- BAM_pileup(b_BAMfile)
  consensus <- VSeqTools::consensus(bpu)

  pu <- BAM_pileup(BAMfile)

  if (any(bpu$pos != pu$pos))
    stop("different positions in base and comparison pileups")
  npos <- length(bpu$pos)

  divergent_pos <- sapply(1:npos, function(i) {
    c_consensus <- consensus[i]
    diverged_freq <- 1 - pu[i,c_consensus]

    if (diverged_freq < min_divergence_freq)
      return (F)
    return (T)
  })

  pos <- bpu$pos[divergent_pos]

  return (pos)
}


########################################################################
# variant_call algorithms

#' Returns variant calls given a VSeqTools_pileup object
#'
#' Given a pileup, variants with variability above min_freq are
#' identified.  This is a simple variant caller which is much faster
#' than VPhaser2, but cannot call low frequency variants with precision
#' since phasing information is not used
#'
#' @param pu a VSeqTools_pileup object returned by BAM_pileup
#' @param min_freq Minimum frequency above which a variant is called.
#'
#' @return
#' A VSeqTools_variant_call object which inherits the data.frame class.
#' VSeqTools_variant_call objects are
#' data.frames with
#' columns pos, A, C, G, T, d, i.  pos gives the position on
#' the reference.  Entries in other columns are logical depending
#' on whether associated variant/deletion/insertion is called
#'
#' Warning:  Currently uses the df data.frame in VSeqTools_pileup
#' object.  Needs to be extended to allow for different strands.
#'
variant_calls_cutoff <- function(pu, min_freq=.03)
{
  if (!is.element("VSeqTools_pileup", class(pu))) {
    stop("variant_calls expected VSeqTools_pileup object")
  }

  nucs <- get_pileup_nt_names()

  df <- pu
  df_nucs <- as.matrix(df[,nucs])

  df_vc <- matrix(as.numeric(df_nucs) > min_freq, nrow=nrow(df_nucs),
                  ncol=ncol(df_nucs))
  df <- cbind(df$pos, df_vc)
  df <- data.frame(df)
  names(df) <- c("pos", nucs)

  class(df) <- c("VSeqTools_variant_call", "data.frame")

  return (df)
}




#' Returns variant calls given a VSeqTools_pileup object
#'
#' Given a pileup, implements a poisson significance test to
#' call variants.
#'
#' @param pu a VSeqTools_pileup object returned by BAM_pileup
#' @param sig Significance level of poisson test
#' @param heavy_tail If T use betabinomial error distribution, otherwise
#' use poisson error distribution.
#'
#' @return
#' A VSeqTools_variant_call object which inherits the data.frame class.
#' VSeqTools_variant_call objects are
#' data.frames with
#' columns pos, A, C, G, T, d, i.  pos gives the position on
#' the reference.  Entries in other columns are logical depending
#' on whether associated variant/deletion/insertion is called
#'
#' @references Wang,C. et al. (2007) Characterization of mutation
#' spectra with ultra-deep pyrosequencing: application to HIV-1 drug
#' resistance. Genome Res., 17, 1195–1201,
#' @references Gerstung et al. (2011) Reliable detection of subclonal
#' single-nucleotide variants in tumour cell populations.  Nature
#' Communications.
#'
variant_calls_sig <- function(pu, sig=.01, heavy_tail=T)
{
  if (!is.element("VSeqTools_pileup", class(pu))) {
    stop("variant_calls expected VSeqTools_pileup object")
  }

  nucs <- get_pileup_nt_names()
  nnucs <- length(nucs)

  df <- pu
  df_nucs <- as.matrix(df[,nucs])

  # for each nucleotide, we determine the error frequency by
  # averaging over position, but ignoring errors greater than
  # 3% to avoid including true variation.
  consensus <- apply(df_nucs, 1, which.max)

  # construct an error rate for each nucleotide
  errors <- lapply(nucs, function(cn) {
    freqs <- as.numeric(pu[,cn])
    # throw out anything > 3% as cutoff
    ind <- freqs < .03

    data.frame(error_freq=freqs[ind], coverage=pu$coverage[ind])
  })
  names(errors) <- nucs

  mean_error_rate <- sapply(errors, function(cdf) {
    p <- sum(cdf$error_freq*cdf$coverage)/sum(cdf$coverage)
    p <- max(p, 1E-5)
    return (p)
  })
  cat("mean_error rate", mean_error_rate, "\n")

  # if heavy_tail=T, fit beta to error_rates
  if (heavy_tail) {
    # fitdistr does not converge so use other method
    #ft <- fitdistr(error_rate, "beta", list(shape1=10, shape2=1/mean_error_rate))
    #alpha <- ft$estimate[1]
    #beta <- ft$estimate[2]

    # parameterize beta(alpha, beta) by (alpha-1) = Rp, (beta-1) = R(1-p)
    # we assume alpha,beta > 1 and p = mean_error_rate.  p will be
    # the max of the beta pdf.  This gives R >= 1/p.

    # using MLE leads to too tight of a distribution.  better is to capture
    # the width of the error rate.  so we focus on the 5% and 95% quantiles
    fit_f <- function(R, quants, p) {
       quants_hat <- qbeta(c(.05, .95), R*p, R*(1-p))
       return (sum((quants-quants_hat)^2))
    }

    bb_par <- Map(function(cdf, p) {
      quants <- quantile(cdf$error_freq, prob=c(.05, .95))

      op_out <- optimize(fit_f, interval=c(1/p, 1000/p), quants, p, maximum=F)
      R <- op_out$minimum
      alpha <- 1 + R*p
      beta <- 1 + R*(1-p)

      bb_p <- alpha/(alpha+beta)
      bb_dispersion <- alpha+beta
      return (list(bb_p=bb_p, bb_dispersion=bb_dispersion))
    }, errors, mean_error_rate)
    names(bb_par) <- nucs

  }

  # go through each position and nucleotide and call variation
  # based on coverage,
  df_vc <- sapply(1:nnucs, function(i) {
    p <- df_nucs[,i]
    coverage <- df$coverage
    successes <- round(p*coverage)
    poisson_mean <- mean_error_rate[i]*coverage

    bb_p <- bb_par[[i]]$bb_p
    bb_dispersion <- bb_par[[i]]$bb_dispersion

    if (heavy_tail) {
      pval <- 1 - pbetabinom(successes, coverage, bb_p, bb_dispersion) +
              dbetabinom(successes, coverage, bb_p, bb_dispersion)
    } else {
      pval <- ppois(successes, poisson_mean, lower.tail=F) +
        dpois(successes, poisson_mean)
    }

    t <- pval < sig
    return (t | (consensus == i))
  })

  df <- cbind(df$pos, df_vc)
  df <- data.frame(df)
  names(df) <- c("pos", nucs)

  class(df) <- c("VSeqTools_variant_call", "data.frame")

  return (df)
}


#############################################################
#' Converts reads to a sequence in the reference space
#'
#' Given a vector of read sequences and their corresponding cigars, returns the sequences
#' in the reference space by removing soft clippings, removing insertions (relative to
#' the reference), and inserting
#' "-" for deletions (relative to the reference)
#'
#' @param seq Read sequences
#' @param cigar Read cigars
#'
#' @return A vector of sequences in the reference space.  Insertions are ignored since they
#' are not in the reference space
#'
create_refspace_seq <- function(seq, cigar)
{
  # creates character vector of sequences in reference space.
  cat("number seq ", length(seq), "\n")

  # seq = vectors of DNA strings
  # cigar = cigar strings

  # extracts the letters from cigar strings ("150M20S" --> c("M", "S")
  ops <- explodeCigarOps(cigar)

  # returns a IRangesList object.  Each IRanges describes the cigar on the reference
  #   meaning that soft clippings and insertions are ignored, width==0
  ref_pos <- cigarRangesAlongReferenceSpace(cigar)
  # describes the cigar on the sequence strings meaning that deletions are ignored
  query_pos <- cigarRangesAlongQuerySpace(cigar)


  seq_paste <- mapply(function(cop, ref_qwidth, query_qwidth,
                               query_start, query_end, cseq, i) {
    # if qwidth > 0 in ref_pos, we check qwidth in query_pos
    #  if query_pos qwidth == 0, we have a deletion, so we insert "---" of length qwidth
    #          otherwhise, we have a match, so we insert substr(seq, start, end)
    str <- sapply(1:length(cop), function(i) {
      if (ref_qwidth[i]==0)  # soft clipping or insertion
        return ("")
      if (query_qwidth[i]==0)
        return (paste(rep("-", ref_qwidth[i]), collapse=""))

      return (substr(cseq, query_start[i], query_end[i]))
    })

    if (i %% 1000 == 0)
      cat("converting seq to ref space:", i, "of", length(seq), "\n")

    return (paste(str, collapse=""))
  }, ops, width(ref_pos), width(query_pos), start(query_pos),
  end(query_pos), as.character(seq), 1:length(seq), SIMPLIFY=T)


  return (seq_paste)
}

