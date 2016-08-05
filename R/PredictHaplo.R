
#' Sample from dirichlet distribution
#' @export
rdirichlet<-function(n,a)
  ## pick n random deviates from the Dirichlet function with shape
  ## parameters a

  # copied from
  # https://stat.ethz.ch/pipermail/r-help/2000-December/009561.html
{
  l<-length(a);
  x<-matrix(rgamma(l*n,a),ncol=l,byrow=TRUE);
  sm<-x%*%rep(1,l);
  sample <- x/as.vector(sm);

  return (sample)
}

get_debug_df <- function()
{
  df <- data.frame(count=c(200, 800,
                           250, 750),
                   p1=c("A", "C" , NA, NA),
                   p2=c(NA, NA, "G", "T"),
                   stringsAsFactors = F)
  names(df) <- c("count", 1:2)
  #df <- df[3:4,]
  return (df)
}

#' Constructor for a PredictHaplo object
#'
#' Implements a global haplotype reconstruction using PredictHaplo algorithm.
#' Given a read table and K haplotypes, constructs the pi, phi, and c vectors
#' discussed in section 3.3.1 of the Prabhakaran et al 2014 paper
#'
#' @param df read table data.frame
#' @param K number of haplotypes to infer
#' @param pi_prior A prior for pi.  If NULL default prior is used.
#' @param gamma parameter for bayesian prior
#' @param alpha parameter for bayesian prior
#' @param iterations number of iterations of Gibbs sampling to preform
#'
#' @return A PredictHaplo object which inherits list and contains the entries
#' pi, phi, cvals, df matching their definitions in Prabhakaran et al.  (cvals instead of c).
#' df is the read table and is stored as part of the PredictHaplo object for convenience of
#' post-processing.
#'
#' @details pi is stored as a numeric vector of length K.  phi is a 3-d array with
#' phi[k,i,l] = probability of drawing nuc l at position i from haplotype k.
#' cvals is a matrix with cvals[k,j] = number of read of type j assigned to haplotype k.
#' Reads of type j are described by the jth row of the read table df.
#'
#' @references
#' Prabhakaran, Rey, Zagordi et al, 2014, IEE/ACM TRANSACTIONS ON
#' COMPUTATIONAL BIOLOGY AND BIOINFORMATICS
#' @export
PredictHaplo <- function(df, K=2, pi_prior=NULL, gamma=2, alpha=4,
                         max.iterations=1000,
                         track_posterior=F)
{
  counts <- df$count
  npos <- ncol(df)-1
  nucs <- c("A", "C", "G", "T")
  nnucs <- length(nucs)
  df_nucs <- select(df, -count)

  # create matrices that identify which reads have a given nucleotide at which position
  df_vec <- as.character(as.matrix(df_nucs))
  df_vec <- ifelse(is.na(df_vec), "x", df_vec)
  m_reads_list <- lapply(nucs, function(cnuc)
    matrix(as.numeric(df_vec == cnuc), nrow=nrow(df), ncol=ncol(df)-1)
  )
  names(m_reads_list) <- nucs

  posterior_prob <- rep(NA, max.iterations)

  # draw from prior
  # pi[k] = frequency of haplotype k
  if (is.null(pi_prior))
    pi_prior <- rep(1/K, K)
  pi <- rdirichlet(1, gamma*pi_prior)

  # phi[k,i,l] = probability of drawing nuc l at position i from haplotype k
  phi <- array(as.numeric(NA), dim = c(K, npos, 4),
               dimnames = list(NULL, names(df_nucs), nucs))
  for (k in 1:K)
    for (i in 1:npos)
      phi[k,i,] <- rdirichlet(1, rep(alpha/4, 4))

  # cvals[k,j] = number of reads from row j in df that are assigned to
  # haplotype k
  cvals <- assign_reads.PredictHaplo(df, pi, phi)

  #fit <- PredictHaplo_fit(list(df=df, phi=phi, pi=pi, cvals=cvals))$total

  for (iter in 1:max.iterations) {

    if (iter %% 50 == 0)
      cat("iteration", iter, "\n")

    # update pi based on assignments
    num_haplo_assignments <- rowSums(cvals)
    pi <- rdirichlet(1, num_haplo_assignments + gamma*pi_prior/K)

    # update phi based on pi and assignments
    # z[[l]][k,i] = num of reads assigned to haplotype k with nuc l at pos i
    z_list <- lapply(nucs, function(cnuc) cvals %*% m_reads_list[[cnuc]])
    for (k in 1:K)
      for (i in 1:npos) {
        z <- sapply(1:nnucs, function(l) z_list[[l]][k,i] )
        phi[k,i,] <- rdirichlet(1, alpha/4 + z)
      }

    # update cvals (assignments) based on pi and phi
    cvals <- assign_reads.PredictHaplo(df, pi, phi)

    if (track_posterior)
      posterior_prob[iter] <- posterior_log_prob.PredictHaplo(df, pi, phi,
                                                            cvals, alpha, gamma)

#     if (iter %% 200 == 0) {
#       nfit <- fit.PredictHaplo(list(df=df, phi=phi, pi=pi, cvals=cvals))$total
#       cat(fit, nfit, "\n")
#       if (nfit > fit - .001)
#         break
#       else
#         fit <- nfit
#     }
  }

  if (!track_posterior) {
    posterior_prob <- posterior_log_prob.PredictHaplo(df, pi, phi,
                                                          cvals, alpha, gamma)
  }

  ph <- list(df=df, phi=phi, pi=pi, cvals=cvals,
             posterior=posterior_prob)

  class(ph) <- c("PredictHaplo", "list")
  return (ph)
}

#' Helper function that calculates probability of read assignments to
#' haplotypes
#'
#' @param df A read table
#' @param pi Numeric vector of haplotype frequencies
#' @param phi A haplotype probability matrix from a PredictHaplo object
#'
#' @return A (K x M) matrix with K equaling the number of haplotypes and
#' M equaling the number of read groups (rows in df).
#'
#' @details Entry k,m of the returned matrix is the probability that
#' a read in read group m is assigned to haplotype k.  Some reads
#' may not be assignable, specifically when a read contains a nucleotide
#' at a position for which all phi have probability 0.  In this case, NaN
#' is returned as the probability.
#' @export
read_assignment_prob_matrix.PredictHaplo <- function(df, pi, phi)
{
  df_nucs <- select(df, -count)
  K <- ncol(pi)

  if (K==1)
    return (matrix(1, nrow=1, ncol=nrow(df_nucs)))

  # for each read type, generate assignment prob given by
  #  pi[k] * \prod_s phi[k, s, r_s]
  # where s is an index over the positions in the read that vary
  # and r_s is the nucleotide at position s of the read
  marginal_prob <- apply(df_nucs, 1, function(nucs) {

    ind <- which(!is.na(nucs))
    pos_vals <- names(df_nucs)[ind]

    if (length(ind)==0)
      return (rep(1/K, K))

    nuc_vals <- nucs[ind]

    phi_prob <- sapply(1:K, function(k) {
      hap_probs <- sapply(1:length(ind), function(s) {
        phi[k,pos_vals[s],nuc_vals[s]]
      })
      return (prod(hap_probs))
    })

    probs <- pi*phi_prob

    return (probs/sum(probs))
  })

  return (marginal_prob)
}

#' Helper function to PredictHaplo that constructs read assignments
#'
#' Computes cvals given reads, pi, and phi
#' @param df A read table
#' @param pi Numeric vector of haplotype frequencies
#' @param phi A haplotype probability matrix from a PredictHaplo object
#'
#' @return A (K x M) matrix with K equaling the number of haplotypes and
#' M equaling the number of read groups (rows in df).
#'
#' @details Entry k,m of the returned matrix is the number of reads
#' in read group m that are assigned to haplotype k.
#' @export
assign_reads.PredictHaplo <- function(df, pi, phi)
{
  read_assignment_prob <-read_assignment_prob_matrix.PredictHaplo(df, pi, phi)

  # now assign the counts as multinomials
  # cvals[k,j] = number of reads from row j in df that are assigned to
  # haplotype k
  cvals <- sapply(1:nrow(df), function(i) {
    c_count <- df$count[i]
    c_prob <- read_assignment_prob[,i]

    rmultinom(1, c_count, c_prob)
  })

  return (cvals)
}

#' Calculate posterior probability of parameter values and assignments
#'
#' @param df Read tabel
#' @param pi,phi,cvals,alpha,gamma Parameter values from a PredictHaplo object
#'
#' @return The log-posterior probability
#' @export
posterior_log_prob.PredictHaplo <- function(df, pi, phi, cvals, alpha, gamma)
{
  K <- length(pi)
  read_assignment_prob <- read_assignment_prob_matrix.PredictHaplo(df, pi, phi)

  assign_prob <- sum(cvals*log(read_assignment_prob))

  pi_prior_prob <- (gamma/K-1)*sum(log(pi))

  phi_prior_prob <- 0
  npos <- ncol(df)-1
  for (k in 1:K) {
    phi_prior_prob_pos <- sapply(1:npos, function(i) {
      (alpha/4-1)*sum(log(phi[k,i,]))
    })
    phi_prior_prob <- phi_prior_prob + sum(phi_prior_prob_pos)
  }

  prob <- assign_prob + pi_prior_prob + phi_prior_prob

  return (prob)
}


#' Return the haplotypes corresponding to a PredictHaplo object
#'
#' @param ph PredictHaplo object
#'
#' @return A (K x P) matrix of A/C/G/T values with rows
#' corresponding to K haplotypes and columns corresponding to the number of
#' positions in the read table on which ph is based.
#'
#' @details The matrix phi[k,,] gives the nucledotide frequencies (cols) for each position
#' (rows).   The corresponding haplotype is formed by taking the consensus.
#' @export
haplotypes.PredictHaplo <- function(ph)
{
  phi <- ph$phi;
  K <- length(ph$pi)

  haps <- sapply(1:K, function(k) {
    m <- phi[k,,]
    ind <- apply(m, 1, which.max)
    colnames(m)[ind]
  })
  if (!is.matrix(haps))
    haps <- matrix(haps, ncol=1)

  rownames(haps) <- names(select(ph$df, -count))

  return (t(haps))
}

#' Measure fit of haplotypes and read assignments
#'
#' Given read assignments to a set of haplotypes, two
#' measures of fit are returned:  a mapping accuracy
#' (the fraction of reads that match their assigned haplotype)
#' and a mapping consistency (the variance of assignment frequencies
#' at each nt position)
#'
#' @param ph A PredictHaplo object
#'
fit.PredictHaplo <- function(ph)
{
  haplos <- haplotypes.PredictHaplo(ph)

  K <- length(ph$pi)
  pi <- ph$pi
  df_nucs <- select(ph$df, -count)
  nread_groups <- nrow(df_nucs)
  npos <- ncol(haplos)
  nhaplos <- nrow(haplos)
  cvals <- ph$cvals

  # for each entry in cval (haplo x read group) determine
  # if read and haplotype match
  cvals_same <- matrix(NA, nrow=nrow(cvals), ncol=ncol(cvals))
  for (k in 1:K)
    for (r in 1:nread_groups) {
      c_haplo <- haplos[k,]
      cread <- as.character(df_nucs[r,])
      cread_ind <- !is.na(cread)
      cvals_same[k,r] <- all(cread[cread_ind] == c_haplo[cread_ind])
    }

  # global accuracy
  mapping_accuracy <-  sum(cvals_same*cvals)/sum(cvals)
  # NEED TO GO THROUGH READ GROUPS AND CONNECT TO POSITION
  mapping_accuracy_v <- sapply(1:npos, function(i) {
    active_read_groups <- which(!is.na(df_nucs[,i]))
    if (length(active_read_groups)==0)
      return (NA)

    nmatch <- sum(cvals_same[,active_read_groups]*cvals[,active_read_groups])
    ntotal <- sum(cvals[,active_read_groups])

    nmatch/ntotal
  })

  # mapping consisteny
  consistency_info <- lapply(1:npos, function(i) {
    hap_counts <- rep(0, K)
    # which read groups cover position i
    active_read_groups <- which(!is.na(df_nucs[,i]))
    if (length(active_read_groups)==0)
      next

    for(arg in active_read_groups)
      hap_counts <- hap_counts + ph$cvals[,arg]

    hap_freqs <- hap_counts/sum(hap_counts)

    return (list(freq=hap_freqs, count=sum(hap_counts)))
  })

  mapping_consistency <- sapply(consistency_info, function(ci) ci$freq)
  if (!is.matrix(mapping_consistency))
    mapping_consistency <- matrix(mapping_consistency, nrow=1)

  consistency_v <- apply(mapping_consistency, 2, function(hat_pi) {
    sum(abs(hat_pi-pi))
  })

  counts <- sapply(consistency_info, function(ci) ci$count)

  colnames(mapping_consistency) <- names(df_nucs)
  names(counts) <- names(df_nucs)
  names(consistency_v) <- names(df_nucs)

  #consistency_scalar <- sum(apply(mapping_consistency, 1, var))

  return (list(accuracy=mapping_accuracy,
               accuracy_v=mapping_accuracy_v,
               consistency_m=mapping_consistency,
               consistency_v=consistency_v,
               counts=counts))

}
