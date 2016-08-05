
#' Note: exanding window approach did not work well.  Locally, we lack information so
#' that sparse solutions may be missed.   The best example is two biallelic unlinked
#' loci.  There is little information as to which solution is best.  A sparse solution
#' may select just two alleles.  But a three allele solution is also possible and fits
#' any marginals.  So better is to treat each locus separately and attempt a non-sparse
#' solution that allows for a modest dimension in the global problem.
#'
#'
RegressHaploFiltered <- function(df, global_rho, max_global_dim=500,
                                 max_local_dim=500,
                                 min_cover=0,
                                 run_optimization=T)
{
  # split the read table into unlinked read tables
  sdf <- split_read_table.RegressHaplo(df, min_cover=min_cover,
                                       max_dim=max_local_dim)

  #sdf <- split_unlinked_loci.read_table(df, sort_loci=F, min_cover=min_cover)


  # for each sdf, do a local pass and keep going until dimension is small
  # enough
  local_rho <- .001
  repeat {
    rh_local <- lapply(sdf, RegressHaploLocal,
                       max_dim=500,
                       min_cover=min_cover,
                       local_rho=local_rho)
    # extract all the haplotypes
    h_local <- lapply(rh_local, get_h.RegressHaplo)

    # form all permutations
    h_consistent <- haplotype_permute.RegressHaplo(h_local)

    if (nrow(h_consistent) <= max_global_dim)
      break
    else
      local_rho <- sqrt(10)*local_rho
  }

  # if not running optimization, then just return the consistent haplotypes
  if (!run_optimization)
    return (list(df=df, h=NA,
                 pi=NA,
                 fit=NA,
                 h_consistent=h_consistent))

  # compute_solution
  solution <- compute_solution.RegressHaplo(df, h_consistent, global_rho)

  return (list(df=df, h=solution$h,
               pi=solution$pi,
               fit=solution$fit,
               h_consistent=h_consistent))
}


#' Apply sparse haplotype algorithm to a read table
#'
#' @param df A read table
#' @param local_rho penalty coefficient
#'
#' @return A list with the following elements
#' \describe{
#' \item{df}{The read table passed}
#' \item{h}{A haplotype matrix computed by the algorithm}
#' \item{pi}{A numeric vector of frequencies for each haplotype (row) in h}
#' \item{fit}{Quality of fit of h,pi to the read table df}
#' \item{h_consistent}{A haplotype matrix of all consistent haplotypes for the read table}
#' }
RegressHaploLocal <- function(df, local_rho, max_dim=500,
                              min_cover=1000,
                              run_optimization=T)
{
  if (nrow(df)==0)
    return (list(df=df, h=NA, pi=NA, fit=NA, h_consistent=NA))

  h_consistent <- consistent_haplotypes_across_loci.read_table(df,
                                                               min_cover = min_cover)

  if (!run_optimization | nrow(h_consistent) > max_dim)
    return (list(df=df, h=NA,
                 pi=NA,
                 fit=NA,
                 h_consistent=h_consistent))

  solution <- compute_solution.RegressHaplo(df, h_consistent, local_rho)

  return (list(df=df, h=solution$h,
               pi=solution$pi,
               fit=solution$fit,
               h_consistent=h_consistent))
}


get_h_consistent.RegressHaplo <- function(rh)
{
  return (rh$h_consistent)
}

get_h.RegressHaplo <- function(rh)
{
  return (rh$h)
}

get_pi.RegressHaplo <- function(rh)
{
  return (rh$pi)
}

get_df.RegressHaplo <- function(rh)
{
  return  (rh$df)
}

get_fit.RegressHaplo <- function(rh)
{
  return (rh$fit)
}

split_read_table.RegressHaplo <- function(df, min_cover, max_dim)
{
  # split based on loci first
  sdf <- split_unlinked_loci.read_table(df, sort_loci=F, min_cover=min_cover)

  # cycle through until all read tables have less than max_dim consistent haplotypes
  # and cover less than 1000 bp
  sdf_pre <- sdf
  repeat {
    sdf_new <- list()
    for (i in 1:length(sdf_pre)) {
      pos <- as.numeric(pos_names.read_table(sdf_pre[[i]]))

      # check length of locus and if not too long, number of consistent haplotypes
      if (max(pos)-min(pos) < 1000)
        hc <- consistent_haplotypes.read_table(sdf_pre[[i]], rm.na=T)
      else # set hc so we split
        hc <- matrix(NA, nrow=max_dim+1, ncol=1)

      if (nrow(hc) <= max_dim)
        sdf_new <- append(sdf_new, sdf_pre[i])
      else {
        sdf_split <- split.read_table(sdf_pre[[i]])
        if (is.null(sdf_split))
          stop("BUG!  read table cannot be split to meet dimension requirements")

        sdf_new <- append(sdf_new, list(sdf_split$df1, sdf_split$df2))
      }
    }

    if (length(sdf_new)==length(sdf_pre))
      break
    else
      sdf_pre <- sdf_new

  }

  return (sdf_new)
}


compute_solution.RegressHaplo <- function(df, h_full, rho)
{
  par <- penalized_regression_parameters.RegressHaplo(df, h_full)
  pi_full <- penalized_regression.RegressHaplo(par$y, par$P, par$M,
                                               rho=rho,
                                               debug=F)

  # keep only significant frequencies
  ind <- which(pi_full > 10^-3)
  pi <- pi_full[ind]
  h <- h_full[ind,,drop=F]

  # once haplotypes are known, run regression with no penalty so
  # we are not biasing frequencies
  print("adjusting haplotype frequencies")
  par <- penalized_regression_parameters.RegressHaplo(df, h)
  if (length(pi) > 1)
    pi <- penalized_regression.RegressHaplo(par$y, par$P, par$M, rho=0)

  fit_vec <- par$y - par$P %*% pi
  fit <- sqrt(sum(fit_vec*fit_vec))
  rownames(h) <- round(100*pi)

  # sort h and pi so that freqs are decreasing
  ind <- order(pi, decreasing = T)
  h <- h[ind,,drop=F]
  pi <- pi[ind]

  return (list(h=h, pi=pi, fit=fit))
}


haplotype_permute.RegressHaplo <- function(h_list)
{
  pos_list <- lapply(h_list, colnames)
  pos <- do.call(c, pos_list)

  # for now remove NA
  h_list <- lapply(h_list, function(h) {
    ind <- apply(h, 1, function(hrow) all(!is.na(hrow)))
    matrix(as.character(h[ind,]), nrow=sum(ind))
  })
  h_list_s <- lapply(h_list, function(h) {
    apply(h, 1, paste, collapse="")
  })

  h_df_s <- expand.grid(h_list_s)
  haps <- apply(h_df_s, 1, function(hvec) {
    hvec_paste <- paste(hvec, collapse="")
    strsplit(hvec_paste, split="")[[1]]
  })
  if (class(haps) != "matrix")
    haps <- matrix(haps, ncol=1)
  else
    haps <- t(haps)

  colnames(haps) <- pos
  return (haps)
}
#
# EM.RegressHaplo <- function(df, h, max.iterations=100)
# {
#   K <- nrow(h)
#   pi <- runif(K)
#   pi <- pi/sum(pi)
#   phi <- phi.RegressHaplo(df, h)
#   # pi <- rep(1/K,K)
#
#   for (i in 1:max.iterations) {
#     #print(i)
#
#     assign <- read_assignment_prob_matrix.PredictHaplo(df,
#                                                        matrix(pi, nrow=1),
#                                                        phi)
#   #  assign <- read_assignment_prob_matrix.Shorah (df, pi, h)
#     count_m <- matrix(df$count, nrow=K,
#                       ncol=nrow(df),
#                       byrow=T)
#     assign_count <- assign*count_m
#
#     pi <- rowSums(assign_count/sum(df$count))
#   }
#
#   return (pi)
# }

#' Returns the matrices and vectors associated with the regression
#' optimization
#'
#' @details The regression solves
#' \min |y - P*pi|^2 + \rho \pi^T*M*\pi
#' subject to:  \sum \pi_i = 1, \pi_i \ge 0.
#'
#' @param df read table
#' @param h haplotype matrix
#'
#' @return the matrices P and M and the vector y as a list.
#' @export
penalized_regression_parameters.RegressHaplo <- function(df, h)
{
  # get position nucleotide counts
  #nucs_mat <- nucs_at_pos.read_table(df)
  K <- nrow(h)

  rf <- read_fit.readFit(df, h, pi=NULL)

  P_list <- lapply(rf, function(crf) crf$P)
  y_list <- lapply(rf, function(crf) crf$sampled_freq)

  P <- do.call(rbind, P_list)
  y <- do.call(c, y_list)

  # penalty matrix, it's concave!
  M <- matrix(1, nrow=K, ncol=K)*(1-diag(K))

  return (list(y=y, P=P, M=M))
}

#' Solves min_pi |y - P*pi|^2 + rho*pi^T*M*pi
#' @export
penalized_regression.RegressHaplo <- function(y, P, M, rho,
                                              debug=F)
{
  K <- nrow(M)

  pi0 <- runif(K)
  pi0 <- pi0/sum(pi0)

  quad_f_pi <- function(pi) {
    z <- y - P %*% pi
    # residual SS plus penalty for stability
    return (sum(z*z) + rho*t(pi) %*% M %*% pi)
  }

  quad_gf_pi <- function(pi) {
    grad_f <- 2*(t(P) %*% P + rho*M) %*% pi - 2*t(P) %*% y

    return (grad_f)
  }

  eq_constraint_pi <- function(pi) return(sum(pi)-1)

  solve <- solnp(pi0, quad_f_pi, eq_constraint_pi, 0,
                 LB=rep(0, K), UB=rep(1, K))

  pi <- solve$pars
  if (debug) browser()

  return (pi)
}

