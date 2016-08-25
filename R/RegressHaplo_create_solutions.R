#' Generate the y vector and P matrix of the RegressHaplo optimization and
#' the associated haplotyes
#'
#' Generates y and P using a call to filter_and_optimize.RegressHaplo and
#' also returns a matrix of consistent haploytpes
#'
#' @return a list with entries y (numeric), P(numeric matrix) and h
#' (a character matrix).
parameters.RegressHaplo <- function(df,
                                 max_global_dim=1200,
                                 max_local_dim=1200,
                                 min_cover=500)
{
  # get the read table
  rh <- filter_and_optimize.RegressHaplo(df,
                             global_rho=NULL,
                             max_global_dim=max_global_dim,
                             max_local_dim=max_local_dim,
                             min_cover=min_cover,
                             run_optimization = F)

  h <- get_h.RegressHaplo(rh)
  par <- penalized_regression_parameters.RegressHaplo(df, h)

  y <- matrix(par$y, ncol=1)
  P <- par$P

  return (list(y=y, P=P, h=h))
}

#' Solve the RegressHaplo optimization repeatedly
#'
#' For fixed y and P, generate solutions for every combination of
#' rho and kk in rho_vals and kk_vals.
#'
#' @return A matrix with each column corresponding to a solution.
#' Within each column, the first four rows are
#' fit, rho, K, and kk followed the computed solution.
#' @export
solutions.RegressHaplo <- function(y, P, num_trials,
                                   rho_vals=c(.1,1,2,5),
                                   kk_vals=2)
{
  K <- ncol(P)

  counter <- 1
  trial_vals <- 1:num_trials
  nsolutions <- length(kk_vals)*length(rho_vals)*num_trials

  cat("TOTAL SOLUTIONS TO BE GENERATED:", nsolutions, "\n")
  out_matrix <- matrix(NA, nrow=K+4, ncol=nsolutions)

  for (rho in rho_vals) {
    for (kk in kk_vals) {
      for (trial in trial_vals) {

        if (trial==1) {
          # here Igor finds the best fitting homogeneous solution
          # by trying all possible homogeneous solutions. Is this
          # worth the time?
          template <- rep(0, K)
          fits <- sapply(1:K, function(i) {
            template[i] <- 1
            res <- y - P %*% template
            fit <- sum(res*res)
            template[i] <- 0
            return (fit)
          })
          template[which.min(fits)] <- 1
          pi0 <- template
        } else if (trial==2) {
          pi0 <- rep(1/K, K)
        } else if ((trial>=3) && (trial<=10)) {
          pi0 <- runif(K)
        } else {
          pi0 <- runif(K)
          pi0 <- pi0/sum(pi0)
        }

        cat("RHO", rho, "kk", kk)
        cat("TRIAL", counter, "OF", nsolutions, "\n")

        reg <- penalized_regression.RegressHaplo(y, P, pi=pi0, rho=rho, kk=kk)
        # first four entries are fit, rho, K, kk
        cK <- sum(reg$pi > 0)
        out_matrix[,counter] <- c(reg$fit, rho, cK, kk, reg$pi)
        counter <- counter + 1
      } # end trial
    } # end kk
  } # end rho

  return (out_matrix)
}

