#' Loads a set of solutions created by solutions.RegressHaplo
#'
#' @param solutions Output from solutions.RegressHaplo
#' @param h Haplotypes used to construct the y and P matrix used
#' in generating solutions
#' @param Igor For now Igor's format and mine are slightly different!
#'
#' @return A list containing the entries df_stats, df_p,
#' read_df, founder_df, h_full, and results.
#' df_stats is a data.frame
#' containing the columns run_id,
#' rho, K, Igor_K, fit.  m_p is a matrix with rows p1, ..., pn;  see
#' \code{\link{parse_Igor_csv_to_df.IgorResults}} for further
#' details.  read_df and founder_df are the read_table for the
#' dataset and founder table for the animal, respectively.
#' h_full is the collection of hapltypes corresponding to the
#' entries of p.  results is a list, with an entry for each row
#' of results_df, i.e. each optimization.  Each entry of results
#' is a list containing h (haplotype matrix), pi (haplotype frequencies),
#' and fit.
RegressHaploSolutions <- function(solutions, h_full, Igor=F)
{
  results_dir <- get_simulated_read_data_dir(sim_num)

  df_list <- parse_solutions.RegressHaploSolutions(solutions)
  df_stats <- df_list$df_stats
  m_p <- as.matrix(df_list$df_p)
  colnames(m_p) <- NULL
  rownames(m_p) <- paste("solution", 1:nrow(m_p), sep="")

  h_full <- read.table(data_files[["h"]], sep=",", header=T, check.names = F,
                       colClasses = "character")
  h_full <- as.matrix(h_full)


  return (list(df_stats=df_stats,
               m_p=m_p,
               h=h_full))
}



#' Parse solutions matrix
#'
#' @details A solutions matrix is a matrix with each column
#' corresponding to a solution.
#' Within each column, the first four rows are
#' fit, rho, K, and kk followed the computed solution.
#'
#' @param filename CSV filename
#' @param animal animal to which this CSV file is associated
#' @param wekk week to which this CSV file is associated
#'
#' @return A data.frame with the first four columns, rho, K, K_Igor, fit,
#' followed by columns of dimension p, giving the optimal result.
#' K and K_Igor will be different because Igor counts an entry > 0 while
#' K counts > .01.
parse_solutions.RegressHaploSolutions <- function(solutions)
{
  rhos <- as.numeric(solutions[2,])
  K <- as.numeric(solutions[3,])
  fits <- as.numeric(solutions[1,])
  kk <- as.numeric(solutions[4,])

  p_matrix <- t(as.matrix(solutions[5:nrow(solutions),]))
  K <- apply(p_matrix, 1, function(p) sum(p > 0))

  df_out <- data.frame(rho=rhos, K=K, kk=kk, fit=fits, p_matrix)

  names(df_out) <- c("rho", "K", "kk", "fit",
                     paste("p_coor", 1:ncol(p_matrix), sep=""))

  df_stats <- df_out[,1:4]
  df_stats <- mutate(df_stats, solution_number=1:nrow(df_stats))
  df_p <- df_out[,5:ncol(df_out)]

  return (list(df_stats=df_stats, df_p=df_p))

}

######################################################
get_K.RegressHaploSolutions <- function(ir)
{
  return (unique(ir$df_stats$K))
}


get_stats_df.RegressHaploSolutions <- function(ir)
{
  return (ir$df_stats)
}


get_h.RegressHaploSolutions <- function(ir)
{
  return (ir$h)
}

get_solutions.RegressHaploSolutions <- function(ir)
{
  return (ir$m_p)
}

##########################################################
plot_fits_vs_K.RegressHaploSolutions <- function(sr, minK=NULL)
{
  df_stat <- get_stats_df.RegressHaploSolutions(sr)
  if (!is.null(minK))
    df_stat <- filter(df_stat, K >= minK)
  plot(df_stat$K, df_stat$fit, xlab="K", ylab="fit")
  return (NULL)
}

plot_fits_vs_rho.RegressHaploSolutions <- function(sr)
{
  df_stat <- get_stats_df.RegressHaploSolutions(sr)
  plot(df_stat$rho, df_stat$fit, xlab="rho", ylab="fit")
  return (NULL)
}

#' Determines the best fit solution
#'
#' @param rhs A RegressHaploSolutions object
#'
#' @details K is given by the number of non-zero entries in pi returned
#' by RegressHaplo
best_K.RegressHaploSolutions <- function(rhs)
{
  df_stats <- get_stats_df.RegressHaploSolutions(rhs)

  df <- select(df_stats, fit, K)
  df_K <- ddply(df, .(K), function(cdf) {
    data.frame(fit=min(cdf$fit), K=cdf$K[1])
  })

  # just in case sort
  ind <- order(df_K$K)
  df_K <- df_K[ind,]

  nk <- nrow(df_K)
  improvement <- c(NA, df_K$fit[-nk]/df_K$fit[-1])
  df_K <- mutate(df_K, improvement=improvement)

  # find the last index at which improvement > 1.1
  improved_ind <- which(improvement > 1.1)
  if (length(improved_ind)==0)
    best_ind <- 1
  else
    best_ind <- max(improved_ind)

  best_K <- df_K$K[best_ind]
  return (best_K)
}

best_fit.RegressHaploSolutions <- function(ir, K_val=NULL)
{
  df_stats <- get_stats_df.RegressHaploSolutions(ir)
  if (is.null(K_val))
    K_val <- best_K.RegressHaploSolutions(ir)

  if (!is.element(K_val, get_K.SimultionResults(ir)))
    stop("K_val is not a K value in the solutions")

  df_stats <- filter(df_stats, K==K_val)

  ind <- which.min(df_stats$fit)
  optim_solution_num <- df_stats$solution_number[ind]

  p <- get_solutions.RegressHaploSolutions(ir)[optim_solution_num,]
  h_full <- get_h.RegressHaploSolutions(ir)
  h_ind <- which(p > 0)

  pi <- p[h_ind]
  haplo <- h_full[h_ind,,drop=F]

  H <- Haplo(haplo, pi)

  fit <- df_stats$fit[optim_solution_num]

  return (list(H=H, fit=fit, solution=optim_solution_num))
}

# #' type is "map", "fit_position", "fit_read"
# plot_best_fit.IgorResults <- function(ir, K_val,
#                                       type="map",
#                                       plotToFile=F)
# {
#   optim <- best_fit.IgorResults(ir, K_val)
#   h <- optim$h
#   pi <- optim$pi
#   fit <- optim$fit
#
#   df <- get_read_df.IgorResults(ir)
#   paired_end <- get_paired_end.IgorResults(ir)
#   animal <- get_animal.IgorResults(ir)
#   week <- get_week.IgorResults(ir)
#
#   if (type=="map") {
#     rm <- readMap(df, h, pi)
#     p <- plot.readMap(rm, df=df, animal=animal, paired_end=paired_end)
#   }
#   else if (type == "fit_position") {
#     rf_position <- readFit(df, h, pi, position=T)
#     p <- plot.readFit(rf_position)
#   }
#   else {
#     rf_read <- readFit(df, h, pi, position=F)
#     p <- plot.readFit(rf_read)
#   }
#
#   p <- p + ggtitle(paste(animal, week, fit, sep=" "))
#
#   results_dir <- get_results_dir.IgorResults(animal, week)
#   outfile <- paste(results_dir, "best_fit_K=", K_val,
#                    "_", type, ".jpeg", sep="")
#
#   if (plotToFile) {
#     jpeg(outfile, width=1500, height=1000)
#     print(p)
#     dev.off()
#
#     return (fit)
#   }
#   else
#     return (p)
# }
#
