
grad_AL.engine <- function(y, P, M, rho, mu, kk, pi)
{
  col1 <- matrix(1, nrow=ncol(P), ncol=1)
  QQ <- 2*(t(P) %*% P + rho*M) + kk*col1 %*% t(col1)
  Py_more <- -2*t(P) %*% y - (mu + kk)*col1

  ind <- which(pi > 0)
  gf <- QQ[,ind] %*% pi[ind] + Py_more

  return (gf)

  all1 <- rep(1, length(pi))
  g <- -2*t(P)%*%(y - P %*% pi) + 2 * rho * M %*% pi - (mu - kk *(sum(pi)-1)) * all1
#  -2*Pmat'*(yvec-Pmat*y) + 2*rSL*qpM*y - mu*ydat + kk*(ydat'*y-1)*ydat;
  return (g)
}


grad_AL_fast.engine <- function(QQ, Py_more, pi)
{
  ind <- which(pi > 0)
  return (QQ[,ind,drop=F] %*% pi[ind] + Py_more)
}

grad_AL_debug.engine <- function(y, P, M, rho, mu, kk, pi)
{
  col1 <- matrix(1, nrow=ncol(P), ncol=1)
  QQ <- 2*(t(P) %*% P + rho*M) + kk*col1 %*% t(col1)
  Py_more <- -2*t(P) %*% y - (mu + kk)*col1

  g <- grad_AL.engine(y, P, M, rho, mu, kk, pi)
  ind <- which(pi > 0)
  gf <- QQ[,ind] %*% pi[ind] + Py_more

#  if (any(abs(g-gf) > 1E-10))
#    browser()
#  else
#    return (gf)

  return (gf)
}

#' Compute the augmented lagrangian
#'
#' @details This function is not used in the optimization and is only for debugging
#' purposes
AL.engine <- function(y, P, M, rho, mu, kk, pi)
{
  al1 <- y - P %*% pi
  sp <- sum(pi)
  al <- t(al1) %*% al1 + rho*t(pi) %*% M %*% pi - mu*(sp-1) + kk^2/2*(sp-1)^2

  return (al)
}


max_Hessian_eigenvalue.engine <- function(P, M, rho, kk)
{
  nn <- nrow(M)
  all1 <- rep(1, nn)
  QQ = 2*(t(P) %*% P + rho*M) + kk*all1 %*% t(all1)

  q <- rep(1, nn)
  noq <- sqrt(nn)
  q <- q/noq;

  noqp <- 1e10;
  while (abs(noq-noqp)/noq > 0.1) {
    q = QQ %*% q;
    noqp <- noq;
    noq <- sqrt(sum(q*q))
    q <- q/noq
  }

  return (noq)
}

max_Hessian_eigenvalue_fast.engine <- function(QQ)
{
  nn <- ncol(QQ)

  q <- rep(1, nn)
  noq <- sqrt(nn)
  q <- q/noq;

  noqp <- 1e10;
  while (abs(noq-noqp)/noq > 0.1) {
    q = QQ %*% q;
    noqp <- noq;
    noq <- sqrt(sum(q*q))
    q <- q/noq
  }

  return (noq)
}

constraint_value.engine <- function(pi)
{
  return (sum(pi)-1)
}

gradient_norm.engine <- function(pi, grad, cutoff)
{
  boundary_coordinates <- (pi < cutoff)

  grad_boundary <- grad[boundary_coordinates]
  grad_other <- grad[!boundary_coordinates]

  grad_boundary_error <- ifelse(grad_boundary < 0, -grad_boundary, 0)
  grad_other_error <- abs(grad_other)

  error <- max(grad_boundary_error, grad_other_error)
  return (error)
}

debug_optimize.engine <- function(rho, mu, kk, fast=T)
{
  # save the initial pi
  P <- as.matrix(read.table("Igor_P1.csv", sep=","))
  y <- as.numeric(read.table("Igor_y1.csv")[,1])

  #set.seed(10)
  pi <- runif(ncol(P))
  write.table(pi, "debug_pi.csv", sep=",", row.names = F, col.names=F)

  if (fast)
    z <- optimize_fast.engine(y, P, rho, pi, mu, kk)
  else
    z <- optimize.engine(y, P, rho, pi, mu, kk)

  ind <- which(z > 0)

  zz <- matrix(z[ind], nrow=1)
  colnames(zz) <- ind

  if (fast)
    write.table(zz, "debug_out_fast.csv", row.names = F, sep=",")
  else
    write.table(zz, "debug_out.csv", row.names = F, sep=",")
  return (NULL)
}


optimize.engine <- function(y, P, rho, pi, mu, kk)
{
  glaveps <- 1E-6

  colnames(P) <- 1:ncol(P)
  dim <- ncol(P)
  # create the M matrix
  M <- matrix(1, nrow=dim, ncol=dim) - diag(dim)
  pi_prev <- pi

  # internal function called to reduce problem dimension as optimization iteratess
  remove_zero_coordinates <- function() {
    return (NULL)
    ind <- pi > 0

    P <<- P[,ind]
    M <<- M[,ind]
    M <<- M[ind,]
    pi <<- pi[ind]
    pi_prev <<- pi_prev[ind]
    return (NULL)
  }

  # OUTER LOOP, RAISES PENALTY AND RECALCULATES LAGRANGE MULTIPLIES
  remove_zero_coordinates()

  grad <- grad_AL.engine(y, P, M, rho, mu, kk, pi)
  outer_loop_counter <- 0
  total_error <- max(abs(constraint_value.engine(pi)),
                  gradient_norm.engine(pi, grad, glaveps))

  # for debug
  total_counter <- 1
  f <- c()
  c_err <- c()
  pi_dim <- c()

  # loops 1000 time!!!!!
  while(total_error > glaveps & outer_loop_counter < 1000) { # outerloop
    outer_loop_counter <- outer_loop_counter + 1

    remove_zero_coordinates()

    t <- 1
    # checkstop
    cat("OUTER LOOP", outer_loop_counter, pi[1], "\n")

    # WHY IS L NOT ESTIMATED IN THE INTERNAL LOOP?
    L <- 1.25*max_Hessian_eigenvalue.engine(P, M, rho, kk)

    # INNER LOOP THAT OPTIMIZES AUG LAGRANGIAN IN pi
    inner_loop_counter <- 0
    grad <- grad_AL.engine(y, P, M, rho, mu, kk, pi)

    repeat {
      inner_loop_counter <- inner_loop_counter + 1
      #debug
      f[total_counter] <- AL.engine(y, P, M, rho, mu, kk, pi)
      c_err[total_counter] <- constraint_value.engine(pi)
      pi_dim[total_counter] <- length(pi)
      total_counter <- total_counter + 1

      grad_error <- gradient_norm.engine(pi, grad, glaveps)
      if (outer_loop_counter==-1)
        cat("INNER LOOP", inner_loop_counter, pi[1], "\n")

      # end debug


#       print(f)
#       print(c_err)
#       browser()

      #  x = y - grad/L;
      pi_new <- pi - grad/L
      # x = max(x,0);
      pi_new <- ifelse(pi_new > 0, pi_new, 0)

      # nt = (1+sqrt(1+4*t^2))/2;
      nt <- (1 + sqrt(1 + 4*t^2))/2

      #y = x + (x - xp)*(t-1)/nt;
      pi <- pi_new + (pi_new - pi_prev)*(t-1)/nt
      #xp = x;
      pi_prev <- pi_new
      #t = nt;
      t <- nt


      remove_zero_coordinates()
      grad <- grad_AL.engine(y, P, M, rho, mu, kk, pi)
      grad_error <- gradient_norm.engine(pi, grad, glaveps)



      # SHOULD INNER LOOP COUNTER BE ALLOWED TO RISE TO 1E10??!!!!
      if (grad_error < .1*total_error | inner_loop_counter > 1E10)
        break

    }

    # augmented lagrangian multiplier update
    mu <- mu - kk*constraint_value.engine(pi);
    # penalty update
    kk <- kk*1.01;
    # update total error
    total_error <- abs(constraint_value.engine(pi)) + gradient_norm.engine(pi, grad, glaveps)
  }

  print(pi)

  browser()

  return (pi)
}

optimize_fast.engine <- function(y, P, rho, pi, mu, kk)
{
  glaveps <- 1E-6

  colnames(P) <- 1:ncol(P)
  dim <- ncol(P)
  # create the M matrix
  M <- matrix(1, nrow=dim, ncol=dim) - diag(dim)

  # if pi is not feasible stop!
  if (any(pi<0))
    stop("entries of pi must all be non-negative")
  pi_prev <- pi

  # matrices for quickly computing gradient and Hessian
  # that do not change during optimization
  col1 <- matrix(1, nrow=dim, ncol=1)
  PPandM <- 2*(t(P) %*% P + rho*M)
  col1col1 <- col1 %*% t(col1)
  Py <- -2*t(P) %*% y

  # matrices for quick computation that do change
  QQ <- PPandM + kk*col1col1
  Py_more <- Py - (mu + kk)*col1
  # gradient = QQ %*% pi - 2*t(P) %*% y - (mu + kk)col1 = QQ %*% pi + Py_more
  # Hessian = QQ


  # OUTER LOOP, RAISES PENALTY AND RECALCULATES LAGRANGE MULTIPLIES
  grad <- grad_AL_fast.engine(QQ, Py_more, pi)
  outer_loop_counter <- 0
  total_error <- max(abs(constraint_value.engine(pi)),
                     gradient_norm.engine(pi, grad, glaveps))

  # for debug
  total_counter <- 1
  f <- c()
  c_err <- c()
  pi_dim <- c()

  # loops 1000 time!!!!!
  while(total_error > glaveps & outer_loop_counter < 1000) { # outerloop
    outer_loop_counter <- outer_loop_counter + 1

    t <- 1
    # checkstop
    cat("OUTER LOOP", outer_loop_counter, "\n")

    # WHY IS L NOT ESTIMATED IN THE INTERNAL LOOP?
    L <- 1.25*max_Hessian_eigenvalue_fast.engine(QQ)

    # INNER LOOP THAT OPTIMIZES AUG LAGRANGIAN IN pi
    inner_loop_counter <- 0
    grad <- grad_AL_fast.engine(QQ, Py_more, pi)
    # debug
    #grad_debug <- grad_AL.engine(y, P, M, rho, mu, kk, pi)
    #grad_debug2 <- grad_AL_debug.engine(y, P, M, rho, mu, kk, pi)



    repeat {
      inner_loop_counter <- inner_loop_counter + 1
      #debug
      f[total_counter] <- AL.engine(y, P, M, rho, mu, kk, pi)
      c_err[total_counter] <- constraint_value.engine(pi)
      pi_dim[total_counter] <- length(pi)
      total_counter <- total_counter + 1

      grad_error <- gradient_norm.engine(pi, grad, glaveps)
      if (outer_loop_counter==3)
        cat("INNER LOOP", inner_loop_counter,
            grad_error, total_error, pi[1], "\n")
      # end debug


      #       print(f)
      #       print(c_err)
      #       browser()

      #  x = y - grad/L;
      pi_new <- pi - grad/L
      # x = max(x,0);
      pi_new <- ifelse(pi_new > 0, pi_new, 0)

      # nt = (1+sqrt(1+4*t^2))/2;
      nt <- (1 + sqrt(1 + 4*t^2))/2

      #y = x + (x - xp)*(t-1)/nt;
      pi <- pi_new + (pi_new - pi_prev)*(t-1)/nt
      #xp = x;
      pi_prev <- pi_new
      #t = nt;
      t <- nt


      grad <- grad_AL_fast.engine(QQ, Py_more, pi)
      # grad_debug <- grad_AL.engine(y, P, M, rho, mu, kk, pi)
      # if (any((grad-grad_debug)>1E-10)) browser()
      #  grad_debug2 <- grad_AL_debug.engine(y, P, M, rho, mu, kk, pi)


      grad_error <- gradient_norm.engine(pi, grad, glaveps)



      # SHOULD INNER LOOP COUNTER BE ALLOWED TO RISE TO 1E10??!!!!

      if (grad_error < .1*total_error | inner_loop_counter > 1E10)
        break

    }

    # augmented lagrangian multiplier update
    mu <- mu - kk*constraint_value.engine(pi);
    # penalty update
    kk <- kk*1.01;
    # update matrices for grad/Hessian
    QQ <- PPandM + kk*col1col1
    Py_more <- Py - (mu + kk)*col1
    # update total error
    total_error <- abs(constraint_value.engine(pi)) + gradient_norm.engine(pi, grad, glaveps)
  }

  return (pi)
}


