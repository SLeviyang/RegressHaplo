
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

