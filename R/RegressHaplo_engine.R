
grad_AL.engine <- function(y, P, M, rho, mu, kk, pi)
{
  all1 <- rep(1, length(pi))
  g <- -2*t(P)%*%(y - P %*% pi) + 2 * rho * M %*% pi - (mu - kk *(sum(pi)-1)) * all1
#  -2*Pmat'*(yvec-Pmat*y) + 2*rSL*qpM*y - mu*ydat + kk*(ydat'*y-1)*ydat;
  return (g)
}

AL.engine <- function(y, P, M, rho, mu, kk, pi)
{
  al1 <- y - P %*% pi
  sp <- sum(pi)
  al <- t(al1) %*% al1 + rho*t(pi) %*% M %*% pi - mu*(sp-1) + kk^2/2*(sp-1)^2

  return (al)
}


max_Hessian_eigenvalue.engine <- function(P, M, kk)
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



optimize.engine <- function(y, P, rho, pi, mu, kk)
{
  glaveps <- 1E-6

  colnames(P) <- 1:ncol(P)
  dim <- ncol(P)
  # create the M matrix
  M <- matrix(1, nrow=dim, ncol=dim) - diag(dim)

  # internal function called to reduce problem dimension as optimization iteratess
  remove_zero_coordinates <- function() {
    ind <- pi > 0
    P <- P[,ind]
    M <- M[,ind]
    M <- M[ind,]
    pi <- pi[ind]
    return (NULL)
  }

  # OUTER LOOP, RAISES PENALTY AND RECALCULATES LAGRANGE MULTIPLIES
  remove_zero_coordinates()
  grad <- grad_AL.engine(y, P, M, rho, mu, kk, pi)
  outer_loop_counter <- 0
  total_error <- abs(constraint_value.engine(pi)) +
                  gradient_norm.engine(pi, grad, glaveps)

  while(total_error > glaveps & outer_loop_counter < 1000) { # outerloop
    outer_loop_counter <- outer_loop_counter + 1

    t <- 1
    # checkstop

    # WHY IS L NOT ESTIMATED IN THE INTERNAL LOOP?
    L <- 1.25*max_Hessian_eigenvalue.engine(P, kk)

    # INNER LOOP THAT OPTIMIZES AUG LAGRANGIAN IN pi
    inner_loop_counter <- 0
    grad <- grad_AL.engine(y, P, M, rho, mu, kk, pi)
    repeat {
      inner_loop_counter <- inner_loop_counter + 1

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

  coor <- as.numeric(colnames(P))
  pi_out <- rep(0, dim)
  pi_out[coor] <- pi

  return (pi)
}
