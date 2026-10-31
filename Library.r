# Probability density function
pdf = function(x, alp, bet) {
  alp * bet * x^(alp - 1) * (1 + x^alp)^(-bet - 1)
}
# Cumulative distribution function
cdf = function(x, alp, bet) {
  1 - (1 + x^alp)^(-bet)
}
# Inverse cumulative distribution function
qnt = function(p, alp, bet) {
  ((1 - p)^(-1 / bet) - 1)^(1 / alp)
}

# Triangular membership function
membership_fun = function(g) {
  
  memb.fun1 = function(w) {
    (w - (1 - f1) * g) / (f1 * g)
  }
  
  memb.fun2 = function(w) {
    ((1 + f2) * g - w) / (f2 * g)
  }
  
  list(
    memb.fun1 = memb.fun1,
    lower1 = (1 - f1) * g,
    upper1 = g,
    memb.fun2 = memb.fun2,
    lower2 = g,
    upper2 = (1 + f2) * g
  )
}

# Generate a doubly Type-II censored Burr XII sample
Gen.Data=function(n, r, s, alp, bet) {
  
  V = numeric(s)
  # Step 1
  V[s] = rbeta(1, shape1 = n - r, shape2 = r + 1)
  
  # Steps 2 and 3
  if (s > r + 1) {
    j = seq.int(r + 1, s - 1)
    Z = runif(length(j))
    V[j] = Z^(1 / (n - s - r + j))
  }
  
  # Step 4: successive products V_s, V_{s-1}V_s, ...
  index = seq.int(s, r + 1)
  U = -expm1(cumsum(log(V[index])))
  
  # Step 5
  E = qnt(U, alp, bet)
  data.frame(E = E)
}
## Fuzzy likelihood function

# Fuzzy log-likelihood
LogLike = function(y, E) {
  
  alp = y[1]
  bet = y[2]
  m = length(E)
  
  if (any(!is.finite(y)) || alp <= 0 || bet <= 0) {
    return(-Inf)
  }
  
  like = numeric(m)
  
  for (j in seq_len(m)) {
    
    membership = membership_fun(E[j])
    
    left = integrate(
      function(w) pdf(w, alp, bet) * membership[[1]](w),
      lower = membership[[2]],
      upper = membership[[3]]
    )$value
    
    right = integrate(
      function(w) pdf(w, alp, bet) * membership[[4]](w),
      lower = membership[[5]],
      upper = membership[[6]]
    )$value
    
    like[j] = left + right
  }
  
  if (any(!is.finite(like)) || any(like <= 0)) {
    return(-Inf)
  }
  
  # Contributions from retained fuzzy observations
  loglike = sum(log(like))
  
  # Left-censored contribution
  if (r > 0) {
    loglike = loglike + r * log(cdf(E[1], alp, bet))
  }
  
  # Right-censored contribution
  nr = n - r - m
  
  if (nr > 0) {
    loglike = loglike -
      bet * nr * log(1 + E[m]^alp)
  }
  
  return(loglike)
}

# Fuzzy likelihood
Like = function(y, E) {
  exp(LogLike(y, E))
}


# Expected complete-data log-likelihood (EM objective)
LogCompLike = function(y, E) {
  
  alp = y[1]
  bet = y[2]
  alp_old = y[3]
  bet_old = y[4]
  
  if (any(!is.finite(y)) || any(y <= 0)) {
    return(-Inf)
  }
  
  m = length(E)
  nr = n - r - m
  
  # Complete-data log-density at proposed parameters
  logpdf = function(w) {
    log(alp) + log(bet) +
      (alp - 1) * log(w) -
      (bet + 1) * log(1 + w^alp)
  }
  
  # Density at current EM estimates
  oldpdf = function(w) {
    pdf(w, alp_old, bet_old)
  }
  
  Q = 0
  
  # Retained fuzzy observations
  for (j in seq_len(m)) {
    
    membership = membership_fun(E[j])
    
    weight1 = function(w) {
      oldpdf(w) * membership[[1]](w)
    }
    
    weight2 = function(w) {
      oldpdf(w) * membership[[4]](w)
    }
    
    denominator = (
      integrate(weight1, membership[[2]], membership[[3]])$value +
        integrate(weight2, membership[[5]], membership[[6]])$value
    )
    
    numerator = (
      integrate(function(w) logpdf(w) * weight1(w),
                membership[[2]], membership[[3]])$value +
        integrate(function(w) logpdf(w) * weight2(w),
                  membership[[5]], membership[[6]])$value
    )
    
    Q = Q + numerator / denominator
  }
  
  # Left-censored observations
  if (r > 0) {
    
    numerator = integrate(
      function(w) logpdf(w) * oldpdf(w),
      lower = 0, upper = E[1]
    )$value
    
    denominator = cdf(E[1], alp_old, bet_old)
    
    Q = Q + r * numerator / denominator
  }
  
  # Right-censored observations
  if (nr > 0) {
    
    numerator = integrate(
      function(w) logpdf(w) * oldpdf(w),
      lower = E[m], upper = Inf
    )$value
    
    denominator = (1 + E[m]^alp_old)^(-bet_old)
    
    Q = Q + nr * numerator / denominator
  }
  
  return(Q)
}

CompLike = function(y, E) {
  exp(LogCompLike(y, E))
}

Prior1 = function(x, a, b) {
  dgamma(x, shape = a, rate = b)
}

Prior2 = function(x, c, d) {
  dgamma(x, shape = c, rate = d)
}

# Log-posterior, up to an additive constant
LogPostDist = function(x, E, a, b, c, d) {
  
  if (any(!is.finite(x)) || any(x <= 0)) {
    return(-Inf)
  }
  
  LogLike(x, E) +
    dgamma(x[1], shape = a, rate = b, log = TRUE) +
    dgamma(x[2], shape = c, rate = d, log = TRUE)
}

# Unnormalized posterior
PostDist = function(x, E, a, b, c, d) {
  exp(LogPostDist(x, E, a, b, c, d))
}

# Functions for the Tierney-Kadane method
H = function(y, E, a, b, c, d) {
  LogPostDist(y, E, a, b, c, d) / n
}

# For the posterior mean of alpha
Ha = function(y, E, a, b, c, d) {
  if (any(!is.finite(y)) || any(y <= 0)) return(-Inf)
  
  H(y, E, a, b, c, d) + log(y[1]) / n
}

# For the posterior mean of beta
Hb = function(y, E, a, b, c, d) {
  if (any(!is.finite(y)) || any(y <= 0)) return(-Inf)
  
  H(y, E, a, b, c, d) + log(y[2]) / n
}

# Gradient of the log-likelihood
dLogLike = function(y, E) {
  numDeriv::grad(function(w) LogLike(w, E), y)
}

# Hessian of the log-likelihood
ddLogLike = function(y, E) {
  numDeriv::hessian(function(w) LogLike(w, E), y)
}

# Log-prior density
Rho = function(y, a, b, c, d) {
  
  if (any(!is.finite(y)) || any(y <= 0)) {
    return(-Inf)
  }
  dgamma(y[1], shape = a, rate = b, log = TRUE) +
    dgamma(y[2], shape = c, rate = d, log = TRUE)
}


NR.fun = function(alp, bet, E) {
  
  # Default output if estimation fails
  result = list(
    Err_NR = TRUE,
    nr.mle.alp = NA_real_,
    nr.ese.alp = NA_real_,
    nr.cl.alp = NA_real_,
    nr.cu.alp = NA_real_,
    nr.cp.alp = NA_real_,
    nr.mle.bet = NA_real_,
    nr.ese.bet = NA_real_,
    nr.cl.bet = NA_real_,
    nr.cu.bet = NA_real_,
    nr.cp.bet = NA_real_,
    varcov = matrix(NA_real_, 2, 2)
  )
  
  # Maximize the log-likelihood
  fit = tryCatch(
    optim(
      par = c(alp, bet),
      fn = function(y) -LogLike(y, E),
      method = "L-BFGS-B",
      lower = c(1e-6, 1e-6),
      hessian = TRUE,
      control = list(maxit = 50000)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(result)
  if (fit$convergence != 0 || !is.finite(fit$value)) {
    return(result)
  }
  
  # Save point estimates
  estimate = fit$par
  result$nr.mle.alp = estimate[1]
  result$nr.mle.bet = estimate[2]
  
  # Check that the observed information is positive definite
  information = fit$hessian
  
  if (any(!is.finite(information))) return(result)
  
  check = tryCatch(
    chol(information),
    error = function(e) NULL
  )
  
  if (is.null(check)) return(result)
  
  # Covariance matrix and standard errors
  varcov = chol2inv(check)
  se = sqrt(diag(varcov))
  
  # 95% confidence intervals using the log transformation
  z = qnorm(0.975)
  lower = exp(log(estimate) - z * se / estimate)
  upper = exp(log(estimate) + z * se / estimate)
  
  # Coverage indicators
  truth = c(alp_true, bet_true)
  coverage = as.integer(lower <= truth & truth <= upper)
  
  result$Err_NR = FALSE
  
  result$nr.ese.alp = se[1]
  result$nr.cl.alp = lower[1]
  result$nr.cu.alp = upper[1]
  result$nr.cp.alp = coverage[1]
  
  result$nr.ese.bet = se[2]
  result$nr.cl.bet = lower[2]
  result$nr.cu.bet = upper[2]
  result$nr.cp.bet = coverage[2]
  
  result$varcov = varcov
  
  return(result)
}


EM.fun = function(alp = alp_ini, bet = bet_ini, E,
                   max_iter = Max_Iter, tol = Min_Dif) {
  
  # Default output
  result = list(
    Err_EM = TRUE,
    em.mle.alp = NA_real_,
    em.ese.alp = NA_real_,
    em.cl.alp = NA_real_,
    em.cu.alp = NA_real_,
    em.cp.alp = NA_real_,
    em.mle.bet = NA_real_,
    em.ese.bet = NA_real_,
    em.cl.bet = NA_real_,
    em.cu.bet = NA_real_,
    em.cp.bet = NA_real_
  )
  
  old = c(alp, bet)
  converged = FALSE
  
  if (any(!is.finite(old)) || any(old <= 0)) {
    return(result)
  }
  
  # EM iterations
  for (iter in seq_len(max_iter)) {
    
    # M-step: conditional distribution remains fixed at old
    fit = tryCatch(
      optim(
        par = old,
        fn = function(new) {
          -LogCompLike(c(new, old), E)
        },
        method = "L-BFGS-B",
        lower = c(1e-6, 1e-6),
        control = list(maxit = 50000)
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) return(result)
    
    if (fit$convergence != 0 ||
        !is.finite(fit$value) ||
        any(!is.finite(fit$par))) {
      return(result)
    }
    
    new = fit$par
    difference = sqrt(sum((new - old)^2))
    old = new
    
    if (difference < tol) {
      converged = TRUE
      break
    }
  }
  
  if (!converged) return(result)
  
  # Save converged point estimates
  estimate = old
  result$em.mle.alp = estimate[1]
  result$em.mle.bet = estimate[2]
  
  # Observed information from the observed-data likelihood
  information = tryCatch(
    numDeriv::hessian(
      func = function(y) -LogLike(y, E),
      x = estimate
    ),
    error = function(e) NULL
  )
  
  if (is.null(information)) return(result)
  if (any(!is.finite(information))) return(result)
  
  # Check positive definiteness
  check = tryCatch(
    chol(information),
    error = function(e) NULL
  )
  
  if (is.null(check)) return(result)
  
  # Standard errors
  varcov = chol2inv(check)
  se = sqrt(diag(varcov))
  
  # 95% confidence intervals using the log transformation
  z = qnorm(0.975)
  lower = exp(log(estimate) - z * se / estimate)
  upper = exp(log(estimate) + z * se / estimate)
  
  # Coverage indicators
  truth = c(alp_true, bet_true)
  coverage = as.integer(lower <= truth & truth <= upper)
  
  result$Err_EM = FALSE
  
  result$em.ese.alp = se[1]
  result$em.cl.alp = lower[1]
  result$em.cu.alp = upper[1]
  result$em.cp.alp = coverage[1]
  
  result$em.ese.bet = se[2]
  result$em.cl.bet = lower[2]
  result$em.cu.bet = upper[2]
  result$em.cp.bet = coverage[2]
  
  return(result)
}
SM.fun = function(alp, bet, E) {
  
  result = list(
    Err.SM = TRUE,
    sm.mle.alp = NA_real_,
    sm.ese.alp = NA_real_,
    sm.cl.alp = NA_real_,
    sm.cu.alp = NA_real_,
    sm.cp.alp = NA_real_,
    sm.mle.bet = NA_real_,
    sm.ese.bet = NA_real_,
    sm.cl.bet = NA_real_,
    sm.cu.bet = NA_real_,
    sm.cp.bet = NA_real_
  )
  
  m = length(E)
  nr = n - r - m
  
  stopifnot(
    m > 0, r >= 0, nr >= 0,
    M >= 1, K >= 1,
    alp > 0, bet > 0
  )
  
  old = c(alp, bet)
  best = NULL
  best_loglike = -Inf
  
  for (iter in seq_len(M)) {
    
    # Draw samples once for this iteration
    draws = weights = vector("list", m)
    
    for (j in seq_len(m)) {
      
      membership = membership_fun(E[j])
      
      lower = membership[[2]]
      center = E[j]
      upper = membership[[6]]
      
      p1 = cdf(lower, old[1], old[2])
      p2 = cdf(upper, old[1], old[2])
      
      draws[[j]] = qnt(
        runif(K, min = p1, max = p2),
        old[1], old[2]
      )
      
      weights[[j]] = ifelse(
        draws[[j]] <= center,
        membership[[1]](draws[[j]]),
        membership[[4]](draws[[j]])
      )
      
      if (any(!is.finite(draws[[j]])) ||
          any(!is.finite(weights[[j]])) ||
          sum(weights[[j]]) <= 0) {
        return(result)
      }
      
      weights[[j]] = weights[[j]] / sum(weights[[j]])
    }
    
    # Left-censored conditional draws
    left = numeric(0)
    
    if (r > 0) {
      left = qnt(
        runif(K, 0, cdf(E[1], old[1], old[2])),
        old[1], old[2]
      )
    }
    
    # Right-censored conditional draws
    right = numeric(0)
    
    if (nr > 0) {
      right = qnt(
        runif(K, cdf(E[m], old[1], old[2]), 1),
        old[1], old[2]
      )
    }
    
    all_draws = c(unlist(draws), left, right)
    
    if (any(!is.finite(all_draws)) || any(all_draws <= 0)) {
      return(result)
    }
    
    # Fixed Monte Carlo approximation to the EM objective
    Q = function(y) {
      
      if (any(!is.finite(y)) || any(y <= 0)) {
        return(-Inf)
      }
      
      logpdf = function(x) {
        log(y[1]) + log(y[2]) +
          (y[1] - 1) * log(x) -
          (y[2] + 1) * log(1 + x^y[1])
      }
      
      value = 0
      
      for (j in seq_len(m)) {
        value = value +
          sum(weights[[j]] * logpdf(draws[[j]]))
      }
      
      if (r > 0) {
        value = value + r * mean(logpdf(left))
      }
      
      if (nr > 0) {
        value = value + nr * mean(logpdf(right))
      }
      
      return(value)
    }
    
    fit = tryCatch(
      optim(
        par = old,
        fn = function(y) -Q(y),
        method = "L-BFGS-B",
        lower = c(1e-6, 1e-6),
        control = list(maxit = 5000)
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) return(result)
    
    if (fit$convergence != 0 || !is.finite(fit$value)) {
      return(result)
    }
    
    old = fit$par
    
    # Retain the iterate with the largest observed log-likelihood
    current_loglike = tryCatch(
      LogLike(old, E),
      error = function(e) -Inf
    )
    
    if (is.finite(current_loglike) &&
        current_loglike > best_loglike) {
      best = old
      best_loglike = current_loglike
    }
  }
  
  if (is.null(best)) return(result)
  
  result$sm.mle.alp = best[1]
  result$sm.mle.bet = best[2]
  
  # Observed information
  information = tryCatch(
    numDeriv::hessian(
      function(y) -LogLike(y, E),
      best
    ),
    error = function(e) NULL
  )
  
  if (is.null(information)) return(result)
  if (any(!is.finite(information))) return(result)
  
  check = tryCatch(
    chol(information),
    error = function(e) NULL
  )
  
  if (is.null(check)) return(result)
  
  se = sqrt(diag(chol2inv(check)))
  
  # Log-transformed 95% confidence intervals
  z = qnorm(0.975)
  lower = exp(log(best) - z * se / best)
  upper = exp(log(best) + z * se / best)
  
  truth = c(alp_true, bet_true)
  coverage = as.integer(lower <= truth & truth <= upper)
  
  result$Err.SM = FALSE
  
  result$sm.ese.alp = se[1]
  result$sm.cl.alp = lower[1]
  result$sm.cu.alp = upper[1]
  result$sm.cp.alp = coverage[1]
  
  result$sm.ese.bet = se[2]
  result$sm.cl.bet = lower[2]
  result$sm.cu.bet = upper[2]
  result$sm.cp.bet = coverage[2]
  
  return(result)
}

# Tierney-Kadane approximation of posterior means
TK.fun = function(alp = alp_ini, bet = bet_ini,
                   E, a, b, c, d) {
  
  result = list(
    Err_TK = TRUE,
    tk.bay.alp = NA_real_,
    tk.bay.bet = NA_real_
  )
  
  # Fit one of the three objectives
  fit_TK = function(fun, initial) {
    
    tryCatch(
      optim(
        par = initial,
        fn = function(y) -fun(y, E, a, b, c, d),
        method = "L-BFGS-B",
        lower = c(1e-6, 1e-6),
        hessian = TRUE,
        control = list(maxit = 50000)
      ),
      error = function(e) NULL
    )
  }
  
  # Check optimization results
  valid_fit = function(fit) {
    
    if (is.null(fit)) return(FALSE)
    
    fit$convergence == 0 &&
      is.finite(fit$value) &&
      all(is.finite(fit$par)) &&
      all(is.finite(fit$hessian)) &&
      all(fit$par > 1e-6)
  }
  
  # Posterior mode
  fit.H = fit_TK(H, c(alp, bet))
  if (!valid_fit(fit.H)) return(result)
  
  # Modes of the modified objectives
  fit.Ha = fit_TK(Ha, fit.H$par)
  fit.Hb = fit_TK(Hb, fit.H$par)
  
  if (!valid_fit(fit.Ha) || !valid_fit(fit.Hb)) {
    return(result)
  }
  
  # Log-determinant of a positive-definite matrix
  logdet = function(mat) {
    
    R = tryCatch(
      chol(mat),
      error = function(e) NULL
    )
    
    if (is.null(R)) return(NA_real_)
    
    2 * sum(log(diag(R)))
  }
  
  ld.H = logdet(fit.H$hessian)
  ld.Ha = logdet(fit.Ha$hessian)
  ld.Hb = logdet(fit.Hb$hessian)
  
  if (any(!is.finite(c(ld.H, ld.Ha, ld.Hb)))) {
    return(result)
  }
  
  # Calculate the approximations on the log scale
  log.alp = 0.5 * (ld.H - ld.Ha) +
    n * (fit.H$value - fit.Ha$value)
  
  log.bet = 0.5 * (ld.H - ld.Hb) +
    n * (fit.H$value - fit.Hb$value)
  
  estimates = exp(c(log.alp, log.bet))
  
  if (any(!is.finite(estimates)) || any(estimates <= 0)) {
    return(result)
  }
  
  result$Err_TK = FALSE
  result$tk.bay.alp = estimates[1]
  result$tk.bay.bet = estimates[2]
  
  return(result)
}

# Lindley approximation of posterior means
# alp and bet must be the observed-data MLEs
LN.fun = function(alp, bet, E, a, b, c, d) {
  
  result = list(
    Err_LN = TRUE,
    ln.bay.alp = NA_real_,
    ln.bay.bet = NA_real_
  )
  
  estimate = c(alp, bet)
  
  if (any(!is.finite(estimate)) || any(estimate <= 0)) {
    return(result)
  }
  
  # Observed-data log-likelihood
  loglike = function(y) {
    LogLike(y, E)
  }
  
  # Hessian of the log-likelihood
  Hess = function(y) {
    numDeriv::hessian(loglike, y)
  }
  
  derivatives = tryCatch({
    
    H = Hess(estimate)
    
    # Derivatives of the four Hessian entries
    third = numDeriv::jacobian(
      func = function(y) as.vector(Hess(y)),
      x = estimate
    )
    
    list(H = H, third = third)
    
  }, error = function(e) NULL)
  
  if (is.null(derivatives)) return(result)
  
  if (any(!is.finite(derivatives$H)) ||
      any(!is.finite(derivatives$third))) {
    return(result)
  }
  
  # Inverse observed information
  check = tryCatch(
    chol(-derivatives$H),
    error = function(e) NULL
  )
  
  if (is.null(check)) return(result)
  
  Sigma = chol2inv(check)
  
  # Third-derivative array:
  # T[i,j,k] = derivative of H[i,j] with respect to parameter k
  T = array(derivatives$third, dim = c(2, 2, 2))
  
  laaa = T[1, 1, 1]
  
  # Average equivalent mixed derivatives
  laab = mean(c(
    T[1, 1, 2],
    T[1, 2, 1],
    T[2, 1, 1]
  ))
  
  labb = mean(c(
    T[1, 2, 2],
    T[2, 1, 2],
    T[2, 2, 1]
  ))
  
  lbbb = T[2, 2, 2]
  
  # Gradient of the log-prior
  rho = c(
    (a - 1) / alp - b,
    (c - 1) / bet - d
  )
  
  saa = Sigma[1, 1]
  sab = Sigma[1, 2]
  sbb = Sigma[2, 2]
  
  # Third-derivative corrections
  v = c(
    laaa * saa + 2 * laab * sab + labb * sbb,
    laab * saa + 2 * labb * sab + lbbb * sbb
  )
  
  # Approximate posterior means
  bayes = as.vector(
    estimate + Sigma %*% rho + 0.5 * Sigma %*% v
  )
  
  if (any(!is.finite(bayes)) || any(bayes <= 0)) {
    return(result)
  }
  
  result$Err_LN = FALSE
  result$ln.bay.alp = bayes[1]
  result$ln.bay.bet = bayes[2]
  
  return(result)
}


MH.fun = function(alp, bet, E, varcov, N, P,
                   a, b, c, d, verbose = FALSE) {
  
  stopifnot(
    alp > 0, bet > 0,
    N == floor(N), P == floor(P),
    P >= 0, N - P >= 2,
    all(c(a, b, c, d) > 0)
  )
  
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Install the coda package first.")
  }
  
  # varcov is the covariance matrix on the ORIGINAL scale
  D = diag(1 / c(alp, bet))
  proposal_cov = D %*% varcov %*% D
  
  if (any(!is.finite(proposal_cov))) {
    stop("The proposal covariance must contain finite values.")
  }
  
  # Cholesky factor; also checks positive definiteness
  R = chol(proposal_cov)
  
  # Run one chain for a specified prior
  run_chain = function(prior) {
    
    current = log(c(alp, bet))
    
    # Posterior density on the log-parameter scale
    log_target = function(z) {
      
      theta = exp(z)
      
      if (any(!is.finite(theta)) || any(theta <= 0)) {
        return(-Inf)
      }
      
      value = LogPostDist(
        theta, E,
        prior[1], prior[2], prior[3], prior[4]
      ) + sum(z)
      
      if (!is.finite(value)) return(-Inf)
      
      return(value)
    }
    
    current_log = log_target(current)
    
    if (!is.finite(current_log)) {
      stop("The initial log-posterior is not finite.")
    }
    
    samples = matrix(NA_real_, nrow = N, ncol = 2)
    colnames(samples) = c("alpha", "beta")
    accepted = 0
    
    for (iter in seq_len(N)) {
      
      # Symmetric normal random walk on the log scale
      proposal = current +
        as.vector(t(R) %*% rnorm(2))
      
      proposal_log = log_target(proposal)
      log_ratio = proposal_log - current_log
      
      if (log(runif(1)) < min(0, log_ratio)) {
        current = proposal
        current_log = proposal_log
        accepted = accepted + 1
      }
      
      # Always record the current state
      samples[iter, ] = exp(current)
      
      if (verbose && iter %% 1000 == 0) {
        cat("Iteration:", iter, "\n")
      }
    }
    
    # Discard burn-in, including the case P = 0
    samples = samples[seq.int(P + 1, N), , drop = FALSE]
    
    means = colMeans(samples)
    posterior_sd = apply(samples, 2, sd)
    
    intervals = coda::HPDinterval(
      coda::mcmc(samples),
      prob = 0.95
    )
    
    lower = intervals[, "lower"]
    upper = intervals[, "upper"]
    
    truth = c(alp_true, bet_true)
    coverage = as.integer(lower <= truth & truth <= upper)
    
    list(
      mh.bay.alp = means[1],
      mh.ese.alp = posterior_sd[1],
      mh.cl.alp = lower[1],
      mh.cu.alp = upper[1],
      mh.cp.alp = coverage[1],
      
      mh.bay.bet = means[2],
      mh.ese.bet = posterior_sd[2],
      mh.cl.bet = lower[2],
      mh.cu.bet = upper[2],
      mh.cp.bet = coverage[2],
      
      mh.accp.rate = accepted / N,
      MH_sample = samples
    )
  }
  
  list(
    MH.Res_non = run_chain(c(0.001, 0.001, 0.001, 0.001)),
    MH.Res_inf = run_chain(c(a, b, c, d))
  )
}

# Safely calculate a mean
safe_mean = function(x) {
  
  x = x[is.finite(x)]
  
  if (length(x) == 0) return(NA_real_)
  
  mean(x)
}


# Summarize one parameter for one method
parameter_summary = function(x, truth, parameter) {
  
  # alpha occupies columns 1:5; beta occupies columns 6:10
  j = if (parameter == "alpha") 1 else 6
  
  estimate = x[, j]
  lower = x[, j + 2]
  upper = x[, j + 3]
  
  valid_estimate = is.finite(estimate)
  
  valid_interval = (
    is.finite(lower) &
      is.finite(upper) &
      lower <= upper
  )
  
  errors = estimate[valid_estimate] - truth
  
  # Recalculate coverage from the interval endpoints
  coverage = (
    lower[valid_interval] <= truth &
      truth <= upper[valid_interval]
  )
  
  data.frame(
    Parameter = parameter,
    N_est = sum(valid_estimate),
    N_interval = sum(valid_interval),
    Estimate = safe_mean(estimate),
    Bias = safe_mean(errors),
    AEB = abs(safe_mean(errors)),
    MSE = safe_mean(errors^2),
    RMSE = sqrt(safe_mean(errors^2)),
    Mean_lower = safe_mean(lower[valid_interval]),
    Mean_upper = safe_mean(upper[valid_interval]),
    Mean_length = safe_mean(
      upper[valid_interval] - lower[valid_interval]
    ),
    ECP = safe_mean(as.numeric(coverage))
  )
}


# Construct a common summary table
Summary.Results = function(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF) {
  
  methods = list(
    MLE = NR_MLE,
    EM = EM_MLE,
    MCEM = SM_MLE,
    LN.NON = LN_NON,
    TK.NON = TK_NON,
    MH.NON = MH_NON,
    LN.INF = LN_INF,
    TK.INF = TK_INF,
    MH.INF = MH_INF
  )
  
  output = vector("list", length(methods))
  
  for (j in seq_along(methods)) {
    
    x = methods[[j]]
    
    alpha = parameter_summary(x, alp_true, "alpha")
    beta = parameter_summary(x, bet_true, "beta")
    
    output[[j]] = data.frame(
      Method = names(methods)[j],
      rbind(alpha, beta),
      row.names = NULL
    )
  }
  
  result = do.call(rbind, output)
  rownames(result) = NULL
  
  return(result)
}


# Display simulation summaries
Pres.fun = function(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF) {
  
  result = Summary.Results(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF
  )
  
  display = result[, c(
    "Method", "Parameter", "N_est", "N_interval",
    "Estimate", "Bias", "MSE", "ECP"
  )]
  
  # Display coverage as a percentage
  display$ECP = 100 * display$ECP
  
  for (name in c("Estimate", "Bias", "MSE", "ECP")) {
    display[[name]] = round(display[[name]], 3)
  }
  
  print(display, row.names = FALSE)
  
  invisible(result)
}


# Save raw results and summary tables
Save.Results = function(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF) {
  
  result = Summary.Results(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF
  )
  
  raw_results = list(
    NR_MLE = NR_MLE,
    EM_MLE = EM_MLE,
    SM_MLE = SM_MLE,
    TK_NON = TK_NON,
    LN_NON = LN_NON,
    MH_NON = MH_NON,
    TK_INF = TK_INF,
    LN_INF = LN_INF,
    MH_INF = MH_INF
  )
  
  settings = list(
    n = n, r = r, s = s,
    alp_true = alp_true,
    bet_true = bet_true
  )
  
  filename = paste0("Simulation_n", n, "_r", r, "_s", s)
  
  # Save full-precision summaries and raw simulation results
  save(
    raw_results, result, settings,
    file = paste0(filename, ".RData")
  )
  
  write.csv(
    result,
    file = paste0(filename, ".csv"),
    row.names = FALSE,
    na = "NA"
  )
  
  # Optional LaTeX table
  if (requireNamespace("knitr", quietly = TRUE)) {
    
    latex = knitr::kable(
      result,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      row.names = FALSE
    )
    
    writeLines(
      as.character(latex),
      con = paste0(filename, ".tex")
    )
  }
  
  invisible(result)
}

# parameters
Select.Set = function(index) {
  
  if (length(index) != 1 || !is.finite(index) ||
      index != floor(index) || index < 1 || index > 18) {
    stop("index must be an integer from 1 to 18.")
  }
  
  # Nine censoring schemes
  schemes = data.frame(
    n = c(25, 25, 25, 50, 50, 50, 100, 100, 100),
    r = c(0, 5, 10, 0, 10, 20, 0, 20, 40),
    s = c(15, 15, 15, 30, 30, 30, 60, 60, 60)
  )
  
  j = (index - 1) %% 9 + 1
  
  settings = as.list(schemes[j, ])
  
  if (index <= 9) {
    parameters = list(
      alp_true = 0.5, bet_true = 1.5,
      a = 1, b = 2, c = 3, d = 2
    )
  } else {
    parameters = list(
      alp_true = 1.5, bet_true = 1,
      a = 3, b = 2, c = 3, d = 3
    )
  }
  
  return(append(settings, parameters))
}

Save.Results = function(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF) {
  
  # Summaries for the current censoring scheme
  current = Summary.Results(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF
  )
  
  current$r = r
  current$s = s
  
  # Separate files for each sample size and parameter combination
  filename = paste0(
    "Results_n", n,
    "_alpha", alp_true,
    "_beta", bet_true
  )
  
  summary_file = paste0(filename, ".rds")
  
  # Add the current scheme to previously saved schemes
  if (file.exists(summary_file)) {
    
    previous = readRDS(summary_file)
    
    # Replace an existing result for the same scheme
    previous = previous[
      !(previous$r == r & previous$s == s), ,
      drop = FALSE
    ]
    
    results = rbind(previous, current)
    
  } else {
    results = current
  }
  
  saveRDS(results, summary_file)
  
  # Save raw results for the current scheme
  save(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF,
    file = paste0(filename, "_r", r, "_s", s, ".RData")
  )
  
  # Format numbers; missing values produce blank cells
  number = function(x) {
    if (length(x) != 1 || !is.finite(x)) return("")
    sprintf("%.3f", x)
  }
  
  # Five columns for one parameter
  cells = function(row) {
    
    interval = ""
    
    if (is.finite(row$Mean_lower) &&
        is.finite(row$Mean_upper)) {
      
      interval = paste0(
        "(", number(row$Mean_lower),
        ", ", number(row$Mean_upper), ")"
      )
    }
    
    c(
      number(row$AEB),
      number(row$RMSE),
      interval,
      number(row$Mean_length),
      number(row$ECP)
    )
  }
  
  # Method order and printed labels
  methods = c(
    "MLE", "EM", "MCEM",
    "LN.NON", "TK.NON", "MH.NON",
    "LN.INF", "TK.INF", "MH.INF"
  )
  
  labels = c(
    "MLE", "EM", "MCEM",
    "LN.NON", "TK.NON", "MH.NON",
    "LN.INF", "TK.INF", "MH.INF"
  )
  
  # Table opening
  latex = c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0(
      "\\caption{Simulation results for classical and Bayesian ",
      "estimation for $n=", n, "$.}"
    ),
    paste0(
      "\\label{tab:sim_n", n,
      "_a", alp_true, "_b", bet_true, "}"
    ),
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{rlrrlrrrrlrr}",
    paste0(
      "&&\\multicolumn{5}{c}{$\\alpha=", alp_true, "$}",
      "&\\multicolumn{5}{c}{$\\beta=", bet_true, "$} \\\\"
    ),
    "\\hline",
    paste0(
      "$(r,s)$ & Method & AEB & RMSE & ACI/HPD & LEN & ECP",
      " & AEB & RMSE & ACI/HPD & LEN & ECP \\\\"
    ),
    "\\hline"
  )
  
  # Available censoring schemes
  schemes = unique(results[, c("r", "s")])
  schemes = schemes[order(schemes$r, schemes$s), , drop = FALSE]
  
  for (i in seq_len(nrow(schemes))) {
    
    block = results[
      results$r == schemes$r[i] &
        results$s == schemes$s[i], ,
      drop = FALSE
    ]
    
    for (j in seq_along(methods)) {
      
      alpha = block[
        block$Method == methods[j] &
          block$Parameter == "alpha", ,
        drop = FALSE
      ]
      
      beta = block[
        block$Method == methods[j] &
          block$Parameter == "beta", ,
        drop = FALSE
      ]
      
      scheme_label = ""
      
      if (j == 1) {
        scheme_label = paste0(
          "\\multirow{9}{*}{(",
          schemes$r[i], ",", schemes$s[i], ")}"
        )
      }
      
      row = c(
        scheme_label,
        labels[j],
        cells(alpha),
        cells(beta)
      )
      
      latex = c(
        latex,
        paste0(paste(row, collapse = " & "), " \\\\")
      )
    }
    
    latex = c(latex, "\\hline")
  }
  
  # Table closing
  latex = c(
    latex,
    "\\end{tabular}%",
    "}",
    "\\end{table}"
  )
  
  writeLines(latex, paste0(filename, ".tex"))
  
  message("LaTeX table saved: ", filename, ".tex")
  
  invisible(results)
}
