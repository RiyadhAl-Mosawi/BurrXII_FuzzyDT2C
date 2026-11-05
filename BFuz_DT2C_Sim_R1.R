########################################################################################################
# This Code was Written by:
#       Dr. Riyadh Rustam Al-Mosawi in 07-08-2025 
#------------------------------------------------------------------------------------
# Estimation the two parameters of Bur XII distribution 
# with double censored fuzzy observations
# We compute the MLE using NR and EM
# We compute the Bayesian estimators using Lindlye's, TK methods and MCMC


#--- The library of helper functions
##-------------------------------
# User-defined functions
##-------------------------------
setwd("C:/Users/hp/Dropbox/Riyadh_AlMosawi/CurrentFolder/SubmittedPapers/BurrFuzzyDT2Cens/Rcode")
source("Library.r")

# Install the package if needed, then load it
load_or_install = function(pkg) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  
  library(pkg, character.only = TRUE)
}



index = 7

settings = Select.Set(index)
list2env(settings, envir = .GlobalEnv)

n=settings$n
r=settings$r
s=settings$s
alp_true=settings$alp_true
bet_true=settings$bet_true
a=settings$a
b=settings$b
c=settings$c
d=settings$d

# Initial estimates
alp_ini = 0.1
bet_ini = 0.1

# Fuzzification
f1 = 0.05
f2 = 0.05

# EM settings
Min_Dif = 0.0001
Max_Iter = 30

# MCMC settings
N = 11000
P = 1000

# Monte Carlo EM settings
M = 100
K = 5

# Simulation replications
iterations = 2000

# Required packages for the revised functions
for (pkg in c("numDeriv", "coda", "knitr")) {
  load_or_install(pkg)
}

set.seed(2021)
start.time = Sys.time()

# Column order used by all methods
columns = c(
  "alp", "sd.alp", "lower.alp", "upper.alp", "coverage.alp",
  "bet", "sd.bet", "lower.bet", "upper.bet", "coverage.bet"
)

# Initialize results: one row for each replication
empty_results = function(mh = FALSE) {
  
  names = columns
  
  if (mh) {
    names = c(names, "acceptance")
  }
  
  matrix(
    NA_real_,
    nrow = iterations,
    ncol = length(names),
    dimnames = list(NULL, names)
  )
}

NR_MLE = empty_results()
EM_MLE = empty_results()
SM_MLE = empty_results()

TK_NON = empty_results()
TK_INF = empty_results()
LN_NON = empty_results()
LN_INF = empty_results()

MH_NON = empty_results(mh = TRUE)
MH_INF = empty_results(mh = TRUE)

# Catch errors without stopping the simulation
safe_run = function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      message(conditionMessage(e))
      NULL
    }
  )
}

# Check a method's error flag
successful = function(result, flag) {
  is.list(result) && identical(result[[flag]], FALSE)
}

# Extract estimates, SDs, intervals, and coverage
extract_results = function(result, prefix, bayes = FALSE) {
  
  estimate = if (bayes) "bay" else "mle"
  
  fields = paste0(
    prefix, ".",
    c(
      paste0(estimate, ".alp"),
      "ese.alp", "cl.alp", "cu.alp", "cp.alp",
      paste0(estimate, ".bet"),
      "ese.bet", "cl.bet", "cu.bet", "cp.bet"
    )
  )
  
  as.numeric(unlist(result[fields], use.names = FALSE))
}

# Extract TK or Lindley point estimates
extract_points = function(result, prefix) {
  
  values = rep(NA_real_, 10)
  
  values[1] = result[[paste0(prefix, ".bay.alp")]]
  values[6] = result[[paste0(prefix, ".bay.bet")]]
  
  values
}

for (i in seq_len(iterations)) {
  
  cat("** Index set : ",index," ** Replication: ", i, "of", iterations, "\n")
  
  Data = safe_run(Gen.Data(n,r,s,alp_true,bet_true))
  
  if (is.null(Data)) next
  
  E = Data$E
  
  if (length(E) == 0 ||
      any(!is.finite(E)) ||
      any(E <= 0) ||
      is.unsorted(E)) {
    message("Invalid generated sample.")
    next
  }
  
  # MLE using L-BFGS-B
  NR.Res = safe_run(NR.fun(alp_ini, bet_ini, E))
  
  if (successful(NR.Res, "Err_NR")) {
    NR_MLE[i, ] = extract_results(NR.Res, "nr")
  }
  
  # EM
  EM.Res = safe_run(EM.fun(alp_ini, bet_ini, E))
  
  if (successful(EM.Res, "Err_EM")) {
    EM_MLE[i, ] = extract_results(EM.Res, "em")
  }
  
  # Monte Carlo EM
  SM.Res = safe_run(SM.fun(alp_ini, bet_ini, E))
  
  if (successful(SM.Res, "Err.SM")) {
    SM_MLE[i, ] = extract_results(SM.Res, "sm")
  }
  
  # TK: non-informative prior
  TK.Res_non = safe_run(
    TK.fun(
      alp_ini, bet_ini, E,
      0.001, 0.001, 0.001, 0.001
    )
  )
  
  if (successful(TK.Res_non, "Err_TK")) {
    TK_NON[i, ] = extract_points(TK.Res_non, "tk")
  }
  
  # TK: informative prior
  TK.Res_inf = safe_run(
    TK.fun(alp_ini, bet_ini, E, a, b, c, d)
  )
  
  if (successful(TK.Res_inf, "Err_TK")) {
    TK_INF[i, ] = extract_points(TK.Res_inf, "tk")
  }
  
  # Lindley requires valid, positive MLEs.
  # Point estimates may exist even if the covariance failed.
  mle_available = (
    is.list(NR.Res) &&
      length(NR.Res$nr.mle.alp) == 1 &&
      length(NR.Res$nr.mle.bet) == 1 &&
      all(is.finite(c(NR.Res$nr.mle.alp, NR.Res$nr.mle.bet))) &&
      all(c(NR.Res$nr.mle.alp, NR.Res$nr.mle.bet) > 0)
  )
  
  if (mle_available) {
    
    LN.Res_non = safe_run(
      LN.fun(
        NR.Res$nr.mle.alp, NR.Res$nr.mle.bet, E,
        0.001, 0.001, 0.001, 0.001
      )
    )
    
    if (successful(LN.Res_non, "Err_LN")) {
      LN_NON[i, ] = extract_points(LN.Res_non, "ln")
    }
    
    LN.Res_inf = safe_run(
      LN.fun(
        NR.Res$nr.mle.alp, NR.Res$nr.mle.bet,
        E, a, b, c, d
      )
    )
    
    if (successful(LN.Res_inf, "Err_LN")) {
      LN_INF[i, ] = extract_points(LN.Res_inf, "ln")
    }
  }
  
  # MH requires the MLE covariance matrix as well
  if (successful(NR.Res, "Err_NR")) {
    
    MH.Res = safe_run(
      MH.fun(
        alp = NR.Res$nr.mle.alp,
        bet = NR.Res$nr.mle.bet,
        E = E,
        varcov = NR.Res$varcov,
        N = N, P = P,
        a = a, b = b, c = c, d = d,
        verbose = FALSE
      )
    )
    
    if (is.list(MH.Res)) {
      
      if (is.list(MH.Res$MH.Res_non)) {
        MH_NON[i, ] = c(
          extract_results(MH.Res$MH.Res_non, "mh", bayes = TRUE),
          MH.Res$MH.Res_non$mh.accp.rate
        )
      }
      
      if (is.list(MH.Res$MH.Res_inf)) {
        MH_INF[i, ] = c(
          extract_results(MH.Res$MH.Res_inf, "mh", bayes = TRUE),
          MH.Res$MH.Res_inf$mh.accp.rate
        )
      }
    }
  }
}

# Number of usable point estimates per method
results = list(
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

usable = sapply(results, function(x) {
  sum(is.finite(x[, 1]) & is.finite(x[, 6]))
})

print(data.frame(
  Method = names(usable),
  Usable = as.integer(usable),
  Failed_or_skipped = iterations - as.integer(usable)
))

# Display results if the function is available
if (exists("Pres.fun", mode = "function")) {
  Pres.fun(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF
  )
}

# Save results
if (exists("Save.Results", mode = "function")) {
  
  Save.Results(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF
  )
  
} else {
  
  save(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF,
    file = "Simulation_Results.RData"
  )
}

cat("Total execution time:\n")
print(Sys.time() - start.time)


Save.Results(
    NR_MLE, EM_MLE, SM_MLE,
    TK_NON, LN_NON, MH_NON,
    TK_INF, LN_INF, MH_INF)
