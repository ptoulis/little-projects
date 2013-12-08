## Stat 240 Homework 2. Analysis of Hirano et. al.
##
## Panos Toulis, panos.toulis@gmail.com
source("terminology.R")
library(coda)

get.compliance.probs <- function(theta, X) {
  #  Computes Ψ = (ψ_it) = P(Ci=t | Xi, θ)
  # 
  # Returns: (N x p) matrix of probabilities.
  #
  N = nrow(X)
  p = ncol(X)
  Pred = X %*% t(get.psi.matrix(theta))
  # if large/small values become a problem uncomment these.
  #Pred[Pred < -100] = -100
  #Pred[Pred > 100] = 100
  Psi = exp(Pred)  # (N x 3) matrix
  r = rowSums(Psi)
  Psi = Psi / r   # now Psi(i, t) = P(Ci=t | Xi, θ)
  CHECK_TRUE(all(!is.na(Psi)), msg="NA values in the Psi matrix..")
  return(Psi)
}

conditional.Ptyz <- function(t, z, Y, X, theta) {
  # Computes conditionl distribution over all units of 
  # Y | C=t, ...
  #
  # Args:
  #   t = compliance type
  #   z =  {0, 1}, treatment assignment
  #   Y = vector of outcomes
  #   X = covariates  (N x p) matrix
  #   Ψ = N x 3  matrix of category probabilities (Ψij = P(Ci = j | X, θ))
  #
  # Returns:
  #   V = (v_i)  where v_i = P(Yi=yi, Zi=z | Ci=t, Xi, θ)
  CHECK_EQ(nrow(X), length(Y))
  CHECK_MEMBER(z, c(0,1), msg="z in {0,1}")
  
  beta = get.beta.vector(t, z, theta)  # 3x1 vector
  preds = X %*% beta
  L1 = logistic(preds)
  # 1. Get predictive probabilities
  Ltz = Y * L1  + (1-Y) * (1-L1)  # (N X 1) where Ltz_i = P(Yi=yi | Ci=t, Zi=z, Xi, θ)
  
  return(Ltz) # ( N x 1) vector s.t. v(i) = P(Yi=1 | Ci=t, Zi=t, Xi, θ) 
}

joint.Ptyz <- function(t, z, Y, X, theta) {
  # Computes joint distribution over all units of 
  # C=t (compliance) and Y (outcomes) given the specified *individual* assignment
  #
  # Args:
  #   t = compliance type
  #   z =  {0, 1}, treatment assignment
  #   Y = vector of outcomes
  #   X = covariates  (N x p) matrix
  #   Ψ = N x 3  matrix of category probabilities (Ψij = P(Ci = j | X, θ))
  #
  # Returns:
  #   V = (v_i)  where v_i = P(Ci=t, Zi=z, Yi=yi | Xi, θ)
  #
  # 1. Get conditional/predictive probabilities
  Ltz = conditional.Ptyz(t, z, Y, X, theta)
  # 2. Get compliance "prior" probabilities.
  Psi = get.compliance.probs(theta, X)
  
  CHECK_EQ(length(Ltz), length(Psi[, t]))
  return(Ltz * Psi[, t]) # ( N x 1) vector s.t. v(i) = P(Yi=1 | Ci=t, Zi=t, Xi, θ) * P(Ci=t|...)
}

get.stratified.probs <- function(Data, theta) {
  # Stratifies by compliance type and then computes data probabilities.
  #
  # Returns: list(Szw, Qtz)
  #
  # Sij = (N x 1) vector where Szw_i = 1 if Zi=z and Wi=1
  # So, S01 are those Zi=0 and Wi=1, i.e. always-takers
  #     S00 are those Zi=0 and Wi=0, so either compliers or never-takers.
  #
  # Qtz = (N x 1) vector with Qtz_i = P(Ci=t, Zi=z, Y=Yobs_i | Χ_i, θ)
  #
  CHECK_Data(Data)
  CHECK_Theta(theta)
  
  Qtz <- function(t, Zi) {
    joint.Ptyz(t, z=Zi, Data$Y, Data$X, theta)
  }
  # Qtz = P(Y, C=t, Z=z | X, θ) = (Nx1) matrix
  # i.e. Qtz_i = P(Yi=ti, Ci=t, Zi=z | Xi, θ)
  Qn0 = Qtz(t=kTYPE_NEVERTAKER, Zi=0)
  Qn1 = Qtz(t=kTYPE_NEVERTAKER, Zi=1)
  
  Qc0 = Qtz(t=kTYPE_COMPLIER, Zi=0)
  Qc1 = Qtz(t=kTYPE_COMPLIER, Zi=1)
  
  Qa0 = Qtz(t=kTYPE_ALWAYSTAKER, Zi=0)
  Qa1 = Qtz(t=kTYPE_ALWAYSTAKER, Zi=1)
  
  W = Data$W
  Z = Data$Z
  
  S00 = (1-Z) * (1-W)
  S10 = Z * (1-W)
  S01 = (1-Z) * W
  S11 = Z * W
  
  CHECK_SETEQ(S00 + S01 + S11 + S10, c(1), msg="disjoint")
  return(list(S00=S00, S10=S10, S01=S01, S11=S11,
              Qn0=Qn0, Qn1=Qn1,
              Qc0=Qc0, Qc1=Qc1,
              Qa0=Qa0, Qa1=Qa1))
}

log.likelihood <- function(theta, Data) {
  # Computes observed log-likelihood
  # Tries to linear-algebrize computations.
  # About 300x  faster than the 'slow' version which iterates over units.
  #
  Cprobs <- get.stratified.probs(Data, theta)
  # Qtz = P(data | θ, C=t, Z=z)
  ll = with(Cprobs, {sum(S00 * log(Qc0 + Qn0) + S10 * log(Qn1) +
                         S01 * log(Qa0) + S11 * log(Qc1 + Qa1))
  })
  return(ll)
}

log.likelihood.complete <- function(theta, compliance, Data) {
  # See the difference with obs log likelihood:
  #
  #
  CHECK_Theta(theta)
  CHECK_EQ(length(Data$Y), length(compliance))
  CHECK_MEMBER(compliance, kTYPES, msg="correct compliance codes")
  # Check compatibility of compliance with data.
  alwaystakers = which(((1-Data$Z) * Data$W) == 1)
  nevertakers = which((Data$Z * (1-Data$W))==1)
  CHECK_MEMBER(alwaystakers, which(compliance==kTYPE_ALWAYSTAKER))
  CHECK_MEMBER(nevertakers, which(compliance==kTYPE_NEVERTAKER))
  
  s = 0
  
  for (t in kTYPES) {
    Ct = as.numeric(compliance==t)
    Qt1= joint.Ptyz(t, 1, Data$Y, Data$X, theta)
    Qt0= joint.Ptyz(t, 0, Data$Y, Data$X, theta)
    Zt1 = Ct * (Data$Z==1)  # units with Zi=1, Ci=t
    Zt0 = Ct * (Data$Z==0)  # units with Zi=0 Ci=t
    
    s = s + sum(Zt1 * log(Qt1) + Zt0 * log(Qt0))
  }
  return(s)
}

log.prior <- function(theta, X) {
  ## Defines the prior as in Hirano et. al.
  # Faster because it is using a vector format for all units.
  s = 0
  N = nrow(X)
  for(t in kTYPES)
    for(Zi in c(0, 1))
      for(y in c(0, 1))
        s = s + sum(log(joint.Ptyz(t=t, z=Zi, Y=rep(y, N), X=X, theta=theta)))
  return(2.5 * s / N)
}

Gibbs.sample.compliance <- function(Data, theta) {
  # Given the data, it will sample compliance types for the units
  #
  # Args:  DATA and THETA
  #
  # Returns:
  #   (N x 1) of compliance types (in kTYPES)
  #    P(C | Z, W, θ)
  # e.g. if θ=0 then for Z=0 and W=0 we should pick equally between Complier + Nevertaker
  #
  N = length(Data$Y)
  C = rep(NA, N)
  
  Cprobs = get.stratified.probs(Data, theta)
  # Qii = first column = complier probabilities
  #       second column = either always-taker or never-taker
  Q11 = matrix(c(Cprobs$Qc1, Cprobs$Qa1), nrow=N, ncol=2)
  Q00 = matrix(c(Cprobs$Qc0, Cprobs$Qn0), nrow=N, ncol=2)
  
  CA = c(kTYPE_COMPLIER, kTYPE_ALWAYSTAKER)
  CN = c(kTYPE_COMPLIER, kTYPE_NEVERTAKER)
  sample.type11 = apply(Q11, 1, function(x) {
    if(all(x==0)) return(sample(CA, size=1))
    sample(CA, size=1, prob=x)
  })
  sample.type00 = apply(Q00, 1, function(x) {
    if(all(x==0)) return(sample(CN, size=1))
    return(sample(CN, size=1, prob=x))
  })
  
  return(with(Cprobs, { S01 * kTYPE_ALWAYSTAKER + S10 * kTYPE_NEVERTAKER + 
                        S11 * sample.type11 + S00 * sample.type00 }))
}

run.mcmc <- function(theta0, Data, 
                     log.density, 
                     niters=1000, proposal.scale=0.8) {
  # Runs an MCMC for the problem.
  #
  # Args:
  #   Data = DATA object (See terminology)
  CHECK_Data(Data)
  # no. of units.
  N = length(Data$Y)
  
  # Compliance chain (N x niters)
  C.chain = matrix(NA, nrow=N, ncol=niters)
  C.chain[, 1] = sample(kTYPES, size=N, replace=T)
  
  # Chain of θ vectors (unpacked θ)
  V.chain = matrix(NA, nrow=length(unpack.theta(theta0)), ncol=niters)
  V.chain[, 1] = unpack.theta(theta0)
  print("Running MCMC")
  pb = txtProgressBar(style=3)
  # Keep chain statistics here.
  chain.stats = list(acc=c(0))
  
  for (t in 2:niters) {
    # 1. Sample the compliance types
    v.old = V.chain[, t-1]
    compliance.old = C.chain[, t-1]
    theta.old = pack.theta.vector(theta.vector=v.old, base.theta=theta0)
    compliance.new = Gibbs.sample.compliance(Data, theta.old)
    
    # 2. Sample new theta (Metropolis steps for each component)
    n0 = length(v.old)
    v.new = v.old
    # Runs intermediate MH steps.
    for(index in 1:n0) {
      v.cand = v.new
      v.cand[index] <- v.cand[index] + proposal.scale * rt(1, df=5)
      theta.cand = pack.theta.vector(theta.vector=v.cand, base.theta=theta0)
      log.acc = min(0, log.density(theta=theta.cand, compliance=compliance.new, Data=Data) - 
                       log.density(theta=theta.old, compliance=compliance.new, Data=Data))
      acc = exp(log.acc)
      if(!is.na(acc)) {
        ## can be NaN if it tries a very weird value (remember we are trying t-proposals)
        chain.stats$acc <- c(chain.stats$acc, acc)
        if(runif(1) < acc) {
          v.new = v.cand
        }
      }
    }
    # Now v.new = new θ vector
    
    m.acc = mean(tail(chain.stats$acc, 200))
    if(m.acc < 0.2) {
      print(sprintf("WARNING: Chain has small acceptance rate %.2f", m.acc))
    } else if(t %in% round(seq(10, niters, length.out=10))) {
      print(sprintf("iter =%d/%d Acceptance rate = %.2f%%", t, niters, 100 * m.acc))
    }
    
    # Updates
    C.chain[, t] <- compliance.new
    V.chain[, t] <- v.new
    setTxtProgressBar(pb, value=t/niters)
  }  # MCMC iterations.
  keep.iters = tail(1:niters, 0.9 * niters)
  print("Statistics of acceptance probability.")
  print(summary(chain.stats$acc))
  return(list(theta0=theta0, Data=Data,
              theta.vector.chain=V.chain[, keep.iters],
              compliance.chain=C.chain[, keep.iters]))
}

mcmc.full.posterior <- function(Data, niters=1000, proposal.scale=0.1) {
  # Runs MCMC with the prior only
  # We call it catalytic because it should match the ITT that is 
  # contained in the DATA itself
  theta0 = perturb.theta(empty.theta(), epsilon=3)
  log.density <- function(theta, compliance, Data) {
   return(log.likelihood.complete(theta, compliance, Data=Data) + 
          log.prior(theta=theta, X=Data$X))
  }
  
  return(run.mcmc(theta0=theta0, 
                  Data=Data,
                  log.density=log.density,
                  niters=niters,
                  proposal.scale=proposal.scale))
}

mcmc.baseline <- function(Data, niters=1000, proposal.scale=0.1) {
  # Runs MCMC with a bad density/proposal function
  # Just to make sure
  theta0 = perturb.theta(empty.theta(), epsilon=3)
  log.density <- function(theta, compliance, Data) {
    return(log(runif(1)))
  }
  
  return(run.mcmc(theta0=theta0, 
                  Data=Data,
                  log.density=log.density,
                  niters=niters,
                  proposal.scale=proposal.scale))
}

mcmc.prior <- function(nsamples=100) {
  # Tests whether the prior gives informative ITT estimates:
  #
  # 1. Sample θ from the prior
  # 2. Sample ITT | θ
  # 3. Iterate
  # 
  # Finally return the samples ITTs. These should be centered around 0
  # since we want the prior not to give much information about the ITTs.
  #
  theta0 = empty.theta()
  v.old = unpack.theta(theta0)
  p = length(v.old)
  v.new = v.old
  X = sample.X(1000)
  itt.samples = list()
  for(t in kTYPES) {
    itt.samples[[t]] = rep(NA, nsamples)
  }
  chain.stats = list(acc=c())
  pb = txtProgressBar(style=3)
  for (t in 1:nsamples) {
    # Gibbs-MH for a new θ
    v.cand = v.new
    for(i in 1:p) {
      v.cand = v.new
      v.cand[i] <- v.cand[i] + rt(1, df=5)
      theta.cand = pack.theta.vector(v.cand, base.theta=theta0)
      theta.new = pack.theta.vector(v.new, base.theta=theta0)
      log.acc = min(0, log.prior(theta=theta.cand, X=X) - 
                       log.prior(theta=theta.new, X=X))
      # print(log.acc)
      acc = exp(log.acc)
      if(!is.na(acc)) {
        ## can be NaN if it tries a very weird value (remember we are trying t-proposals)
        chain.stats$acc <- c(chain.stats$acc, acc)
        if(runif(1) < acc) {
          v.new = v.cand
        }
      }
    }
    # Done with sampling θ.new
    # v.new has the sampled parameter vector
    v.old = v.new
    theta.new = pack.theta.vector(v.new, base.theta=theta0)
    itt = sample.ITT.theta(theta=theta.new, X=X, nsamples=1)
    for(j in 1:length(itt)) {
      itt.samples[[j]][t] <- itt[[j]][1]
    }
    setTxtProgressBar(pb, value=t/nsamples)
  }
  print("Acceptance summary")
  print(summary(chain.stats$acc))
  return(itt.samples)
}

##  Sample ITT
sample.ITT.theta <- function(theta, X, nsamples=100) {
  N = nrow(X)
  Psi = get.compliance.probs(theta, X)
  ITT = list()
  for(t in kTYPES) {
    ITT[[t]] = rep(0, nsamples)
  }
  
  for(i in 1:nsamples) {
    compliance = apply(Psi, 1, function(x) sample(kTYPES, size=1, prob=x))
    CHECK_EQ(length(compliance), N)
    for(t in kTYPES) {
      It = as.numeric(compliance==t)
      Nt = sum(It)
      Yt1 = rbinom(N, size=1, prob=logistic(X %*% get.beta.vector(t, z=1, theta=theta)))
      Yt0 = rbinom(N, size=1, prob=logistic(X %*% get.beta.vector(t, z=0, theta=theta)))
      CHECK_EQ(length(Yt1), N)
      estimate = ifelse(Nt==0, 0, sum(It * (Yt1-Yt0)) / Nt)
      ITT[[t]][i] <- estimate
    }
  }
  return(ITT)
}

sample.ITT.mcmc <- function(mcmc.out, use.obs=T) {
  # Given a matrix of θ, Compliance, sample the ITT effects.
  # 1. For every θj, Cj (j=iteration):
  #    1a. Sample Y1, Y0 through β_tz (z=1, 0 resp.) for every t
  #    1b. Then impute Y1, Y0 for type t 
  #    1c. Compute ITT effects only for those of type t (Cj_i==t)
  #
  # Returns:
  #   LIST(n, c, a)  for nevertakers, compliers and alwaystakers respectively.
  compliance.matrix = mcmc.out$compliance.chain
  theta.vector.matrix = mcmc.out$theta.vector.chain
  Data = mcmc.out$Data
  theta0 = mcmc.out$theta0
  
  CHECK_Data(Data)
  warning("No test for sample.iTT")
  N = length(Data$Y)
  niters = ncol(compliance.matrix)
  CHECK_EQ(ncol(theta.vector.matrix), niters)
  CHECK_EQ(dim(compliance.matrix), c(N, niters))
  
  # return object.
  out <- list(N=c(), C=c(), A=c())
  add.value <- function(compliance.type, value) {
    name = kComplianceTypes[compliance.type]
    out[[name]] <- c(out[[name]], value)
    return(out)
  }
  
  for(iter in 1:ncol(compliance.matrix)) {
    v.iter <- theta.vector.matrix[, iter]
    C.iter <- compliance.matrix[, iter]
    theta.iter = pack.theta.vector(theta.vector=v.iter, base.theta=theta0)
    
    for(t in kTYPES) {
      # Ct = (0, 1, 1, 0, 0..) where 1 if i is of type t
      Ct = as.numeric(C.iter == t)
      # Nt = #units of type t
      Nt = sum(Ct)
      if(Nt==0) {
        out <- add.value(t, 0)
        next
      }

      # impute potential outcomes
      bt1 = get.beta.vector(t=t, z=1, theta=theta.iter)
      preds1 = Data$X %*% bt1
      # Sample Y | Ci=t, Zi=1, Xi, θ)
      Y1.sample = rbinom(N, size=1, prob=logistic(preds1))
      # Sample Y | Ci=t, Zi=0, Xi, θ
      bt0 = get.beta.vector(t=t, z=0, theta=theta.iter)
      preds0 = Data$X %*% bt0
      Y0.sample = rbinom(N, size=1, prob=logistic(preds0))
      # encouraged
      Z1 = as.numeric(Data$Z == 1)
      Yobs = Data$Y
      #  Impute potential outcomes
      # Yi(Zi=1, Ci=t) = {
      #     0  if Ci != t
      #        or
      #    Yobs_i if Z=1 
      #    Y1_i  if Z=0
      Y1.imputed = Y1.sample
      Y0.imputed = Y0.sample
      if(use.obs) {
        Y1.imputed = Ct * (Yobs * Z1 + Y1.sample * (1-Z1))
        Y0.imputed = Ct * (Yobs * (1-Z1) + Y0.sample *  Z1)
      }
      out <- add.value(t, sum(Y1.imputed - Y0.imputed) / Nt)
    }
  }
  return(out)
}

summarize.ITT.samples <- function(ITT, theta) {
  par(mfrow=c(3,1))
  X = sample.X(1000)
  for(t in 1:length(ITT)) {
    print(sprintf("Summary for compliance %s", t))
    print(summary(ITT[[t]]))
    print("Heuristic derivation:")
    m1 = mean(logistic(X %*% get.beta.vector(t, z=1, theta=theta)))
    m0 = mean(logistic(X %*% get.beta.vector(t, z=0, theta=theta)))
    print(sprintf("Heuristic difference (ITT effect) = %.3f", m1-m0))
    hist(ITT[[t]], main=sprintf("Compliance type %s", kComplianceTypes[t]))
  }
}
