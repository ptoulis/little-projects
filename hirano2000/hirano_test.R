## Tests for hirano.R
source("hirano.R")
# Controls the "example.theta" values.
kTestIntercepts = c(0.2, -1.2, 0.25, 5.1, -2.1, 0.776)
kTestSlope1 <- -1.2
kTestSlope2 <- 1.9
kTestPsi.c <- c(0, 0, 0)
kTestPsi.a <- c(0, 0, 0)

example.theta <- function() {
  ## For the tests, fix values for the parameters.
  #
  # Creates an example parameter vector.
  # Only sets β_c1 = complier, Z=1
  #      and  ψ_c = complier
  #
  theta <- empty.theta()
  # Set slope for X1
  cnt = 1
  for(t in kTYPES)
    for(z in kTREATS) {
      theta$beta[t, z, 1] <- kTestIntercepts[cnt]
      theta$beta[t, z, 2] <- kTestSlope1
      theta$beta[t, z, 3] <- kTestSlope2
      cnt = cnt + 1
    }
  # ψ_c
  theta$psi[kTYPE_COMPLIER, ] <- kTestPsi.c
  theta$psi[kTYPE_ALWAYSTAKER, ] <- kTestPsi.a
  
  CHECK_Theta(theta)
  return(theta)
}

test.get.compliance.probs <- function() {
  theta <- empty.theta()
  X = sample.X(nUnits=10)
  Psi = get.compliance.probs(theta=theta, X=X)
  CHECK_TRUE(all(ncol(X) * Psi == 1))
  # Test 2.
  theta$psi[,] <- 100
  Psi = get.compliance.probs(theta=theta, X=X)
  CHECK_TRUE(all(ncol(X) * Psi == 1))
  # Test 3.
  # Compliers coefficients are big and X > 0
  # so we only get compliers.
  theta <- empty.theta()
  theta$psi[kTYPE_COMPLIER, ] <- 100
  X = abs(X)
  Psi = get.compliance.probs(theta=theta, X=X)
  CHECK_TRUE(all(Psi[, kTYPE_COMPLIER] > 0.99))
}
  
test.get.stratified.probs <- function() {
  Data = list() 
  #     (C,N)    A   N  (C,A) A  (C, N) 
  Data$Z = c(0,  0,  1,  1,   0,   0)
  Data$W = c(0,  1,  0,  1,   1,   0)
  Data$Y = c(1,  0,  1,  0 ,  1 ,  0)
  
  N = length(Data$Z)
  nTypes = length(kTYPES)
  Data$X = abs(sample.X(N)) # X  > 0
  
  #  β_tz = 0  except for β_c1 = (10, 10, 10)   -- big values.
  # So this means that Yi=1 is more probable for Ci=Complier, Zi=1
  # anything else is equal chance.
  theta = empty.theta()
  theta$same.slopes = F
  theta$beta[kTYPE_COMPLIER, 1, ] <- rep(0, kXdim)
  theta$beta[kTYPE_COMPLIER, 2, ] <- rep(100, kXdim)
  theta$psi[kTYPE_ALWAYSTAKER, ] <- 1
  Psi = get.compliance.probs(theta=theta, X=Data$X)
  
  Cprobs = get.stratified.probs(Data=Data, theta=theta)
  CHECK_EQ(which(Cprobs$S01==1), c(2, 5))
  CHECK_EQ(which(Cprobs$S11==1), c(4))
  CHECK_NEAR(Cprobs$Qn1, Psi[, kTYPE_NEVERTAKER] * (1/2))
  CHECK_NEAR(Cprobs$Qa0, Psi[, kTYPE_ALWAYSTAKER] * (1/2))
  # Compliers by definition above have Y=1 a.s when Z=1
  CHECK_NEAR(Cprobs$Qc1, Psi[, kTYPE_COMPLIER] * Data$Y)
  CHECK_NEAR(Cprobs$Qc0, Psi[, kTYPE_COMPLIER] * (1/2))
}

test.sample.compliance <- function() {
  theta = empty.theta()
  nUnits = 100
  D = gen.data(nUnits, theta)
  nsamples = 1000
  
  count.types <- function(M, t) {
    apply(M, 2, function(x) length(which(x==t)))
  }
  
  # 1. All never-takers
  D$Z = rep(1, nUnits)
  D$W = rep(0, nUnits)
  C = matrix(replicate(nsamples, { sample.compliance(D, theta) }),
             nrow=length(D$Y), ncol=nsamples)
  
  type.n = count.types(C, kTYPE_NEVERTAKER)
  CHECK_SETEQ(type.n, c(nUnits))
  
  # 2. All always-takers
  D$Z = rep(0, nUnits)
  D$W = rep(1, nUnits)
  C = matrix(replicate(nsamples, { sample.compliance(D, theta) }),
             nrow=length(D$Y), ncol=nsamples)
  
  type.a = count.types(C, kTYPE_ALWAYSTAKER)
  type.c =  count.types(C, kTYPE_COMPLIER)
  CHECK_SETEQ(type.a, c(nUnits))
  CHECK_SETEQ(type.n, c(nUnits))
  
  # 3. Compliers or always takers
  D$Z = rep(1, nUnits)
  D$W = rep(1, nUnits)
  C = matrix(replicate(nsamples, { sample.compliance(D, theta) }),
             nrow=length(D$Y), ncol=nsamples)
  type.a = count.types(C, kTYPE_ALWAYSTAKER)
  type.c = count.types(C, kTYPE_COMPLIER)
  CHECK_SETEQ(type.a + type.c, c(nUnits), "Only C and A")
  # simple statistical test.
  CHECK_TRUE(prop.test(x=type.a, n=rep(nUnits, nsamples), conf.level=0.99)$p.value > 0.01, msg="C or A with equal prob.")
}

test.pack.theta <- function() {
  # Tests the packing/unpacking of the theta.
  theta1 <- example.theta()
  theta2 <- empty.theta()
  
  v1 = unpack.theta(theta1)

  CHECK_EQ(v1, c(kTestIntercepts, kTestSlope1, kTestSlope2, kTestPsi.c, kTestPsi.a))
  CHECK_EQ(get.beta.intercepts(theta1)[, 1], v1[c(1, 3, 5)])
  CHECK_EQ(get.beta.slopes(theta1),  v1[7:8])
  v2 = unpack.theta(theta2)

  t1 = pack.theta.vector(v1, theta1)
  t2 = pack.theta.vector(v2, theta2)
  
  CHECK_EQ_THETA(t1, theta1)
  CHECK_EQ_THETA(t2, theta2)
}

test.enable.restriction <- function() {
  theta1 <- example.theta()
  theta2 <- enable.exclusion.restriction(restriction="A", theta=theta1)
  theta2 <- enable.exclusion.restriction(restriction="N", theta=theta2)
  
  v1 = unpack.theta(theta1)
  v2 = unpack.theta(theta2)
  CHECK_EQ(v2, c(kTestIntercepts[-c(1, 5)], kTestSlope1, kTestSlope2, kTestPsi.c, kTestPsi.a))
  t1 = pack.theta.vector(v1, theta1)
  t2 = pack.theta.vector(v2, theta2)

  CHECK_EQ_THETA(t1, theta1)
  CHECK_EQ_THETA(t2, theta2)
  
}

test.log.likelihood <- function(loglik.or.prior="loglik",
                             ntrials=100, 
                             nDataSize=1000, 
                             epsilon=0.1) {
  # Suggested runs:
  # 1. test.log.density("loglik", ntrials=200, nDataSize=3000, epsilon=0.1)
  # 2. test.log.density("prior", ntrials=200, nDataSize=3000, epsilon=0.1)
  #
  # If everything is OK, then 1. should not get many log-likelihood improvements
  # This means that the LL is maximized around the ground truth (one used to generate the data)
  # Also 2. (prior) should yield many improvements because we need the prior to be
  # non-informative.
  #
  # We should expect the following:
  #   1. Increasing nTrials, will get improved estimates (assessments)
  #   2. Increasing nDataSize will improve sharpness of likelihood and
  #      so it will be harder to get improvements over the ground-truth
  #   3. Increasing "epsilon" increases the step size and so we will get less
  #      successful trials. Smaller sizes will be more successful that they will
  #     yield smaller %improvements. So, there is a trade-off one needs to explore.
  #
  theta <- example.theta()
  theta <- enable.exclusion.restriction("A", theta)
  theta <- enable.exclusion.restriction("N", theta)
  theta$psi[3, ] <- runif(3, min=0, max=1.5)  # add some more variation
  D = gen.data(nUnits=nDataSize, theta=theta)
  # Define neg-log likelihood
  
  fn <- function(v) -log.likelihood(pack.theta.vector(v, base.theta=theta), Data=D)
  if (loglik.or.prior == "prior")
    fn <- function(v) -log.prior(theta=pack.theta.vector(v, base.theta=theta), X=D$X)
  
  v0 = unpack.theta(theta)
  # uncomment this if you want to test a random point
  # e.g.  test.log.density("loglik", ntrials=200, nDataSize=5000, epsilon=0.1)
  # should give vectors that mostly increase the likelihood.
  # v0 = runif(length(v0), min=-2, max=-2)
  # warning("testing neighboring points of a random point")
  
  ll0 = fn(v0)
  print(sprintf("K=%d Ground-truth Neg Log-Lik=%.2f", length(v0), ll0))
  
  # keeps track which trial improves the ll
  found.better = c()
  # best improvement achieved (%)
  best.improvement = 0
  # best ll
  best.ll = 0
  
  pb <- txtProgressBar(style=3)
  
  for(i in 1:ntrials) {
    # 1. Try simple perturbation
    v1 = v0 + runif(length(v0), min=-epsilon, max=epsilon)
    ll = fn(v1)
   
    if(ll < ll0)  {
      found.better[i] <- 1
      best.improvement = 100 * abs( (ll-ll0) / ll0)
      best.ll = ll
    } else {
      found.better[i] <- 0
    }
    setTxtProgressBar(pb, value= i/ntrials)
  }
  
  better.prop = 100 * mean(found.better)
  # If we test log-likelihood, then we need to make sure 
  # that around the ground-truth small perturbations do not increase
  # the likelihood. Note: The "OK"/"FAIL" flags are arbitrary.
  success = ifelse(better.prop < 10, "OK", "FAIL")
  # If we are testing the prior, we need to make sure, we improve the likelihood
  # and thus that the prior is not informative.
  if (loglik.or.prior=="prior")
    success = ifelse(better.prop > 40, "OK", "FAIL")
  
  print(sprintf("Total %.2f%% of random perturbations around the ground truth had better LogL. This is [%s]", 
                better.prop, success))
  # Even if you get a "FAIL" it is only important to see how much the LL is increased.
  # That is, even if many θ values get better likelihood, if the improvement 
  # is not great, this is OK when testing the log-likelihood.
  print(sprintf("LL0=%.3f LL1=%.2f  Improvement=%.2f%%", ll0, best.ll, best.improvement))
  if(success=="FAIL" && best.improvement < 1) {
    print(sprintf("Best improvement < 1%%. [FAIL] tag is unreliable. Try increasing sample size or epsilon"))
  }
}

test.run.mcmc.catalytic <- function() {
  theta <- example.theta()
  D = gen.data(500, theta=theta)
  for (t in kTYPES) {
    print(sprintf("ITT effect for %s = %.3f", kComplianceTypes[t], D$ITT[[t]]))
  }
  mcmc.out = run.mcmc.catalytic(Data=D, niters=10000)
  ITT = sample.ITT(mcmc.out)
  par(mfrow=c(3, 1))
  for(t in names(ITT)) {
    s = sprintf("MCMC: ITT effect for %s", t)
    samples = ITT[[t]]
    print(s)
    print(summary(samples))
    print(sprintf("mean=%.2f se=%.2f", mean(samples), boot(samples)))
    hist(samples, main=s)
  }
}

test.sample.ITT.theta <- function() {
  theta0 = empty.theta()
  X = sample.X(1000)
  samp = sample.ITT.theta(theta0, X=X, nsamples=1000)
  CHECK_MU0(samp[[2]], 0, msg="Complier")
  CHECK_MU0(samp[[3]], 0, msg="alwaystaker effect=0")
  theta0$beta[kTYPE_ALWAYSTAKER, 1, ] <- -10
  theta0$beta[kTYPE_ALWAYSTAKER, 2, ] <- 10
  samp = sample.ITT.theta(theta0, X=X, nsamples=1000)
  CHECK_MU0(samp[[3]], 1, msg="always taker effect==1")
}