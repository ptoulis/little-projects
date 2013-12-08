# terminology file.
#
#   θ = {beta, psi, restrictions} = parameter vector
#     β = {β_tz} = (p x 1) vector of coefficients for type t, assignment z
#     ψ = {ψ_t} = (p x 1) vector of coeffs for type t
#     restrictions = {A:T or F,  N: T or F}
#     where R = exclusion restrictions. If R$A=T, then there are 
#       exclusion restrictions for the always-takers (A)
# We represent θ (THETA) as a composite LIST(beta, psi)
# each being a multi-dim array.
# In the paper (Hirano et. al.) they impose several constraints
# on the individual parameters. These are checked in CHECK_Theta()
#
# TYPES are compliance types. kTYPE_* has all defined types.
# TREATS are treatment arms. Here is is ACTIVE (Z=1) or CONTROL(Z=0)
# 
rm(list=ls())
source("../../r-toolkit/checks.R")

# Dimension of the covariates.
kXdim <- 3
# Compliance
#                   1     2    3
kComplianceTypes <- c("N", "C", "A")
kTYPES <- c(1, 2, 3)
kTYPE_NEVERTAKER <- 1
kTYPE_COMPLIER <- 2
kTYPE_ALWAYSTAKER <- 3

# Treatments
kTREATS <- c(1, 2)
kTreatmentArms <- c("Z0", "Z1")
kTREAT_CONTROL <- 1
kTREAT_ACTIVE <- 2
get.treatment.index <- function(zi.assignment) {
  CHECK_MEMBER(zi.assignment, c(0, 1))
  return(zi.assignment + 1)
}


sample.X <- function(nUnits, verbose=F) {
  # Samples covariates for a unit
  X <- matrix(runif(kXdim * nUnits, min=0, max=3), nrow=nUnits, ncol=kXdim)
  X[, 1] <- 1  # set to 1.
  return(X)
}

CHECK_Theta <- function(theta) {
  # Checks the THETA object
  # How we check:
  #   1. Check it has properties (beta, psi, restrictions)
  #   2. Check the slopes are equal β_tz1 = const.   β_tz2 =  const.
  #   3. Check the exclusion restrictions.
  CHECK_MEMBER(names(theta), c("beta", "psi", "restrictions", "same.slopes"))
  nC = length(kComplianceTypes)
  nZ = length(kTreatmentArms)
  CHECK_EQ(1:nC, kTYPES)  # 1:3 == 1:3 (kTYPES is defined in terminology for convenience)
  CHECK_EQ(1:nZ, kTREATS)  # 1:2 == 1:2 
  # dimensions
  CHECK_EQ(dim(theta$beta), c(nC, nZ, kXdim))
  CHECK_EQ(dim(theta$psi), c(nC, kXdim))
  # ψ_n = 0
  CHECK_EQ(theta$psi[kTYPE_NEVERTAKER,], rep(0, kXdim), "psi_n=0")
  x = theta$beta[kTYPE_COMPLIER, kTREAT_ACTIVE, 2]  # slope 1
  y = theta$beta[kTYPE_COMPLIER, kTREAT_ACTIVE, 3]  # slope 2
  # slopes should be equal
  if(theta$same.slopes) {
    for (t in kTYPES)
      for(z in kTREATS) {
        CHECK_EQ(theta$beta[t, z, 2], x, "Check slopes for p=1")
        CHECK_EQ(theta$beta[t, z, 3], y, "Check slopes for p=2")
      }
  }
  
  # Check restrictions
  CHECK_MEMBER(names(theta$restrictions), c("A", "N"))
  if(theta$restrictions$A) {
    CHECK_EQ(theta$beta[kTYPE_ALWAYSTAKER, 1, 1], 
             theta$beta[kTYPE_ALWAYSTAKER, 2, 1],
             msg="Exclusion restriction for Always-takers")
  }
  if(theta$restrictions$N) {
    CHECK_EQ(theta$beta[kTYPE_NEVERTAKER, 1, 1], 
             theta$beta[kTYPE_NEVERTAKER, 2, 1],
             msg="Exclusion restriction for Never-takers")
  }
}

empty.theta <- function() {
  # Create an empty THETA parameter object (all set to 0)
  # Parameters are:
  #  ψ_t = (,,) px1 vector of coefficients for types
  #  β_tz = (,,,)  p x 1 vector of coefficients for
  #          , type=t, assignment Z=z
  # params for types (n=fixed at zero)
  psi <- array(0, dim=c(length(kComplianceTypes),
                        kXdim))
  # params for outcomes
  # β_tz = (px1) vector
  betas = array(0, dim=c(length(kComplianceTypes),
                         length(kTreatmentArms),
                         kXdim))
  
  restrictions = list(A=FALSE, N=FALSE)
  theta = list(beta=betas, psi=psi,
               restrictions=restrictions,
               same.slopes=T)
  CHECK_Theta(theta)
  return(theta)
}

# Access functions to parameters.
get.beta.vector <- function(t, z, theta) {
  # Returns β_tz
  CHECK_MEMBER(t, kTYPES)
  CHECK_MEMBER(z, c(0, 1))
  return(theta$beta[t, z + 1, ])
}

get.beta.intercepts <- function(theta) {
  CHECK_Theta(theta)
  return(theta$beta[,,1])
}

get.beta.slopes <- function(theta) {
  CHECK_Theta(theta)
  return(c(theta$beta[1, 1, 2], theta$beta[1, 1, 3]))
}

get.restrictions <- function(theta) {
  CHECK_Theta(theta)
  return(c(theta$restrictions$N, theta$restrictions$A))
}

get.psi.matrix <- function(theta) {
  return(theta$psi)
}

get.psi.vector <- function(t, theta) {
  # Returns ψ_t
  CHECK_MEMBER(t, kTYPES)
  return(theta$psi[t, ])
}

logistic <- function(x) {
  x[x < -100] = -100
  x[x > 100] = 100
  exp(x) / (1+exp(x))
}

logit <- function(x) {
  CHECK_INTERVAL(x, 0, 1)
  if(x == 0) return(-Inf)
  if(x == 1) return(Inf)
  return(log(x / (1-x)))
}

CHECK_Data <- function(Data) {
  # Checks the DATA object.
  # DATA = {Z, W, Y, X}
  CHECK_MEMBER(names(Data), c("Z", "W", "X", "X0", "Y", "ITT", "C"), msg="Data members")
  N = length(Data$Z)
  CHECK_EQ(length(Data$Y), N, msg="Length of Yobs")
  CHECK_EQ(length(Data$W), N, msg="Length of Wobs")
  CHECK_EQ(nrow(Data$X), N, msg="Length of X")
  CHECK_MEMBER(c(Data$Y, Data$W, Data$Z), c(0, 1))
}

normalize.data <- function(D) {
  for(i in 2:ncol(D$X)) {
    mu = mean(D$X[, i])
    se = sd(D$X[, i])
    D$X[,i] = (D$X[, i] - mu) / se
  }
  return(D)
}

## Data
read.data <- function() {
  CHECK_TRUE(file.exists("flu240.txt"))
  data = read.table("flu240.txt", header=T)
  data$const <- rep(1, nrow(data))
  X = as.matrix(subset(data, select=c("const", "age", "copd")))
  Y = data$outcome
  Z = data$treatment.assigned
  W = data$treatment.received
  D = list(X=X, Y=Y, W=W, Z=Z, ITT=list())
  D$X0 <- X
  D = normalize.data(D)
  return(D)
}

gen.data <- function(nUnits=1000, theta) {
  # Generates a DATA object given the input
  #
  # Args:
  #   nUnits = #units to be sampled
  #   theta = parameters
  #
  # Returns:
  #   A DATA object (Y, X, W, Z)
  N = nUnits
  CHECK_EQ(nUnits %% 2, 0, msg="no. of units should be even.")
  Z = sample(c(rep(0, N/2), rep(1, N/2)))
  W = c()
  X = sample.X(nUnits)
  Y = c()
  Y1 = c()
  Y0 = c()
  Compl = c()  # compliance
  ComplianceProbs <- get.compliance.probs(theta=theta, X)
  for(i in 1:N) {
    Zi = Z[i]
    # 1. Sample covariates
    Xi = X[i, ]
    # 2. Sample compliance
    ci.probs <- ComplianceProbs[i, ]
    Ci = sample(kTYPES, size=1, prob=ci.probs)
    Wi = NA
    if(Ci==kTYPE_COMPLIER) {
      Wi = Zi
    } else if(Ci==kTYPE_ALWAYSTAKER) {
      Wi = 1
    } else if(Ci==kTYPE_NEVERTAKER) {
      Wi = 0
    }
    
    # 3. Sample Yi
    b = get.beta.vector(t=Ci, z=Zi, theta=theta)
    Yi.obs <- rbinom(1, size=1, prob=logistic(sum(b * Xi)))
    
    Zi.counter = 1-Zi
    b.counter = get.beta.vector(t=Ci, z=Zi.counter, theta=theta)
    Yi.counter <- rbinom(1, size=1, prob=logistic(sum(b.counter * Xi)))
    
    # 4. Load
    W[i] <- Wi
    X[i, ] <- Xi
    Y[i] <- Yi.obs
    Y1[i] <- ifelse(Zi==1, Yi.obs, Yi.counter)
    Y0[i] <- ifelse(Zi==1, Yi.counter, Yi.obs)
    Compl[i] <- Ci
  }
  
  ITT = list()
  for(t in kTYPES) {
    ct = (Compl==t)
    Nt = sum(ct)
    ITT[[t]] <- ifelse(Nt==0, 0, sum(ct * (Y1 - Y0)) / Nt) 
  }
  
  Data = list(Z=Z, W=W, X=X, Y=Y, ITT=ITT, C=Compl)
  CHECK_Data(Data)
  return(Data)  
}

perturb.theta <- function(theta, epsilon=0.1) {
  v0 = unpack.theta(theta=theta)
  v1 = v0 + runif(length(v0), min=-epsilon, max=epsilon)
  return(pack.theta.vector(v1, base.theta=theta))
}

CHECK_EQ_THETA <- function(theta1, theta2) {
  CHECK_Theta(theta1)
  CHECK_Theta(theta2)
  for(i in 1:kXdim)
    CHECK_TRUE(all(theta1$beta[,,i] == theta2$beta[,,i]))
  CHECK_TRUE(all(theta1$psi == theta2$psi))
  CHECK_EQ(get.restrictions(theta1), get.restrictions(theta2))
}

unpack.theta <- function(theta) {
  # Takes an object THETA (parameters)
  # and unpacks into a vector. 
  # This will be useful for MCMC, or EM etc.
  # Unpack(theta) = (β_intercepts,  slopes,  ψ)
  # We know that:
  #   ψ = (ψc, ψa) = always  regardless of restruitions
  #   slopes = (x1, x2) = same for all (t, z) pairs.
  #   intercepts  =   6x1  if no restrictions
  #                   5x1 if 1 restriction
  #                   4x1 if both restrictions
  # Returns a vector.
  CHECK_Theta(theta)
  all.intercepts <- c() 
  for(t in kTYPES) {
    for(z in c(0, 1)) {
      beta.tz0 = get.beta.vector(t, z, theta)[1]
      all.intercepts <- c(all.intercepts, beta.tz0)
    }
  }
  inter.index= c(1, 2, 3, 4, 5, 6)  # be explicit
  if(theta$restrictions$N) {
    inter.index = setdiff(inter.index, c(2)) # intercepts are the same for Never-takers
  }
  if(theta$restrictions$A) {
    inter.index = setdiff(inter.index, c(6))  # intercepts are the same for Always-takers
  }
  
  intercepts = all.intercepts[inter.index]
  slopes = c(theta$beta[1, 1, 2], theta$beta[1, 1, 3])
  psis = c(get.psi.vector(t=kTYPE_COMPLIER, theta),
           get.psi.vector(t=kTYPE_ALWAYSTAKER, theta))
  
  theta.vec <- c(intercepts, slopes, psis)
  return(theta.vec)
}

pack.theta.vector <- function(theta.vector, base.theta) {
  # Takes a vector and packs it into a THETA params object.
  #
  # Should be that pack(unpack(θ)) = θ
  # Returns: θ (THETA) paramaters
  CHECK_Theta(base.theta)
  
  theta <- empty.theta()
  theta$restrictions = base.theta$restrictions
  
  if(theta$restrictions$N) {
    theta.vector <- append(theta.vector, theta.vector[1], after=1)
  }
  if(theta$restrictions$A) {
    theta.vector <- append(theta.vector, theta.vector[5], after=5)
  }
  
  i = 1
  for(t in kTYPES) {
    for(z in kTREATS) {
      beta.tz0 = theta.vector[i]
      theta$beta[t, z, 1] <- beta.tz0
      theta$beta[t, z, 2] <- theta.vector[6 + 1]
      theta$beta[t, z, 3] <- theta.vector[6 + 2]
      i <- i +1
    }
  }
  
  theta$psi[kTYPE_NEVERTAKER, ] <- rep(0, kXdim)
  theta$psi[kTYPE_COMPLIER, ] <- theta.vector[seq(6+3, 6+3+kXdim-1)]
  theta$psi[kTYPE_ALWAYSTAKER, ] <- tail(theta.vector, kXdim)
  CHECK_Theta(theta)
  return(theta)
}

enable.exclusion.restriction <- function(restriction, theta) {
  # Sets the exclusion restrictions.
  # Basically sets the flag (T or F) which will be used from pack/unpack functions.
  CHECK_MEMBER(restriction, c("N", "A"))
  CHECK_Theta(theta)
  # there is no going-back after this
  theta$restrictions[[restriction]] <- TRUE
  t = which(kComplianceTypes==restriction)
  # set intercepts to be equal
  theta$beta[t, 1, 1] = theta$beta[t, 2, 1]
  # Now the entire vector should be equal
  CHECK_EQ(get.beta.vector(t, z=0, theta=theta),
           get.beta.vector(t, z=1, theta=theta),
           msg="After enabliing exclusion-restriction")
  CHECK_Theta(theta)
  return(theta)
}

boot <- function(x) sd(replicate(1000, { mean(sample(x, size=length(x), replace=T)) }))