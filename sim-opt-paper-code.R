############# Rinott ############# 

Rinotth <- function(k, n0, pstar, conf, rep)
{
  # last update 7/13/2015
  # function to return estimate of Rinott h and upper confidence bound on it
  # k = number of systems
  # n0 = first-stage sample size
  # pstar = 1-alpha value (PCS)
  # rep = number of replications to use for estimate
  # conf = confidence level on upper bound
  #
  # set.seed(13)
  Z <- matrix(rnorm((k-1) * rep), ncol = k-1)
  Y <- matrix(rchisq((k-1) * rep, n0 - 1.), ncol = k-1)
  C <- matrix(rchisq(rep, n0 - 1.), ncol = 1.)
  Cmat <- matrix(C, ncol = k-1, nrow = rep, byrow = F)
  denom <- sqrt((n0 - 1.) * (1./Y + 1./Cmat))
  H <- sort(apply(Z * denom, 1., max))
  Hstar <- quantile(H, pstar)
  upper <- ceiling(pstar * rep + qnorm(conf) * sqrt(pstar * (1. - pstar) *
                                                      rep) + 0.5)
  Hupper <- H[upper]
  list(h = Hstar, UCB = Hupper)
}

Rinott <- function(k, alpha, n0, delta){
  # implements Rinott's procedure
  # k = number of systems
  # 1-alpha = desired PCS
  # n0 = first-stage sample size
  # delta = indifference-zone parameter
  # note: uses 99% UCB for Rinott's h
  h <- Rinotth(k, n0, 1-alpha, 0.99, 10000)$UCB
  Ybar <- NULL
  Vars <- NULL
  Ns <- NULL
  # loop through the k systems
  for (x in 1:k){
    Y <- MySim(x, n0)
    S2 <- var(Y)
    N <- ceiling(h^2*S2/delta^2)
    if (N > n0){
      Y <- c(Y, MySim(x, N-n0))
    }
    Ybar <- c(Ybar, mean(Y))
    Vars <- c(Vars, S2)
    N <- max(N, n0)
    Ns <- c(Ns, N)
  }
  list(Best = which.max(Ybar), Ybar = Ybar, Var = Vars, N = Ns)
}


############# NSGS ############# 

NSGS <- function(k, alpha, n0, delta){
  # implements NSGS
  # k = number of systems
  # 1-alpha = desired PCS
  # n0 = first-stage sample size
  # delta = indifference-zone parameter
  # note: uses 99% UCB for Rinott's h
  h <- Rinotth(k, n0, 1-alpha/2, 0.99, 10000)$UCB
  # first apply subset selection
  Yall <- NULL
  for (x in 1:k){
    Yall <- cbind(Yall, MySim(x, n0))
  }
  Ybar <- apply(Yall, 2, mean)
  S2 <- apply(Yall, 2, var)/n0
  tval <- qt((1-alpha/2)^(1/(k-1)), df = n0-1)
  Subset <- 1:k
  for (i in 1:k){
    for (j in 1:k){
      if (Ybar[i] < (Ybar[j]-tval*sqrt(S2[i] + S2[j]))){
        Subset[i] <- 0
        break
      }
    }
  }
  # next apply Rinott to the survivors
  Survivors <- (1:k)[Subset != 0]
  # loop through the survivors
  Ybars <- NULL
  Ns <- NULL
  if (length(Survivors) > 1){
    for (x in Survivors){
      N <- ceiling(h^2*S2[x]/delta^2)
      Y <- Yall[,x]
      if (N > n0){
        Y <- c(Y, MySim(x, N-n0))
      }
      Ybars <- c(Ybars, mean(Y))
      N <- max(N, n0)
      Ns <- c(Ns, N)
    }
    Best = Survivors[which.max(Ybars)]
  }
  else{
    Best <- Survivors
    Ybars <- Ybar[Survivors]
    Ns <- n0
  }
  list(Best = Best, Survivors = Survivors, Ybar = Ybars, N = Ns)
}
