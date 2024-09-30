#Main script
library(igraph)
library(netseg)

#Output
n0 <- 100
foclist <- matrix(data = 0, nrow = n0, ncol = 3)

mechanism <- c("dolr", "dolrp")

#Simulation specs
nodetype <- list(c("S", "T"),c("SU", "SP", "TU", "TP"),c("S", "T"),c("SU", "SP", "TU", "TP"),c("S", "T"),c("SU", "SP", "TU", "TP"))
nodenum <- 100
eom <- 1
invec <- 1

opvec <- rep(3,nodenum)
#typevec <- sample(nodetype[[1]], nodenum, replace = TRUE)

#Simulation parameters
rho <- c(0.25,0.5,0.75)
bel <- c(0.25,0.5,0.75)

#output prep
resout <- matrix(data = 0, nrow = 9, ncol = 11)

rcount <- expand.grid(rho,bel)

for (cc in 1:nrow(rcount)) {
  
  #graph generation
  rph <- sample_gnp(nodenum, rcount[cc,1], directed = TRUE, loops = FALSE)
  
  adjrph <- as.matrix(as_adjacency_matrix(rph))
  
  #Deciding strategy
  focstrat <- which(adjrph[1,] == 1)
  
  if(length(focstrat)<=5){
    ridx <- sample(nrow(adjrph), 4)
    ridx <- ridx[! ridx %in% focstrat]
    adjrph[1,ridx] = 1
  }
  
  #Updating the graph after adding random edges
  rph <- graph_from_adjacency_matrix(adjrph, mode = "directed")
  
  #Graph properties
  props <- c(mean_distance(rph, directed = TRUE), mean(degree(rph, mode = "in")), mean(degree(rph, mode = "out")))
  
  ###DOLR Strategies###
  #Deciding on strategies
  belieftrust <- rcount[cc,2]
  
  #Number of trusting
  trustnode <- floor(belieftrust*length(focstrat))
  skepnode <- length(focstrat) - trustnode
  
  #Deciding on Focus
  if (trustnode > 1){
    x <- min(sample(c(1:(trustnode-1)), 1),min(trustnode,skepnode))
  }
  
  if(trustnode <=1){
    x <- 1
  }
  
  #Focus 1 - DOLR and Strategy
  focusnodes1 <- max(trustnode,skepnode) + x
  
  #Focus 2 - DOLR and Strategy
  focusnodes2 <- min(trustnode,x-1) + min(skepnode,x-1) + 1
  
  #Creating Strategy vectors
  focivec <- c(focusnodes1,focusnodes1,focusnodes2,focusnodes2,focusnodes1,focusnodes2)
  mechvec <- c(mechanism[1],mechanism[2],mechanism[1],mechanism[2],mechanism[1],mechanism[2])
  
  #Rinott#
  rinot <- vector(mode = "numeric", length = 10)
  skepspread <- vector(mode = "numeric", length = 10)
  for (l in 1:10) {
    k <- 6
    alpha <- 0.05
    n0 <- 100
    h <- Rinotth(k, n0, 1-alpha, 0.99, 10000)$UCB
    Ybar <- NULL
    Vars <- NULL
    Ns <- NULL
    delta <- 10
    for (s in 1:k) {
      if(s<=4){
        typevec <- sample(nodetype[[s]], nodenum, replace = TRUE)
        for (n in 1:n0) {
          eom <- 1
          invec <- eom
          opvec <- rep(3,nodenum)
          foclist[n,] <- random_gsimul_graph(eom,invec,rph,adjrph,nodetype[[s]],typevec,nodenum,mechvec[s],sample(focstrat,focivec[s]))
        }
        Y <- foclist[,2]
        S2 <- var(Y)
        N <- ceiling(h^2*S2/delta^2)
        if (N > n0){
          foxlist <- matrix(data = 0, nrow = N-n0, ncol = 3)
          for (n in 1:(N-n0)) {
            eom <- 1
            invec <- eom
            opvec <- rep(3,nodenum)
            foxlist[n,] <- random_gsimul_graph(eom,invec,rph,adjrph,nodetype[[s]],typevec,nodenum,mechvec[s],sample(focstrat,focivec[s]))
          }
          Y <- c(Y, foxlist[,2])
        }
        Ybar <- c(Ybar, mean(Y))
        Vars <- c(Vars, S2)
        N <- max(N, n0)
        Ns <- c(Ns, N)
      }
      else{
        typevec <- sample(nodetype[[s]], nodenum, replace = TRUE)
        for (n in 1:n0) {
          eom <- 1
          invec <- eom
          opvec <- rep(3,nodenum)
          foclist[n,] <- random_gsimul_graph(eom,invec,rph,adjrph,nodetype[[s]],typevec,nodenum,mechvec[s],sample(focstrat,focivec[s]))
        }
        Y <- foclist[,3] ##Collecting Broadcast results
        S2 <- var(Y)
        N <- ceiling((h^2)*S2/delta^2)
        if (N > n0){
          foxlist <- matrix(data = 0, nrow = N-n0, ncol = 3)
          for (n in 1:(N-n0)) {
            eom <- 1
            invec <- eom
            opvec <- rep(3,nodenum)
            foxlist[n,] <- random_gsimul_graph(eom,invec,rph,adjrph,nodetype[[s]],typevec,nodenum,mechvec[s],sample(focstrat,focivec[s]))
          }
          Y <- c(Y, foxlist[,3])
        }
        Ybar <- c(Ybar, mean(Y))
        Vars <- c(Vars, S2)
        N <- max(N, n0)
        Ns <- c(Ns, N)
      }
    }
    sol <- list(Best = which.max(Ybar), Ybar = Ybar, Var = Vars, N = Ns)
    rinot[l] <- sol$Best
    skepspread[l] <- length(which(typevec == "S" | typevec == "SU" | typevec == "SP"))
  }
  
  
  ##NSGS###
  # implements NSGS
  # k = number of systems
  # 1-alpha = desired PCS
  # n0 = first-stage sample size
  # delta = indifference-zone parameter
  # note: uses 99% UCB for Rinott's h
  h <- Rinotth(k, n0, 1-alpha/2, 0.99, 10000)$UCB
  # first apply subset selection
  Yall <- NULL
  nsgs <- vector(mode = "numeric", length = 10)
  skepspread2 <- vector(mode = "numeric", length = 10)
  for (l in 1:10) {
    k <- 6
    alpha <- 0.05
    n0 <- 100
    h <- Rinotth(k, n0, 1-alpha, 0.99, 10000)$UCB
    Ybar <- NULL
    Vars <- NULL
    Yall <- NULL
    Ns <- NULL
    delta <- 10
    
    for (s in 1:k){
      typevec <- sample(nodetype[[s]], nodenum, replace = TRUE)
      for (n in 1:n0) {
        eom <- 1
        invec <- eom
        opvec <- rep(3,nodenum)
        foclist[n,] <- random_gsimul_graph(eom,invec,rph,adjrph,nodetype[[s]],typevec,nodenum,mechvec[s],sample(focstrat,focivec[s]))
      }
      if(s <=4){
        Yall <- cbind(Yall, foclist[,2])
      }
      else{
        Yall <- cbind(Yall, foclist[,3])
      }
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
    if(length(Survivors) > 1){
      for (x in Survivors){
        N <- ceiling(h^2*S2[x]/delta^2)
        Y <- Yall[,x]
        if(N > n0){
          typevec <- sample(nodetype[[s]], nodenum, replace = TRUE)
          for (n in 1:n0) {
            eom <- 1
            invec <- eom
            opvec <- rep(3,nodenum)
            foclist[n,] <- random_gsimul_graph(eom,invec,rph,adjrph,nodetype[[s]],typevec,nodenum,mechvec[s],sample(focstrat,focivec[s]))
          }
          if(s <=4){
            Y <- c(Y, foclist[,2])
          }
          else{
            Y <- c(Y, foclist[,3])
          }
        }
        Ybars <- c(Ybars, mean(Y))
        N <- max(N, n0)
        Ns <- c(Ns, N)
      }
      Best = Survivors[which.max(Ybars)]
    }
    if(length(Survivors) <= 1){
      Best <- Survivors
      Ybars <- Ybar[Survivors]
      Ns <- n0
    }
    
    soln <- list(Best = Best, Survivors = Survivors, Ybar = Ybars, N = Ns)
    nsgs[l] <- soln$Best
    skepspread2[l] <- length(which(typevec == "S" | typevec == "SU" | typevec == "SP"))
  }
  
  #Capturing results
  broadper <- (sum(rinot>=5) + sum(nsgs>=5))*5
  focper <- (20-(sum(rinot>=5) + sum(nsgs>=5)))*5
  resout[cc,] <- c(rcount[cc,1],rcount[cc,2],broadper,focper,mean(skepspread),mean(skepspread2),focusnodes1,focusnodes2,props) 
  colnames(resout) <- c("Edge-likelihood", "Belief", "Broadcastpercent", "Focuspercent", "PercentSkepticalOverall1", "PercentSkepticalOverall2",  "Focus1", "Focus2", "mean_dist", "mean_in_deg", "mean_out_deg")
}




