rand_graph_sim_fun <- function(gg,adjmat){
  #Simulation specs
  nodetype <- c("SU", "SP", "TU", "TP")
  nodenum <- 20
  mechanism <- c("dolr", "dolrp")
  opvec <- rep(3,nodenum)
  typevec <- sample(nodetype, nodenum, replace = TRUE)
  
  #graph generation
  rph <- as.directed(sample_grg(nodenum,0.8), mode = "random")
  
  adjrph <- as.matrix(as_adjacency_matrix(rph))
  
  #Deciding strategy
  focstrat <- which(adjrph[1,] == 1)
  
  if(length(focstrat)<=2){
    ridx <- sample(nrow(adjrph), 4)
    ridx <- ridx[! ridx %in% focstrat]
    adjrph[1,ridx] = 1
  }
  
  #Updating the graph after adding random edges
  rph <- graph_from_adjacency_matrix(adjrph, mode = "directed")
  
  ###DOLR Strategies###
  #Deciding on strategies
  belieftrust <- 0.5
  
  #Number of trusting
  trustnode <- floor(belieftrust*length(focstrat))
  skepnode <- length(focstrat) - trustnode
  
  #Focus 1 - DOLR and Strategy
  focusnodes1 <- max(trustnode,skepnode) + x
  
  #Focus 2 - DOLR and Strategy
  focusnodes2 <- min(trustnode,x-1) + min(skepnode,x-1) + 1
  
  #random_gsimul_graph(eom,invec,rph,adjrph,nodetype,nodenum,mechanism[2],focusnodes)
  
  #rand_opinion_update(invec,eom,adjrph,typevec,opvec,"B",rph,focusnodes,"dolr")
  
  #Output
  foclist <- matrix(data = 0, nrow = 100, ncol = 3)
  for (n in 1:100) {
    eom <- 1
    invec <- eom
    opvec <- rep(3,nodenum)
    foclist[n,] <- random_gsimul_graph(eom,invec,rph,adjrph,nodetype,nodenum,mechanism[2],focusnodes)
  }
  
}