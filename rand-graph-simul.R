random_gsimul_graph <- function(eom,invec,gg,adjrph,nodetype,typevec,nodenum,mechanism,focusnodes){
  typevec <- typevec
  #Centrality Measures
  #Degree.Directed <- degree(rph)
  #Indegree <- degree(rph, mode="in")
  #Outdegree <- degree(rph, mode="out")
  
  #Eig <- eigen_centrality(rph)$vector
  
  #Hub <- hub_score(rph)$vector
  #Authority <- authority_score(rph)$vector
  
  #Closeness <- closeness(rph, mode = "all")
  
  #reach1 <- ego_size(rph, 1)/(vcount(rph)-1)
  
  #Betweenness <- betweenness(rph)
  
  #centralities <- cbind(mean(na.omit(Degree.Directed)), mean(na.omit(Eig)), mean(na.omit(Hub)), mean(na.omit(Authority)), mean(na.omit(Closeness)), mean(na.omit(reach1)), mean(na.omit(Betweenness)))
  
  #Homophily Measures
  #fhom <- coleman(rph, "label")
  #fhom2 <- ei(rph, typevec, directed = TRUE)
  
  ##Section for Focus##
  
  for (k in 1:10) {
    if(rand_opinion_update(invec,eom,adjrph,typevec,opvec,"F",rph,focusnodes,mechanism)[[3]] != "Finished"){
      
      linvec <- rand_opinion_update(invec,eom,adjrph,typevec,opvec,"F",rph,focusnodes,mechanism)[[2]]
      lopvec <- rand_opinion_update(invec,eom,adjrph,typevec,opvec,"F",rph,focusnodes,mechanism)[[1]]
      
      invec <- linvec
      opvec <- lopvec
    }
  }
  
  #Output from Focus
  fopvec <- opvec
  
  
  ##Section for Broadcast##
  eom <- 1
  invec <- eom
  opvec <- rep(3,nodenum)
  
  for (k in 1:10) {
    if(rand_opinion_update(invec,eom,adjrph,typevec,opvec,"B",rph,focusnodes,mechanism)[[3]] != "Finished"){
      
      linvec <- rand_opinion_update(invec,eom,adjrph,typevec,opvec,"B",rph,focusnodes,mechanism)[[2]]
      lopvec <- rand_opinion_update(invec,eom,adjrph,typevec,opvec,"B",rph,focusnodes,mechanism)[[1]]
      
      invec <- linvec
      opvec <- lopvec
    }
  }
  
  #Output from Broadcast
  bopvec <- opvec
  
  #Counting number of 1 in focus
  nfop <- sum(fopvec==1)
  
  #Counting number of 1 in broadcast
  nbop <- sum(bopvec==1)
  
  return(c(length(focusnodes),nfop,nbop))
  
}