##Function to update node opinions given an input##
opinion_update <- function(invec,eom,adjmat,typevec,opvec,act,gg,nact) { # create a function with the name my_function
  #eom is the starting node
  #Non-eom nodes will be using loops
  #This function uses three matrices/vectors: First, adjacency matrix itself. Second, it uses a typevecrix where each entry is a character label according to type of the node
  #Third, it uses an opinion vector which stores the opinion of each node, post update.
  if(sum(invec) == eom){
    immn <- neighbors(gg, invec)
    ##Starting loop for immediate neighbors of eom
    if(act == "B"){
      for (j in 1:length(immn)) {
        if(typevec[immn[j]] == "S"){
          opvec[immn[j]] = 0
        }
        else{
          opvec[immn[j]] = 1
        }
      }
    }
    else{
      for (j in 1:length(nact)) {
        if(typevec[nact[j]] == "S"){
          opvec[nact[j]] = 0
        }
        else{
          opvec[nact[j]] = 1
        }
      }
      immn <- nact
    }
   
  }
  else{
    nodeunion <- list()
    for (j in 1:length(invec)) {
      vee <- neighbors(gg, invec[j])
      nodeunion <- c(nodeunion, vee)
    }
    nunion <- unique(unlist(nodeunion))
    
    if(any(nunion) == 1){
      nunion <- nunion[! nunion %in% c(1)]
    }
    
    ##Sticky or Loose opinion update -- this code snippet determines##
    upd_idx <- which(opvec<=1) ##Nodes that have opinions updated to 0 or 1 from previous rounds
    nunion <- nunion[! nunion %in% upd_idx] #Removing the nodes that have opinions
    if(length(nunion) > 0){
      eig_val_vec <- centr_eigen(gg, directed = TRUE)$vector[invec]
      up_node <- invec[which.max(eig_val_vec)]
      
      up_val <- opvec[up_node]
      opvec[nunion] <- up_val
      
      immn <- nunion
    }
    else{
      immn <- nunion
      return(return(list(opvec,immn,c("Finished"))))
    }
    
  }
  
  return(list(opvec,immn,c("In progress")))
  
} 