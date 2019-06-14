######

# The code is divided in two parts:
# 1. The two main functions lead_abs() and add_abstraction()
# 2. The helping functions

######


library(pcalg)


#### THE TWO MAIN FUNCTIONS ####


## Leading function ("the leading part")
## takes as input the pdag and all the groups (all Zs)
## loops over the groups to decide which are pairs (i.e. which have connected edges)
## decides the order of directing edges between pairs
## first, the pairs are directed that have 1 or more directed edges
## last, all pairs that have no directed edges are directed
## all pairs are sent to the add_abstraction function to direct the edges
## output: the new pdag, with all directed edges in it
lead_abs <- function(pdag, abs_groups){
  original_pdag <- pdag
  undirected_list <- list()
  for(i in 1:length(abs_groups)){
    cat("\nworking on abs_group ", i)
    for(z in i:length(abs_groups)){
      if((z+1) <= length(abs_groups)){
        absList <- find_edge_between_groups(pdag, abs_groups[[i]], abs_groups[[z+1]])
        if(length(absList$abs_group1) > 0){
          cat("\nfound connecting edges")
          group1 <- noquote(gsub("[^0-9]","", absList$abs_group1))
          group2 <- noquote(gsub("[^0-9]","", absList$abs_group2))
          dir_list <- check_directed_edges(pdag, group1, group2)
          if(!is.null(dir_list)){
            cat("\ndirected edges found between groups")
            pdag <- add_abstraction(pdag, absList$abs_group1, absList$abs_group2)
          }
          else {
            cat("\nno directed edges found between groups")
            undirected_list <- append(undirected_list, absList)
          }
        }
      }
    }
  }
  if(length(undirected_list) > 0){
    cat("\nstart with undirected list")
    print(length(undirected_list))
    pdag_changed <- FALSE
    repeat{
      for(i in seq(1, length(undirected_list), by=2)){
        new_pdag <- add_abstraction(pdag, undirected_list[i]$abs_group1, undirected_list[i+1]$abs_group2)
        if(all(new_pdag != "both_directions_possible")){
          pdag_changed <- TRUE
          pdag <- new_pdag
          undirected_list[i] <- NULL
          undirected_list[i] <- NULL
        }
      }
      if(length(undirected_list) == 0 || pdag_changed == FALSE)
        break
      else
        pdag_changed <- FALSE
    }
  }
  # check for cycles on the highlevel!
  cycle_check <- check_highlevel_cycles(pdag, abs_groups)
  if(cycle_check){
    cat("\nERROR: directing the edges creates a cycle on the high-level. The original pdag is returned")
    return(original_pdag)
  }
  else {
    cat("\nno cycles found on the high-level. Directing the edges as such has been approved.")
    return(pdag)
  }
}


## Directing leading function ("the directing part")
## directs the edges between a pair of clusters
## there are 3 options when an abstraction is added:
## 1) there is no directed edge in the abstraction. This means both directions must be tried
## 2) there is 1 directed edge in the abstraction. This means all edges must be pointed that way
## 3) there are more than 1 directed edges in the abstraction. This means that first, it must
## be checked whether they all point in the same direction. Then, it must be checked if there
## are undirected edges left. If so, they must point in the same direction as well.
## input: a pair of clusters (abs_group1 and abs_group2)
## output: the updated pdag with directed edges

add_abstraction <- function(pdag, abs_group1, abs_group2){
  abs_group1 <- as.numeric(noquote(gsub("[^0-9]","", abs_group1)))
  abs_group2 <- as.numeric(noquote(gsub("[^0-9]","", abs_group2)))

  dir_list <- check_directed_edges(pdag, abs_group1, abs_group2)
  v_structures <- number_v_nodes(pdag)

  if (length(dir_list) == 0){
    #this means there are no directed edges between the pair
    #we try both directions of the edges

    ## pointing the one way
    pdag1 <- change_direction(pdag, abs_group1, abs_group2)
    dir1 <- is_direction_possible(pdag1, v_structures)

    ## pointing the other way
    pdag2 <- change_direction(pdag, abs_group2, abs_group1)
    dir2 <- is_direction_possible(pdag2, v_structures)

    if(dir1 && dir2){
      #this means that both directions are possible
      #we return the message "both_directions_possible"
      return("both_directions_possible")
    }
    if(dir1 && !dir2){
      #this means that the direction abs_group1 --> abs_group2 is the only possible direction
      pdag <- apply_mec_rules(pdag1)
      return(pdag)
    }
    if(!dir1 && dir2){
      #this means that the direction abs_group1 <- abs_group2 is the only possible direction
      pdag <- apply_mec_rules(pdag2)
      return(pdag)
    }
    if(!dir1 && !dir2){
      #this means that both directions are not possible
      #the original pdag is returned, without any directed edges
      return(pdag)
    }

  }
  else if (length(dir_list) == 2){
    #this means that there is one directed edge between the pair of clusters
    #we convert the other edges the same way
    group1_check <- check_direction(pdag, dir_list, abs_group1)
    if(!(group1_check)){
      q <- abs_group1
      abs_group1 <- abs_group2
      abs_group2 <- q
    }
    pdag <- change_direction(pdag, abs_group1, abs_group2)
    check <- is_direction_possible(pdag, v_structures)
    if(check){
      pdag <- apply_mec_rules(pdag)
      return(pdag)
    }
    else{
      cat("ERROR: directing the edges is not possible. Check if something is wrong with the graph.")
      return(pdag)
    }

  }
  else {
    #this means that multiple directed edges are found between the pair of clusters
    #we need to check if they are directed the same way, and direct the other edges
    leftrightlist <- list_dir_edges(dir_list)

    result <- match_dir_edges(leftrightlist$left, leftrightlist$right, abs_group1)

    if(result == "error"){
      cat("ERROR: the edges are not directed the same way. An abstraction is impossible!")
      return(pdag)
    }
    else if(result == "left"){
      #the directed edges point the same way
      for(i in abs_group1){
        if (!(i %in% leftrightlist$left)){
          #not all edges are directed yet
          group1_check <- check_direction(pdag, dir_list, abs_group1)
          if(!(group1_check)){
            q <- abs_group1
            abs_group1 <- abs_group2
            abs_group2 <- q
          }
          pdag <- change_direction(pdag, abs_group1, abs_group2)
          check <- is_direction_possible(pdag, v_structures)
          if(check){
            pdag <- apply_mec_rules(pdag)
            return(pdag)
          }
          else{
            cat("ERROR: directing the edges is not possible. Check if something is wrong with the graph.")
            return(pdag)
          }
        }
        else{
          #all edges are directed already. No extra edges need to be directed.
          return(pdag)
        }
      }
    }
    else if(result == "right"){
      #the directed edges point the same way
      for(i in abs_group1){
        if (!(i %in% leftrightlist$right)){
          #not all edges are directed yet
          group1_check <- check_direction(pdag, dir_list, abs_group1)
          if(!(group1_check)){
            q <- abs_group1
            abs_group1 <- abs_group2
            abs_group2 <- q
          }
          pdag <- change_direction(pdag, abs_group1, abs_group2)
          check <- is_direction_possible(pdag, v_structures)
          if(check){
            pdag <- apply_mec_rules(pdag)
            return(pdag)
          }
          else{
            cat("ERROR: directing the edges is not possible. Check if something is wrong with the graph.")
            return(pdag)
          }
        }
        else{
          #all edges are directed already. No extra edges need to be directed.
          return(pdag)
        }
      }
    }
  }
}


#### HELP FUNCTIONS ####


## NOTE: code from pcalg package
## function to produce the pdag from the skeleton
## the pdag is a tabular representation of the mec
## input parameter: skeleton (from pcalg)
## output: pdag of mec
find_pattern <- function(skel){
  res <- skel
  g <- as(skel, "matrix")
  p <- as.numeric(dim(g)[1])
  pdag <- g
  ind <- which(g == 1, arr.ind = TRUE)
  for (i in seq_len(nrow(ind))) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    allz <- setdiff(which(g[y, ] == 1), x)
    for (z in allz) {
      if (g[x, z] == 0 && !(y %in% skel@sepset[[x]][[z]] ||
                            y %in% skel@sepset[[z]][[x]])) {
        pdag[x,y] <- pdag[z,y] <- 1
        pdag[y,x] <- pdag[y,z] <- 0
      }
    }
  }
  return(pdag)
}


## function to find the edges between two groups of variables
## input: the pdag, the first group and the second group
## output: two lists in one list absList:
## absList$abs_group1 = in order, the variables that are connected to group_2
## absList$abs_group2 = in order, the variables that are connected to group_1
## so both lists are in order and define which variable is connected to which
find_edge_between_groups <- function(pdag, group1, group2){
  t = 0
  abs_group1 <- list()
  abs_group2 <- list()
  for(i in group1){
    conn_list = find_conn_edges(pdag, i)
    for(z in conn_list){
      if(z %in% group2){
        t = t+1
        abs_group1[t] <- i
        abs_group2[t] <- z
      }
    }
  }
  absList <- list("abs_group1" = abs_group1, "abs_group2" = abs_group2)
  return(absList)
}


## function to find the variables connected to a certain variable
## input: the pdag, and the variable for which you want to find the connected edges
## for example: x2 is connected to x1, x4, x5
## output: list of the connected edges
find_conn_edges <- function(pdag, row){
  conn_list = list()
  t = 0
  z = 0
  for(i in pdag[row,]){
    t = t+1
    if(i == 1){
      z = z+1
      conn_list[z] = colnames(pdag)[t]
    }
  }
  return(conn_list)
}


## function to check for directed edges between two groups of variables
## input: the two groups of variables
## output: the list of directed edges (if any) in one list
## if no directed edges are found, return is NULL
check_directed_edges <- function(pdag, abs_group1, abs_group2) {
  ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
  dir_edge <- FALSE
  ind_no <- list()
  t <- 1

  for (i in seq_len(nrow(ind))){
    if(ind[i,1] %in% abs_group1){
      if(ind[i,2] %in% abs_group2){
        dir_edge = TRUE
        ind_no[[t]] <- ind[i,1]
        t <- t+1
        ind_no[[t]] <- ind[i,2]
        t <- t+1
      }
    }
    else if (ind[i,2] %in% abs_group1){
      if(ind[i,1] %in% abs_group2){
        dir_edge = TRUE
        ind_no[[t]] <- ind[i,1]
        t <- t+1
        ind_no[[t]] <- ind[i,2]
        t <- t+1
      }
    }
  }

  if(dir_edge){
    return(ind_no)
  }
  else{
    return(NULL)
  }
}


## function to check the number of v-structures in a graph
## input: the pdag
## output: the number of v-structures
## this is useful for comparing the first MEC with the second MEC with abstraction:
## if the number of v-structures is not equal, the given abstraction is not valid
number_v_nodes <- function(pdag) {
  ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
  v_nodes <- ind[duplicated(ind[,2]),2]
  no_v_nodes <- length(v_nodes)
  return(no_v_nodes)
}


## function that changes all undirected edges into directed edges
## the edges are directed as abs_group1 --> abs_group2
## input: the pdag and the two abstraction groups
## output: the directed pdag
change_direction <- function(pdag, abs_group1, abs_group2){
  for(i in seq_along(abs_group1)){
    a <- abs_group1[i]
    b <- abs_group2[i]

    if(pdag[a,b] == 1 & pdag[b,a] == 1){
      pdag[b,a] <- 0
    }
  }
  return(pdag)
}


## this function checks whether a certain direction of the edges is possible
## it first check whether the direction alone causes any more v-structures
## then, it applies the MEC-Rules of the pcalg algorithm, and checks whether
## this causes any more v-structures
## if both don't cause more v-structures, the direction is allowed
## input: the pdag and the number of v-structures in the original MEC
## output: TRUE if direction is allowed, FALSE if direction is not allowed
is_direction_possible <- function(pdag, v_structures){
  possible_direction <- TRUE
  pdag_v <- number_v_nodes(pdag)
  if(v_structures != pdag_v){
    possible_direction <- FALSE
  }
  if(possible_direction){
    ## point other edges according to MEC-rules
    pdag <- apply_mec_rules(pdag)
    pdag_v <- number_v_nodes(pdag)
    if(v_structures != pdag_v){
      possible_direction <- FALSE
    }
  }
  return(possible_direction)
}


## NOTE: code from pcalg pacakge
## function that applies the 3 rules of PC algorithm
## input parameter: pdag
## output: alternated pdag according to the rules
## this function checks if edges must be directed according to the rules of mec
apply_mec_rules <- function(pdag){
  g <- as(skel, "matrix")
  p <- as.numeric(dim(g)[1])

  old_pdag <- matrix(0, p, p)

  ## Rule 1
  while (!all(old_pdag == pdag)) {
    old_pdag <- pdag
    ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      indC <- which((pdag[b, ] ==1 & pdag[, b] ==1) &
                      (pdag[a, ] == 0 & pdag[, a] ==0))
      if (length(indC) > 0) {
        pdag[b, indC] <- 1
        pdag[indC, b] <- 0
      }
    }

    ## Rule 2
    ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      indC <- which((pdag[a, ] == 1 & pdag[, a] == 0) &
                      (pdag[, b] == 1 & pdag[b, ] == 0))
      if (length(indC) > 0) {
        pdag[a, b] <- 1
        pdag[b, a] <- 0
      }
    }

    ## Rule 3
    ind <- which((pdag == 1 & t(pdag) ==1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      indC <- which((pdag[a, ] == 1 & pdag[, a] == 1) &
                      (pdag[, b] == 1 & pdag[b, ] == 0))
      if (length(indC) >= 2) {
        g2 <- pdag[indC, indC]
        if (length(g2) <= 1) {
          g2 <- 0
        }
        else {
          diag(g2) <- rep(1, length(indC))
        }
        if (any(g2 == 0)) {
          pdag[a, b] <- 1
          pdag[b, a] <- 0
        }
      }
    }
  }
  return(pdag)
}


## function to check the direction of the variables in dir_list
## assumptions: the variables in the abs_groups are directed equally,
## and there is at least 1 directed edge in abs_group1
## (this has been checked before the function is called)
## input: the pdag, the dir_list (from check_directed_edges function) and abs_group1
## if the direction is abs_group1 --> abs_group2, return group1 = TRUE
## if the direction is abs_group1 <-- abs_group2, return group1 = FALSE
check_direction <- function(pdag, dir_list, abs_group1){
  a <- dir_list[[1]]
  b <- dir_list[[2]]
  if(pdag[a,b] == 1 & pdag[b,a] == 0){
    #this means the edge is directed as a -> b
    if (a %in% abs_group1){
      return(group1 <- TRUE)
    }
    else{
      return(group1 <- FALSE)
    }
  }
  if(pdag[b,a] == 1 & pdag[a,b] == 0){
    #this means the edge is directed as a <- b
  }
  if (b %in% abs_group1){
    return(group1 <- FALSE)
  }
  else{
    return(group1 <- TRUE)
  }
}


## function to put the results of check_directed_edges function into two lists:
## one list with all the nodes found on the left side, one with nodes on the right side
## the functions thus splits up the list found in check_directed_edges
## input: dir_list (from check_directed_edges function)
## output: two lists: left_list and right_list
list_dir_edges <- function(dir_list){
  left_list <- list()
  right_list <- list()
  l <- 1
  r <- 1
  for (i in seq(1, length(dir_list), 2)){
    left_list[l] <- dir_list[[i]]
    l <- l+1
  }
  for (i in seq(2, length(dir_list), 2)){
    right_list[r] <- dir_list[[i]]
    r <- r+1
  }
  leftright_list <- list("left" = left_list, "right" = right_list)
  return(leftright_list)
}


## function to check if variables in leftlist are in either absgroup1 or absgroup2
## this is necessary, otherwise there has been made a mistake, in the given absgroups
## input: leftlist and rightlist (from list_dir_edges function), and abs_group1
## if there is a mistake, then leftlist is both in absgroup1 and absgroup2
## the function then returns "error"
## else, the function returns "left" when leftlist is in absgroup1
##        this means an abstraction of absgroup1 --> absgroup2
## the function return "right" when leftlist is in absgroup2
##        this means an abstraction of absgroup1 <-- absgroup2
match_dir_edges <- function(leftlist, rightlist, absgroup1){
  result <- "error"
  group1 <- FALSE
  group2 <- FALSE
  for (i in leftlist){
    if(i %in% absgroup1){
      group1 <- TRUE
    }
    else{
      group2 <- TRUE
    }
  }

  if(group1 && group2){
    return(result)
  }
  if(group1){
    for(i in rightlist){
      if(i %in% absgroup1){
        group1 <- FALSE
      }
      else{
        group2 <- TRUE
      }
    }
    if(group1 && group2){
      ##in this case, leftlist is in absgroup1, rightlist is in absgroup2
      result <- "left"
      return(result)
    }
  }
  else if(group2){
    for(i in rightlist){
      if(i %in% absgroup1){
        group1 <- TRUE
      }
      else{
        group2 <- FALSE
      }
    }
    if(group1 && group2){
      ##in this case, leftlist is in absgroup2, rightlist is in absgroup1
      result <- "right"
      return(result)
    }
    else{
      return(result)
    }
  }
}


## function to check whether there are cycles on the high-level
## input: low-level pdag and abs_groups
## output: value of cycle_check:
## if TRUE: a cycle exists on the high-level
## if FALSE: no cycle exists on the high-level
check_highlevel_cycles <- function(pdag, abs_groups){
  highlevel_pdag <- construct_highlevel(pdag, abs_groups)
  cycle_check <- find_cycle(highlevel_pdag)
  return(cycle_check)
}


## function to construct the high-level pdag
## that is: the table in which the edges between the clusters are defined
## input: the low-level pdag and abs_groups
## output: the high-level pdag
construct_highlevel <- function(pdag, abs_groups){
  p <- length(abs_groups)
  highlevel_pdag <- matrix(0, p, p)
  for(i in 1:length(abs_groups)){
    for(z in 1:length(abs_groups)){
      absList <- find_edge_between_groups(pdag, abs_groups[[i]], abs_groups[[z]])
      if(length(absList$abs_group1) > 0){
        group1 <- noquote(gsub("[^0-9]","", absList$abs_group1))
        group2 <- noquote(gsub("[^0-9]","", absList$abs_group2))
        dir_list <- check_directed_edges(pdag, group1, group2)
        if(!is.null(dir_list)){
          outgoing_edge <- check_direction(pdag, dir_list, group1)
          if(outgoing_edge){
            highlevel_pdag[i, z] <- 1
          }
        }
      }
    }
  }
  return(highlevel_pdag)
}


## function to check whether there is a cycle in the model
## input: pdag
## output: cycle value which is TRUE if a cycle is found, FALSE otherwise
find_cycle <- function(pdag){
  cycle <- FALSE
  ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
  for (i in seq_len(nrow(ind))) {
    a <- ind[i, 1]
    b <- ind[i, 2]
    c <- which((ind[,1]==b))

    while(length(c) > 0){
      b <- ind[c,2]
      if(a %in% b){
        break
      }
      c <- which((ind[,1]==b))
    }
    if(a %in% b){
      cycle <- TRUE
    }
  }
  return(cycle)
}
