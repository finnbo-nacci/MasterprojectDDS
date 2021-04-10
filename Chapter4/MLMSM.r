source('example_DD.r')  # load data of Dixit&Dill example 

simulation <- function(init_state, K, length){
  # Input parameters
  # init_state: initial state
  # K: matrix of transition probabilities
  # length: number of time steps to simulate 
  
  # Initialize
  n <- dim(K)[1]  # number of states
  trajectory <- c(init_state)
  state <- init_state
  # Main loop
  for(time_step in 1:length){
    p <- K[state,] # transition probabilities given current state
    state <- sample.int(n, size=1 , prob=p) # sample a state according to p
    trajectory <- c(trajectory, state) # add the state to the trajectory
  }
  return(trajectory)
}


transition_index <- function(jump,n){
  from <- jump[1]; to <- jump[2]
  return((to-1)*n + from)
}


count_transitions <- function(traj, n){
  #count_matrix <- matrix(0, nrow=n, ncol=n)
  jumps <- cbind(traj[1:(length(traj)-1)], traj[2:(length(traj))])
  indices <- apply(jumps, 1, transition_index, n)
  counts <- sapply(indices, function(i,n){cm<-numeric(n**2); cm[i]<-1; return(cm)}, n)
  count_matrix <- matrix(rowSums(counts), nrow=n, ncol=n)
  return(count_matrix)
}


rev_ML <- function(counts){
  # initialize
  n <- dim(counts)[1]
  x <- counts + t(counts)
  rx <- rowSums(x)
  rc <- rowSums(counts)
  # loop
  for(iter in 1:10){
    for(i in 1:n){
      x[i,i] <- counts[i,i]*(rx[i]-x[i,i])/(rc[i]-counts[i,i])
    }
    rx <- rowSums(x)
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        a <- rc[i]-counts[i,j]+rc[j]-counts[i,j]
        b <- rc[i]*(rx[j]-x[i,j])+rc[j]*(rx[i]-x[i,j]) - 
          (counts[i,j]+counts[j,i])*(rx[i]+rx[j]-2*x[i,j])
        d <- -(counts[i,j]+counts[j,i])*(rx[i]-x[i,j])*(rx[j]-x[i,j])
        x[i,j] <- x[j,i] <- (-b+sqrt(b**2 - 4*a*d))/(2*a)
        rx <- rowSums(x)
      }
    }
  }
  return(x/rx)
}


get_pi <- function(MSM){
  pi_MSM <- eigen(t(MSM))$vectors[,1]
  return(pi_MSM/sum(pi_MSM))
}


ML_est <- function(n, traj){
  # Estimates a ML MSM from a single trajectory
  ##  INPUT; n: number of states; traj: trajectory 
  ## OUTPUT; MSM_ML; MSM_ML_DB
  
  ### COUNT TRANSITIONS
  counts <- count_transitions(traj, n)  # construct countmatrix.

  ### ESTIMATION
  MSM_ML <- (counts/rowSums(counts))
  MSM_ML_DB <- rev_ML(counts)           # estimate time reversible ML MSM
  
  ### STATIONARY DISTRIBUTIONS
  if(!any(c(is.nan(MSM_ML),is.nan(MSM_ML_DB)))){
  pi_ML    <- get_pi(MSM_ML)
  pi_ML_DB <- get_pi(MSM_ML_DB)
  }
  else pi_ML <- pi_ML_DB <- NA
  
  #return(MSM_ML_DB)
  return(list(MSM_ML=MSM_ML, MSM_MLDB=MSM_ML_DB, CountMatrix=counts, 
              pi_ML=pi_ML, pi_ML_DB=pi_ML_DB))
}


MD_data_long <- simulation(1, K, 2000)
traj_lengths <- seq(1,length(MD_data_long),25)
models <- lapply(traj_lengths, function(traj_length){
                                return(ML_est(n, MD_data_long[1:traj_length]))})
MLMSMerrors <- sapply(models, function(model) return(sum(abs(model$MSM_ML-K))))
DBMSMerrors <- sapply(models, function(model) return(sum(abs(model$MSM_MLDB-K))))

