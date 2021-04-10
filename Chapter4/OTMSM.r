source('example_DD.r')  # load data of Dixit&Dill example 

#### MODEL INFERENCE ####
OT_inference <- function(x, y, pi, C, C_avg, g=.1,
                         max_iter=1000, convergence_threshold =10^-7){
  # x, y : n by m array of m initial and final distr.s
  # pi   : stationary distribution
  # C    : cost matrix (amino acid contacts broken)
  # C_avg: avg of observable
  # g    : learning rate
  
  ## Initialization
  n <- dim(x)[1]; m <- dim(x)[2]
  a <- 10
  l <- replicate(m, numeric(n)+1)
  r <- array(apply(x, 2, "/", pi), dim=c(n,m))
  
  # Convergences
  l_conv <- 0  
  r_conv <- 0
  a_conv <- 0
  
  # Data to be stored
  a_n <- numeric(max_iter)
  k_n <- array(dim=c(max_iter, n*n))
  l_n <- array(dim=c(max_iter, n*m))
  rho_n <- array(dim=c(max_iter, n))
  
  ## Helper functions
  sum_OT_LMs <- function(l,r){
    # sums the OT lagrange multipliers, i.e.:
    # computes lr = sum_i(l^i_a r^i_b + l^i_b r^i_a)
    lr <- 0
    for(j in 1:dim(r)[2]){
      lr <- lr + as.vector(l[,j])%o%as.vector(r[,j])
    }
    return(lr)
  }
  
  findRho <- function(W, pi, n, cth, max_i){
    # Finds Lagrange Multiplier rho to normalize row sums of k_ab
    # fix sum_b W_ab rho_b = pi_a/rho_a for one a at the time
    rho <- numeric(n)+.1
    rho_old <- numeric(n)+1
    i <- 0
    while(all(c(any(abs(rho_old-rho)>cth), i<max_i))){
      j <- i %% n
      if(j==0){rho_old<-rho}
      rho[j+1] <- pi[j+1]/sum(rho*W[j+1,])
      i <- i+1
    }
    return(rho)
  }
  
  ## Main loop
  for(iter in 1:max_iter){  
    # W_ab for given a and l:
    # compute W (eq. 17 in Dixit and Dill)
    W <- exp(-0.5*((C + t(C))/a + 
                     (sum_OT_LMs(l,r)+t(sum_OT_LMs(l,r)))))
    
    # Lagrange multiplier rho:
    rho <- findRho(W,pi,n,1e-15,1000)
    
    # Transition probabilities
    # k_ab= (rho_a rho_b)/pi_a W_ab 
    k <- array((rho/pi) %o% rho, dim=c(n,n)) * W
    
    l <- l+g*(apply(x, 2,'%*%', k) - y)
    a <- a-g*(sum(pi %*% (k*C)) - C_avg)
    k_n[iter,] <- as.vector(k); a_n[iter] <- a
    l_n[iter,] <- l; rho_n[iter,] <- rho
    
    #convergence
    if(iter > 20){
      if(sum(abs(l_n[iter,]-l_n[iter-1,])) <
         convergence_threshold && l_conv==0){
        print("l converged")
        l_conv <- 1
      }
      if(sum(abs(rho_n[iter,]-rho_n[iter-1,])) <
         convergence_threshold && r_conv==0){
        print("rho converged")
        r_conv <- 1
      }
      if(sum(abs(a_n[iter] - a_n[iter-1])) <
         convergence_threshold  && a_conv==0){
        print("alpha converged")
        a_conv <- 1
      }
      if(l_conv && r_conv && a_conv){
        print(iter)
        break
      }
    }
  }
  return(list(k=k, k_all=k_n, l=l_n, 
              alpha=a_n, rho=rho_n))
}

i  <- 1 # initial state of transient dynamics
x_ <- array(x[,i], dim=c(n,length(i))); y_ <-array(y[,i], dim=c(n,length(i)))
OT_inference(x_, y_, pi, C, C_avg, g=.1,
                         max_iter=1000, convergence_threshold =10^-7)

  
