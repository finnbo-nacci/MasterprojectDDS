#### EXAMPLE DATA DIXIT & DILL ####
# Data of transition probabilities and observable
trans_prob <- c(0.0123, 0.9537, 0.0169, 0.0171, 
                0.0477, 0.6155, 0.0881, 0.2487, 
                0.0042, 0.4406, 0.3433, 0.2119, 
                0.0002, 0.0663, 0.0113, 0.9222)
aa_contacts <- c( 0,10,20,50, 
                  10, 0,10,40, 
                  20,10, 0,30, 
                  50,40,30, 0)
n <- sqrt(length(trans_prob))
K <- matrix(trans_prob,  ncol=n, nrow=n, byrow=TRUE) 
C <- matrix(aa_contacts, ncol=n, nrow=n, byrow=TRUE)

# Stationary distribution i.e. dominant eigenvector
pi <- eigen(t(K))$vectors[,1]
pi <- pi/sum(pi)               # Normalize the vector

# Compute the average cost
C_avg <- sum(pi %*% (K*C))  # <M> in D&D

# transient dynamics info
xU <- xI1 <- xI2 <- numeric(n)
xU[1] <- 1; xI1[2] <- 1; xI2[3] <- 1
yU <- xU %*% K; yI1 <- xI1 %*% K; yI2 <- xI2 %*% K
x <- array(c(xU, xI1, xI2),dim=c(n,3))
y <- array(c(yU, yI1, yI2),dim=c(n,3))
