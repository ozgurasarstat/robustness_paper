
 sim_mod1 <- function(nsubj = 100, 
                      ntime = 5,
                      alpha = matrix(c(1, 0.5), ncol = 1), 
                      Sigma = matrix(c(1, 0.5, 0.5, 1), ncol = 2), 
                      sigmasq = 4,
                      dof = 3.5
                      ){
   
   id <- rep(1:nsubj, each = ntime)
   
   t <- rep(0:(ntime - 1), nsubj)
   
   x <- dmat <- cbind(1, t)
   
   dlist <- dmat %>% data.frame %>% split(id) %>% lapply(as.matrix)
   
   d <- do.call(magic::adiag, dlist)
   
   Ustar_mat <- rmvnorm(nsubj, rep(0, 2), Sigma)
   Zstar <- rnorm(nsubj * ntime, 0, sqrt(sigmasq))
   
   V_indv <- rgamma(nsubj, shape = dof/2, rate = dof/2)
   V_ext  <- rep(V_indv, each = ntime)
   
   U_mat <- Ustar_mat
   U_mat[, 1] <- U_mat[, 1]/sqrt(V_indv)
   U_mat[, 2] <- U_mat[, 2]/sqrt(V_indv)
   
   U <- as.numeric(t(U_mat))
   
   Z <- Zstar/sqrt(V_ext)
   
   Y <- x %*% alpha + d %*% U + Z
   
   data <- data.frame(id = id, y = Y, t = t)
   mget(c("data", "U_mat"))  
   
 }

 