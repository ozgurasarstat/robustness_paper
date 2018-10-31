mod1_pred <- function(object,
                      fixed,
                      random,
                      id_pred = NULL,
                      id_var,
                      time_var,
                      data,
                      type = "Smoothing",
                      controls = list(n_sim = 500,
                                      probs = c(0.025, 0.975))){

  
  if(length(controls) < 2){
    controls_full <- list(n_sim = 500,
                          probs = c(0.025, 0.975))
    for(i in 1:2){
      if(!(names(controls_full)[i] %in% names(controls))){
        controls[names(controls_full)[i]] <- controls_full[i] 
      }
    }
    
  }
  

  summary_object <- summary(object)

  fixed_est   <- matrix(summary_object$coef[, 1], ncol = 1)
  Sigma_est   <- summary_object$theta
  sigmasq_est <- summary_object$scale
  dof_est     <- summary_object$family[2] %>% as.character %>% str_extract("\\-*\\d+\\.*\\d*") %>% as.numeric

  #data <- data[order(data[, id_var], data[, time_var]), ]
  
  idlist <- data[, id_var]

  id_x <- cbind(idlist, model.matrix(fixed, data))
  id_y <- cbind(idlist, model.frame(fixed, data)[, 1])
  id_d <- cbind(idlist, model.matrix(random, data))

  if(is.null(id_pred)){
    id_pred <- unique(idlist)
  }

  #### BELOW IS FOR SMOOTHING

  U_pred <- Ystar_pred <- list()

  for(i in 1:length(id_pred)){

    x_i <- id_x[id_x[, 1] == id_pred[i], -1, drop = FALSE]
    y_i <- id_y[id_y[, 1] == id_pred[i], -1, drop = FALSE]
    d_i <- id_d[id_d[, 1] == id_pred[i], -1, drop = FALSE]

    n_i <- nrow(x_i)

    marg_resid_i  <- y_i - x_i %*% fixed_est
    deltasq_mid_i <- solve(d_i %*% Sigma_est %*% t(d_i) + sigmasq_est * diag(n_i))
    deltasq_i     <- t(marg_resid_i) %*% deltasq_mid_i %*% marg_resid_i

    V_giv_Y_shape <- (dof_est + n_i)/2
    V_giv_Y_rate  <- (dof_est + deltasq_i)/2

    U_giv_Y_and_V_mean_1 <- (Sigma_est %*% t(d_i)) %*% deltasq_mid_i %*% marg_resid_i
    U_giv_Y_and_V_var_1  <- Sigma_est - (Sigma_est %*% t(d_i)) %*% deltasq_mid_i %*% (d_i %*% Sigma_est)

    U_giv_Y_and_V_sim <- NULL
    Ystar_pred_i <- list()

    for(j in 1:controls$n_sim){
      V_sim_j <- rgamma(1, V_giv_Y_shape, V_giv_Y_rate)
      U_giv_Y_and_V_sim_j <- rmvnorm(1, U_giv_Y_and_V_mean_1, U_giv_Y_and_V_var_1/V_sim_j)

      Ystar_pred_i[[j]] <- x_i %*% fixed_est + d_i %*% t(U_giv_Y_and_V_sim_j)

      U_giv_Y_and_V_sim <- rbind(U_giv_Y_and_V_sim, U_giv_Y_and_V_sim_j)

      }

    U_pred[[i]] <- U_giv_Y_and_V_sim
    Ystar_pred[[i]] <- Ystar_pred_i

  }

  names(U_pred) <- names(Ystar_pred) <- id_pred
  
  pred.data <- mod1_pred_res(x = Ystar_pred, 
                             n_sim = controls$n_sim, y = id_y[, 2],
                             idlist = idlist, 
                             time_var = data[, time_var], 
                             probs = controls$probs)
  
  mget(c("U_pred", "Ystar_pred", "pred.data"))

}

