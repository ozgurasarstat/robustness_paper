
 mod1_pred_res <- function(x, y, n_sim, idlist, time_var){
   
   nsubj <- length(x)
   
   out <- NULL
   
   for(i in 1:nsubj){
     
     out <- rbind(out, 
                  do.call(cbind, x[[1]]) %>% 
                    apply(1, function(x) c(mean(x), median(x), quantile(x, probs = c(0.025, 0.975)))) %>% t
                  )
     
   }
   
   out <- data.frame(idlist, time_var, y, out)
   names(out) <- c("id", "time", "observed", "mean", "median", "2.5%", "97.5%")
   rownames(out) <- NULL
   return(out)
   
 }
 