
 mod1_pred_res <- function(x, y, n_sim, idlist, time_var, probs){
   
   nsubj <- length(x)
   
   out <- NULL
   
   for(i in 1:nsubj){
     
     out <- rbind(out, 
                  do.call(cbind, x[[i]]) %>% 
                    apply(1, function(x) c(mean(x), median(x), quantile(x, probs = probs))) %>% t
                  )
     
   }
   
   out <- data.frame(idlist, time_var, y, out)
   names(out) <- c("id", "time", "observed", "mean", "median", "lower", "upper")
   rownames(out) <- NULL
   return(out)
   
 }
 