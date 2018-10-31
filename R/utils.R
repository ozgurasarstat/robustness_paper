 covered <- function(true, lower, upper){
   ifelse(true > lower & true < upper, 1, 0) %>% sum %>% '/'(length(true))
 }
 
 interval_length <- function(lower, upper, fun = list(mean)){
    fun[[1]](upper - lower)
 }
 