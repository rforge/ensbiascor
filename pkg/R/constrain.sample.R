

## Get kernel density percentiles:
# ---------------------------------
get.kernel.density.percentiles <- function(cur.obs, probs=c(0.05, 0.25, 0.5, 0.75, 0.95), n.sample=1000) {
  
  kernel.bandwidth = tryCatch(dpik(x=cur.obs,scalest="stdev"), error =
                                function(e) NA)  
  if (is.na(kernel.bandwidth)) {
    return(rep(NA, length(probs)))
  } else {
    kernel.bandwidth = dpik(x=cur.obs,scalest="stdev")
    kernel.sample = sapply(X= sample(x=cur.obs, size=n.sample, replace = T), FUN=function(x) rnorm(n=1, mean=x, sd=kernel.bandwidth))
    kernel.percentiles = quantile(x=kernel.sample, probs=probs)
    
    # kernel.percentiles[which(kernel.percentiles < 0)] <- 0
    # kernel.percentiles[which(kernel.percentiles > 1)] <- 1
    
    return(kernel.percentiles)  
  }
}


