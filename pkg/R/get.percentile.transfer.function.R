## Sebastian Sippel
# 05.09.2016

#' @title A function to derive a percentile transfer function 
#' @export
#' @description Calculates a transfer function between percentile of observed / simulated distributions 
#' @usage get.percentile.transfer.function(obs.grid.cell, mod.grid.cell, n.sample=10000)
#' @param obs.grid.cell A numeric vector of observed values
#' @param mod.grid.cell A numeric vector of simulated values
#' @param n.sample Sample size
#' @details 
#' Calculates a transfer function between percentile of observed / simulated distributions. The transfer function is based on Hermite splines and a Gaussian kernel is fitted to observed data.
#' @return XYZ
#' @references Sippel, S., Otto, F. E. L., Forkel, M., Allen, M. R., Guillod, B. P., Heimann, M., Reichstein, M., Seneviratne, S. I., Thonicke, K. & Mahecha, M. D. (2016) A novel bias correction methodology for climate impact simulations. Earth System Dynamics, 7, 71-88. doi:10.5194/esd-7-71-2016.
#' @author Sebastian Sippel
get.percentile.transfer.function <- function(obs.grid.cell, mod.grid.cell, n.sample=10000) {
  
  # 1. estimate Gaussian Kernel and estimate percentiles:
  kernel.bandwidth = dpik(x=obs.grid.cell)
  kernel.sample = sapply(X= sample(x=obs.grid.cell, size=n.sample, replace = T), FUN=function(x) rnorm(n=1, mean=x, sd=kernel.bandwidth))
  kernel.percentiles = quantile(x=kernel.sample, probs=seq(0, 1, 0.001))
  
  # 4. Fit kernel over model distribution and determine percentiles:
  mod.kernel.bandwidth = dpik(x=mod.grid.cell)
  mod.sample = sort(sapply(X=sample(x=mod.grid.cell, size=10^5, replace=T), FUN=function(x) rnorm(n=1, mean = x, sd = mod.kernel.bandwidth)), decreasing = F)
  mod.sample.percentiles = quantile(mod.sample, probs=seq(0, 1, 0.0001))
  
  # 5. make splinefun from percentiles to model percentiles:
  test = sapply(X=kernel.percentiles, FUN=function(x) which(abs(x - mod.sample.percentiles) == min(abs(x - mod.sample.percentiles)))[1] / length(mod.sample.percentiles))
  
  # 6. return splines to determine ranks:
  return(splinefun(x = seq(0, 1, 0.001), y = test, method="monoH.FC"))
}

