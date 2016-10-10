## Sebastian Sippel
# 05.09.2016
# require(KernSmooth)



# Derive Kernel for illustration:
#' @export
#' @author Sebastian Sippel
get.kernel.percentiles <- function(data, kernel.quantiles = seq(0, 1, 0.001)) {
  kernel.bw = KernSmooth::dpik(x=data)
  kernel.sample = sapply(X= sample(x=data, size=10000, replace = T), FUN=function(x) rnorm(n=1, mean=x, sd=kernel.bw))
  kernel.percentiles = quantile(x=kernel.sample, probs= kernel.quantiles)
  return(kernel.percentiles)
}

# kernel.sample = sapply(X= sample(x=Tair_Europe_ERAI.JJA.area.original, size=1000, replace = T), FUN=function(x) rnorm(n=1, mean=x, sd=dpik(x=Tair_Europe_ERAI.JJA.area.original)))
# kernel.quantile.seq = seq(0, 1, 0.001)
# kernel.percentiles = quantile(x=kernel.sample, probs= kernel.quantile.seq)
# model.percentile = (which(abs(sort(unlist(HadRM3P.field16.JJA.means) - 273.15, decreasing = F) - median(kernel.percentiles)) == min(abs(sort(unlist(HadRM3P.field16.JJA.means) - 273.15, decreasing = F) - median(kernel.percentiles)))) / length(unlist(HadRM3P.field16.JJA.means)))[1]




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
  kernel.bandwidth = KernSmooth::dpik(x=obs.grid.cell)
  kernel.sample = sapply(X= sample(x=obs.grid.cell, size=n.sample, replace = T), FUN=function(x) rnorm(n=1, mean=x, sd=kernel.bandwidth))
  kernel.percentiles = quantile(x=kernel.sample, probs=seq(0, 1, 0.001))
  
  # 4. Fit kernel over model distribution and determine percentiles:
  mod.kernel.bandwidth = KernSmooth::dpik(x=mod.grid.cell)
  mod.sample = sort(sapply(X=sample(x=mod.grid.cell, size=10^5, replace=T), FUN=function(x) rnorm(n=1, mean = x, sd = mod.kernel.bandwidth)), decreasing = F)
  mod.sample.percentiles = quantile(mod.sample, probs=seq(0, 1, 0.0001))
  
  # 5. make splinefun from percentiles to model percentiles:
  test = sapply(X=kernel.percentiles, FUN=function(x) which(abs(x - mod.sample.percentiles) == min(abs(x - mod.sample.percentiles)))[1] / length(mod.sample.percentiles))
  
  # 6. return splines to determine ranks:
  return(splinefun(x = seq(0, 1, 0.001), y = test, method="monoH.FC"))
}


## ORIGINAL UNTOUCHED FUNCTION (cf. ESD paper):
# Sebastian Sippel
# 04.10.2016
# resample.ensemble  <- function(mod.grid.cell, transfer.function, sample.size = 1000) {
# mod.sorted.idx = ceiling(transfer.function(x=runif(n=sample.size)) * length(mod.grid.cell))
#  hadrm3p.idx = sapply(X=sort(mod.grid.cell, decreasing=F)[mod.sorted.idx], FUN=function(idx) which(idx == mod.grid.cell)[1])
# return(hadrm3p.idx)
# }




# Resampling ensembles: NEW FUNCTION
# ---------------------------------------------
# Stopping criterion / parameters:
# ---
# i. How often can any ensemble member be chosen? OK -> bootstrap.resampling
# ii. Which tail of the ensemble should be resampled? -> OK quantile.range
# iii. "search.radius" for each ensemble member. -> OK search.radius

#' @title A function to resample a large ensemble from a spline-based transfer function 
#' @export
#' @description Resamples a large ensemble from a spline-based transfer function
#' @usage resample.ensemble(obs.grid.cell, mod.grid.cell, n.sample=10000)
#' @param obs.grid.cell A numeric vector of observed values
#' @param mod.grid.cell A numeric vector of simulated values
#' @param n.sample Sample size
#' @details 
#' Calculates a transfer function between percentile of observed / simulated distributions. The transfer function is based on Hermite splines and a Gaussian kernel is fitted to observed data.
#' @return XYZ
#' @references Sippel, S., Otto, F. E. L., Forkel, M., Allen, M. R., Guillod, B. P., Heimann, M., Reichstein, M., Seneviratne, S. I., Thonicke, K. & Mahecha, M. D. (2016) A novel bias correction methodology for climate impact simulations. Earth System Dynamics, 7, 71-88. doi:10.5194/esd-7-71-2016.
#' @author Sebastian Sippel
resample.ensemble  <- function(mod.grid.cell, transfer.function, 
                               replacement = T, sample.size = 1500, 
                               quantile.range = c(0, 1), bootstrap.resampling = NULL, 
                               search.radius = NA) {
    
  bootstrap.resampling.samples = list()
  
  if (is.null(bootstrap.resampling)) bootstrap.resampling = 1
  
  for (r in 1:bootstrap.resampling) {
    if (replacement == T) {
      hadrm3p.idx = resample.ensemble_replacement(mod.grid.cell=mod.grid.cell, transfer.function=transfer.function, sample.size=sample.size, quantile.range=quantile.range)
    } else {
      hadrm3p.idx = resample.ensemble_no_replacement(mod.grid.cell=mod.grid.cell, transfer.function=transfer.function, sample.size=sample.size, quantile.range=quantile.range, search.radius=search.radius)
    }
    bootstrap.resampling.samples[[r]] <- hadrm3p.idx
  }
  
  if (bootstrap.resampling == 1) return(unlist(bootstrap.resampling.samples))
  return(bootstrap.resampling.samples) 
}



# resampling of ensemble with replacement:
#' @title A function to derive a percentile transfer function 
#' @export
resample.ensemble_replacement <- function(mod.grid.cell, transfer.function, 
                                           sample.size = 1500, quantile.range = c(0, 1)) {
  
  # i. sample from uniform distribution:
  unif.sample = runif(n=sample.size) * diff(quantile.range) + quantile.range[1]
  
  # ii. get quantiles in resampled ensemble:
  mod.sorted.idx = ceiling(transfer.function(x=unif.sample) * length(mod.grid.cell))
  
  # iii. Translate into ensemble indices in the original ensemble:
  mod.grid.cell.sorted.idx = sort.int(mod.grid.cell, decreasing=F, index.return=T)$ix
  hadrm3p.idx = mod.grid.cell.sorted.idx[mod.sorted.idx]
  # hadrm3p.idx = sapply(X=sort(mod.grid.cell, decreasing=F)[mod.sorted.idx], FUN=function(idx) which(idx == mod.grid.cell)[1])
  
  return(hadrm3p.idx)
}

# 
#' @title A function to derive a percentile transfer function 
#' @export
resample.ensemble_no_replacement <- function(mod.grid.cell, transfer.function, 
                                             sample.size = 1500, quantile.range = c(0, 1),
                                             search.radius) {
    
    # i. sample from uniform distribution:
    unif.sample = runif(n=sample.size) * diff(quantile.range) + quantile.range[1]
    # hist(unif.sample)
    
    # ii. get quantiles in resampled ensemble:
    mod.sorted.idx = ceiling(transfer.function(x=unif.sample) * length(mod.grid.cell))
    
    # iii. Prepare for loop
    mod.sorted.idx.new = NA
    mod.grid.cell.sorted.idx = sort.int(mod.grid.cell, decreasing=F, index.return=T)$ix # indices for sorting grid sells
    mod.grid.cell.sorted = mod.grid.cell[mod.grid.cell.sorted.idx] # sort grid cells
    mod.grid.cell.sorted.selected = mod.grid.cell.sorted  # vector to label ensemble members that have been selected

    # iv. run through each resampled ensemble index, grab closest member that is within search space:
    for (i in 1:sample.size) {
      print(i)
  
      # Ensemble members that fulfill criterion:
      idx = which(mod.grid.cell.sorted > (mod.grid.cell.sorted[mod.sorted.idx[i]] - search.radius) & mod.grid.cell.sorted < (mod.grid.cell.sorted[mod.sorted.idx[i]] + search.radius))
      
      # Distance in temperature units from selected individual:
      dist.vec = mod.grid.cell.sorted.selected[idx] - mod.grid.cell.sorted[mod.sorted.idx[i]]
      
      # Choose minimum distance and put into resampling vector:
      if (all(is.na(dist.vec))) { 
        print("No ensemble member found any more!")
        break }
      
      cur.idx = idx[which.min(abs(dist.vec))]
      # put into resampling vector and label "NA" those guys that have been selected already:
      mod.sorted.idx.new[i] = cur.idx
      mod.grid.cell.sorted.selected[cur.idx] <- NA
    }
    
  # v. Translate into ensemble indices in the original ensemble:  
    hadrm3p.idx = rep(NA, times=sample.size)
    hadrm3p.idx[1:length(mod.sorted.idx.new)] <- mod.grid.cell.sorted.idx[mod.sorted.idx.new]
    # hadrm3p.idx = sapply(X=sort(mod.grid.cell, decreasing=F)[mod.sorted.idx.new], FUN=function(idx) which(idx == mod.grid.cell)[1])
    
  return(hadrm3p.idx)
}

# test_no_replacement = resample.ensemble_no_replacement(mod.grid.cell=mod.data$Tair, transfer.function=cur.fun, sample.size=8000, quantile.range=c(0,1), search.radius=0.5)
# test_replacement = resample.ensemble_replacement(mod.grid.cell=mod.data$Tair, transfer.function=cur.fun, sample.size=10000, quantile.range=c(0,1))






## Function to provide quality control for ensemble resampling:
# --------------------------------------------------------------
# Input: 
# transfer.function: spline function that had been fitted with "get.percentile.transfer.function"
# percentile.ranges: the intervals (of the original data percentiles) in which a test is performed

# Value:
# Fraction of original ensemble members that are in the respective bins.

#' @title A function to derive a percentile transfer function 
#' @export
quality.control.ensemble <- function(transfer.function, percentile.ranges = seq(from=0, to=1, by = 0.1) ) {
  return(round(sapply(X=1:(length(percentile.ranges)-1), FUN=function(idx) transfer.function(percentile.ranges[idx+1]) - transfer.function(percentile.ranges[idx])), 2))
}






