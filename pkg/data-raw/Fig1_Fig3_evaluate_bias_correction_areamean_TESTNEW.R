# ---------------------------------------
# Bias correction & evaluation script
# This script is to evaluate the resampling based bias correction over Europe
# ---------------------------------------

# Sebastian Sippel
# 06.07.2015


library(raster)
library(ncdf4)
library(KernSmooth)
library(sROC)


source("/Users/ssippel/code/tools/convert.to.eurocentric.R")
source("/Users/ssippel/code/tools/frenchcolormap.R")
source("/Users/ssippel/code/tools/SSA_tools.R")
source("/Users/ssippel/code/bias_correction/old_code/resampling_bias_correction.R")
source("/Users/ssippel/code/tools/SSA_tools.R")



# setwd("/Users/ssippel/projects/in-progress/LPJEns/data/HadRM3P_oRiginal/")
hadrm3p.extent = brick("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data/CEurope_cutout/1985/hadam3p_eu_ins0_1985.nc", varname="field16") + 0
hadrm3p.mask = subset(crop(brick("/Users/ssippel/data/grid/0d50_monthly/EOBS/Tair_Europe_components/Tair.1985-2010.monthly.seas.var.trig.reg.2harmonics_EOBS.nc") + 0, y=extent(hadrm3p.extent)), 1)


# Function Definitions
# ------------------------------------------------
# functions to compute indices for subsets and split into list with different years:
compute.HadRM3P.idx <- function(start.year, NR.years = 26) {
  c((1986 - start.year) * 12) : c((1986 - start.year) * 12 + NR.years*12-1)
}

# aux function to read, crop and subset rasterbrick:
read.OBS <- function(file.name, start.year, mask = hadrm3p.mask, idx = NA)  {
  if (any(is.na(idx))) {
    temp.raster = mask(subset(crop(brick(file.name) + 0, y=extent(mask)), subset=compute.HadRM3P.idx(start.year=start.year)), mask=mask)
  } else {
    temp.raster = mask(subset(crop(brick(file.name) + 0, y=extent(mask)), subset=idx), mask=mask)
  }
  projection(temp.raster) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")  # original data is longitude/latitude
  return(temp.raster) 
}

# Define function for trigonometric correction:
trigonometric.correction.HadRM3P <- function(trend.list, seas.list, remainder.list) {
  return( lapply(1:length(remainder.list), FUN=function(year.idx) lapply(remainder.list[[year.idx]], FUN=function(x) x + trend.list[[year.idx]] + seas.list[[year.idx]])) )
}


get.JJA.mean <- function(original.raster, idx.extract = c(6,7,8)) {
  
  # get seasonal means
  test.raster = brick(sapply(0 : (nlayers(original.raster) / 12 - 1), FUN=function(idx) {
    return(calc(subset(original.raster, subset = (idx * 12) + idx.extract), fun=mean)) }))
  return(test.raster)
}


# detrend each grid cell:
get.trend.comp <- function(original.raster.JJA, pad.series = c(30, 30)) {
  # linear padding of the JJA series:
  monthly_0d50.res.pad = calc(x= original.raster.JJA, fun= function(time.series) pad.linear(x=time.series, pad.series=pad.series, linear.pred=c(30, 30)) )
  # 3. extract trend of the remainder:
  monthly_0d50.trend.pad = calc(x=monthly_0d50.res.pad, fun=function(x) SSA.rasterbrick(data=x, borders.wl=list(c(31,Inf)), M = 62, n.comp = 20))
  monthly_0d50.trend = subset(monthly_0d50.trend.pad, subset=c(c(pad.series[1] + 1) : c(nlayers(x=monthly_0d50.trend.pad) - pad.series[2])))
  
  return(monthly_0d50.trend)
}


get.seas.mean.areamean.HadRM3P <- function(HadRM3P.list, nr.samples.per.year = 10, stat=mean, mon.idx = c(7,8,9)) {
  
  list.out = lapply(1:length(HadRM3P.list), FUN=function(year.idx) {
    sapply(X=HadRM3P.list[[year.idx]][1:nr.samples.per.year], FUN=function(x) {
      test = mean(cellStats(x=x, stat)[mon.idx])
      return(test)
    })
  })
  return(unlist(list.out))
}


cut.kernel.density = function(data, probs = c(0.01, 0.99), zero.mean = T) {
  data = na.omit(c(data))
  quantile.est = quantile(x=data, probs=probs)
  data[which(data > quantile.est[1] & data < quantile.est[2])]
  
  if (zero.mean == T) data = data - mean(data)
  
  return(data)
}






# Read bias corrected data:
# ------------------------------------------------
# no single grid cell evaluation (or: later!)
# make example for areamean constraint:

# load resampled rasters from ERA-Interim (areamean):
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data/HadRM3P_bias_corrected/univariate/Tair/HadRM3P_probabilistic_resampling_area/HadRM3P.field16.ERAI.prob_resampling_areamean.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data/HadRM3P_bias_corrected/univariate/Tair/HadRM3P_probabilistic_resampling_area/HadRM3P.field16.ERAI.prob_resampling_areamean_INFO.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/precip/HadRM3P_probabilistic_resampling_area/HadRM3P.field90.ERAI.prob_resampling_areamean.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/SWdown/HadRM3P.field203.ERAI.prob_resampling_areamean.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/LWdown/HadRM3P.field205.ERAI.prob_resampling_areamean.RData")

# load original HadRM3P data:
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_original_monthly/HadRM3P.field16.original.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_original_monthly/HadRM3P.field90.original.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_original_monthly/HadRM3P.field203.original.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_original_monthly/HadRM3P.field205.original.RData")

# load ERA-Interim ISIMIP corrected data:
# ---------------------------
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/Tair/HadRM3P_ISIMIP/HadRM3P.field16.ISIMIP.ERAI.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/precip/HadRM3P_trigonometric_correction/HadRM3P.field90.seascor_6harmonics.ERAI.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/SWdown/HadRM3P.field203.ERAI.ISIMIP.RData")
load("/Users/ssippel/projects/published/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/LWdown/HadRM3P.field205.ERAI.ISIMIP.RData")



# load monthly means of ERA-Interim and of the ensemble:
brick("/Volumes/BGI/people/ssippel/data/grid/0d50_monthly/Oxpewwes/components/Oxpewwes_ensmean.nc", varname = "field16")


# load ERA-Interim data that was used for bias-correction:
# ---------------------------------------------------------
ERAI.Tair =   -273.15 + mask(crop(convert.to.eurocentric(brick("/Users/ssippel/data/grid/0d50_monthly/ERAInterim_V2/Tair/ERAI.Tair.720.360.nc")), hadrm3p.extent), hadrm3p.mask)
ERAI.Precip = 1000 * mask(crop(convert.to.eurocentric(brick("/Users/ssippel/data/grid/0d50_monthly/ERAInterim_V2/precip/ERAI.Precip.720.360.nc")), hadrm3p.extent), hadrm3p.mask)
ERAI.SWdown = 10^6 / 24 / 3600 * mask(crop(convert.to.eurocentric(brick("/Users/ssippel/data/grid/0d50_monthly/ERAInterim_V2/ERAI.SWdown.720.360.nc")), hadrm3p.extent), hadrm3p.mask)
ERAI.LWdown = 10^6 / 24 / 3600 * mask(crop(convert.to.eurocentric(brick("/Users/ssippel/data/grid/0d50_monthly/ERAInterim_V2/ERAI.LWdown.720.360.nc")), hadrm3p.extent), hadrm3p.mask)


HadRM3P.field16.ERAI.prob_resampling_areamean




# 1. GET ERA-INTERIM OBSERVATIONS
# ---------------------------------------------------------

# iv. ERAInterim_V2: areamean
Tair_Europe_ERAI.original = read.OBS(file.name="/Users/ssippel/data/grid/0d50_monthly/ERAInterim_V2/Tair_Europe_components/Tair.1979-2014.original_ERAI.nc", idx = c(1:432))
Tair_Europe_ERAI.JJA.area.original = cellStats(get.JJA.mean(original.raster=Tair_Europe_ERAI.original, idx.extract= c(6,7,8)), stat="mean")
# Tair_Europe_ERAI.JJA.area.trend = filterTSeriesSSA(series=pad.linear(x=Tair_Europe_ERAI.JJA.area.original, pad.series=c(30, 30), linear.pred = c(30, 30)), borders.wl=list(c(31,Inf)), M=62, n.comp=20)$dec.series[1,(1:length(Tair_Europe_ERAI.JJA.area.original) + 30)]

# get JJA means for Precip, SWdown, LWdown:
# ---------------------------------------------------------
Precip_Europe_ERAI.original = cellStats(x=get.JJA.mean(original.raster=ERAI.Precip, idx.extract=c(6,7,8)), stat="mean") * 3
SWdown_Europe_ERAI.original = cellStats(x=get.JJA.mean(original.raster=ERAI.SWdown, idx.extract=c(6,7,8)), stat="mean")
LWdown_Europe_ERAI.original = cellStats(x=get.JJA.mean(original.raster=ERAI.LWdown, idx.extract=c(6,7,8)), stat="mean")

#

# 2a. estimate seasonal means for HadRM3P, including trends:
# ------------------------------------------------------
# make a list of 50 randomly selected files and # get seasonal estimates for HadRM3P: THIS IS THE TRAINING DATA FOR THE RESAMPLING STEP
HadRM3P.field16.list = unlist(lapply(X = 1:26, FUN=function(year.idx) HadRM3P.field16.original[[year.idx]][1:500] ), recursive=F)
mon.idx.list = list(c(1,2,3), c(4,5,6), c(7,8,9), c(10,11,12))
HadRM3P.field16.original_seas = lapply(X=mon.idx.list, FUN=function(idx) unlist(lapply(X=HadRM3P.field16.list, FUN=function(x) cellStats(x=calc(x=subset(x, subset= idx), fun=mean), stat="mean") - 273.15)) )

# 2a. estimate areameans:
JJA.idx = c(7,8,9)
HadRM3P.field16.JJA.means = lapply(X = 1:26, FUN=function(year.idx) {
  print(year.idx)
  return( sapply(HadRM3P.field16.original[[year.idx]], FUN=function(x) cellStats(x=calc(x=subset(x, subset= JJA.idx), fun=mean), stat="mean")) )
} )

HadRM3P.field90.JJA.means = lapply(X = 1:26, FUN=function(year.idx) {
  print(year.idx)
  return( sapply(HadRM3P.field90.original[[year.idx]], FUN=function(x) cellStats(x=calc(x=subset(x, subset= JJA.idx), fun=function(x) mean(x)), stat="mean") * 3) )
} )

HadRM3P.field203.JJA.means = lapply(X = 1:26, FUN=function(year.idx) {
  print(year.idx)
  return( sapply(HadRM3P.field203.original[[year.idx]], FUN=function(x) cellStats(x=calc(x=subset(x, subset= JJA.idx), fun=function(x) mean(x)), stat="mean")) )
} )

HadRM3P.field205.JJA.means = lapply(X = 1:26, FUN=function(year.idx) {
  print(year.idx)
  return( sapply(HadRM3P.field205.original[[year.idx]], FUN=function(x) cellStats(x=calc(x=subset(x, subset= JJA.idx), fun=function(x) mean(x)), stat="mean")) )
} )




# Save data for individual testing of resampling procedure:
years = 1985:2010
OBS.Tair = Tair_Europe_ERAI.JJA.area.original[7:32]
OBS.Precip = Precip_Europe_ERAI.original[7:32]
OBS.SWdown = SWdown_Europe_ERAI.original[7:32]
OBS.LWdown = LWdown_Europe_ERAI.original[7:32]

obs.data = data.frame(years, OBS.Tair, OBS.Precip, OBS.SWdown, OBS.LWdown)
names(obs.data) = c("Year", "Tair", "Precip", "SWdown", "LWdown")

# get model data for different variables:
years = c(sapply(X=1985:2010, FUN=function(x) rep(x, 800)))
length(years)
umid=substr(x=unlist(HadRM3P.field16.file.list), start=89, stop=92)
length(umid)

mod.Tair = (unlist(HadRM3P.field16.JJA.means))
mod.Precip = (unlist(HadRM3P.field90.JJA.means)) * 24 * 3600 * 30
mod.LWdown = (unlist(HadRM3P.field205.JJA.means))
mod.SWdown = (unlist(HadRM3P.field203.JJA.means))
length(mod.SWdown)

# qqplot(x=OBS.LWdown, y=mod.LWdown)
mod.data = data.frame(umid[1:20798], years[1:20798], mod.Tair, mod.Precip, mod.LWdown, mod.SWdown)
names(mod.data) = c("ID", "Year", "Tair", "Precip", "SWdown", "LWdown")

# Save obs/model data:

# What are variables that are needed?

# 1. Obs + mod for calibration of transfer function (26 years)
# 2. large ensemble for resampling


# "Didactic" ILLUSTRATION of resampling bias correction:
# -----------------------------------------
# a. QQ-plot
qqplot(x=Tair_Europe_ERAI.JJA.area.original, y=unlist(HadRM3P.field16.JJA.means) - 273.15)
lines(x=c(0,100), c(0,100), col="red")



# Make evaluation plots:
# ----------------------------------
# Fig. 1a: Areamean (demonstration of correction): Demonstration of methodology
# Fig. 1b: Show that probabilistic resampling improves multivariate simulation


# Appendix: 
# 1. improvement of precipitation on different time scales
# 1. Annual sums of precipitation
# 2. Maps of areamean, gridcell, original
 

# FIG 1. plot Gaussian kernel demonstration for areamean:
# ---------------------------------------------------------

setwd("/Users/ssippel/projects/in-progress/bias_correction/code_new/02_methods/")

pdf("Fig1a_illustration_resampling_CEurope_areamean.pdf", width=9, height=5)
par(mar=c(4,4,4,2), mfrow=c(1,2), las = 1)

kernel.sample = sapply(X= sample(x=Tair_Europe_ERAI.JJA.area.original, size=1000, replace = T), FUN=function(x) rnorm(n=1, mean=x, sd=dpik(x=Tair_Europe_ERAI.JJA.area.original)))
kernel.quantile.seq = seq(0, 1, 0.001)
kernel.percentiles = quantile(x=kernel.sample, probs= kernel.quantile.seq)
model.percentile = (which(abs(sort(unlist(HadRM3P.field16.JJA.means) - 273.15, decreasing = F) - median(kernel.percentiles)) == min(abs(sort(unlist(HadRM3P.field16.JJA.means) - 273.15, decreasing = F) - median(kernel.percentiles)))) / length(unlist(HadRM3P.field16.JJA.means)))[1]

plot(ecdf(Tair_Europe_ERAI.JJA.area.original[7:33]), bty='n', 
     xlab = "Temperature, JJA mean [°C]", ylab="Cumulative Density Function", xlim=c(14,25), main="Central Europe, areamean", yaxt="n")
axis(side=2, las=1)
lines(x = kernel.percentiles, y=kernel.quantile.seq, col="blue", lwd=2)
lines(ecdf(unlist(HadRM3P.field16.JJA.means) - 273.15), col="red", lwd = 2)

lines(x = c(0, median(kernel.percentiles)), y = c(0.5, 0.5), lty=2, lwd = 2)
lines(x = c(median(kernel.percentiles), median(kernel.percentiles)), y = c(0.5, model.percentile), lty=2, lwd=2)
lines(x = c(0, median(kernel.percentiles)), y = c(model.percentile, model.percentile), lty=2, lwd=2)
legend("bottomright", c("ERA-Interim", "Gaussian Kernel, Obs", "Gaussian Kernel, Mod"), bty='n', col=c("black", "blue", "red"), lty=c(NA, 1, 1), pch=c(16, NA,NA), lwd= c(2, 2, 2), cex=0.9)

# ii. Plot the non-parametric transfer function (taken from the bias correction that was conducted!):

plot(x = kernel.percentiles, y = HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](kernel.quantile.seq) * 100, xlab="Observations [°C]", ylab="Ensemble Percentile", bty="n",
     ylim=c(0, 100))
lines(x = kernel.percentiles, y = HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](kernel.quantile.seq) * 100, col="red")
axis(side=3, at=kernel.percentiles[seq(1, 1001, 100)], labels=seq(0, 100, 10))
mtext(text="Obs. Kernel Percentiles", side=3, line=2)

lines(x = c(median(kernel.percentiles), median(kernel.percentiles)), y = c(100, HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](0.5) * 100), lty=2, lwd=2)
lines(x = c(0,median(kernel.percentiles)), y = c(HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](0.5) * 100, 100 * HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](0.5)), lty=2, lwd=2)
legend("bottomright", c("Hermite Splines"), bty='n', lty=1, col="red")

par(new=T, mar=c(0,0,0,0), oma=c(0,0,0,0), font=2, mfrow=c(1,1))
plot.new()
legend(x=-0.06,y=1.05, "a", bty='n')
legend(x=0.5,y=1.05, "b", bty='n')
dev.off()







## Alternative figure using new illustration:



setwd("/Users/ssippel/projects/in-progress/bias_correction/code_new/02_methods/")

pdf("Fig1a_illustration_resampling_CEurope_areamean_NEW.pdf", width=9, height=5)
par(mar=c(4,4,4,4), mfrow=c(1,2), las = 1)

kernel.sample = sapply(X= sample(x=Tair_Europe_ERAI.JJA.area.original, size=1000, replace = T), FUN=function(x) rnorm(n=1, mean=x, sd=dpik(x=Tair_Europe_ERAI.JJA.area.original)))
kernel.quantile.seq = seq(0, 1, 0.001)
kernel.percentiles = quantile(x=kernel.sample, probs= kernel.quantile.seq)
model.percentile = (which(abs(sort(unlist(HadRM3P.field16.JJA.means) - 273.15, decreasing = F) - median(kernel.percentiles)) == min(abs(sort(unlist(HadRM3P.field16.JJA.means) - 273.15, decreasing = F) - median(kernel.percentiles)))) / length(unlist(HadRM3P.field16.JJA.means)))[1]

temp = seq(15.5, 20.5, 0.5)
obs.perc = sapply(X=temp, FUN=function(x) 100 * seq(0, 1, 0.001)[which(abs(x - kernel.percentiles) == min(abs(x - kernel.percentiles)))] )
mod.perc = HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](obs.perc/100) * 100

plot(ecdf(Tair_Europe_ERAI.JJA.area.original[7:33]), bty='n', 
     xlab = "Temperature, JJA mean [°C]", ylab="Cumulative Density Function", xlim=c(14,25), main="Central Europe, areamean", yaxt="n")
axis(side=2, las=1)
lines(x = kernel.percentiles, y=kernel.quantile.seq, col="blue", lwd=2)
lines(ecdf(unlist(HadRM3P.field16.JJA.means) - 273.15), col="red", lwd = 2)

lines(x = c(0, median(kernel.percentiles)), y = c(0.5, 0.5), lty=2, lwd = 2)
lines(x = c(median(kernel.percentiles), median(kernel.percentiles)), y = c(0.5, model.percentile), lty=2, lwd=2)
lines(x = c(0, median(kernel.percentiles)), y = c(model.percentile, model.percentile), lty=2, lwd=2)
legend("bottomright", c("ERA-Interim", "Gaussian Kernel, Obs", "Gaussian Kernel, Mod"), bty='n', col=c("black", "blue", "red"), lty=c(NA, 1, 1), pch=c(16, NA,NA), lwd= c(2, 2, 2), cex=0.9)

# ii. Plot the non-parametric transfer function (taken from the bias correction that was conducted!):
plot(x = seq(0, 1, 0.001) * 100, y = HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](kernel.quantile.seq) * 100, xlab="Observations, Percentiles", ylab="Model Ensemble, Percentiles", bty="n",
     ylim=c(0, 100))
lines(x = seq(0, 1, 0.001) * 100, y = HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](kernel.quantile.seq) * 100, col="red")
axis(side=3, at= obs.perc, labels = temp)
axis(side=4, at= mod.perc, labels = temp)

sapply(1:length(obs.perc), FUN=function(idx) lines(x = rep(obs.perc[idx], 2), y = c(mod.perc[idx],100), lty=3, lwd = 1))
sapply(1:length(obs.perc), FUN=function(idx) lines(x = c(obs.perc[idx], 100), y = c(mod.perc[idx],mod.perc[idx]), lty=3, lwd = 1))

# get median example in here: (!!)
lines(x = c(50, 50), y = c(100, HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](0.5) * 100), lty=2, lwd=2)
lines(x = c(50, 100), y = rep(HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](0.5) * 100, 2), lty=2, lwd=2)

# get median example for resampling:
lines(x = c(50, 50), y = c(0, HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](0.5) * 100), lty=2, lwd=2, col="darkred")
lines(x = c(0, 50), y = rep(HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]](0.5) * 100, 2), lty=2, lwd=2, col="darkred")

# mtext for axis side 3+4 labels:
mtext(text="Observations, Tair [°C]", side=3, line=2)
mtext(text="Model Ensemble, Tair [°C]", side=4, line=3, las=3)


legend("topleft", c("Hermite Splines"), bty='n', lty=1, col="red")
par(new=T, mar=c(0,0,0,0), oma=c(0,0,0,0), font=2, mfrow=c(1,1))
plot.new()
legend(x=-0.06,y=1.05, "a", bty='n')
legend(x=0.5,y=1.05, "b", bty='n')

dev.off()



# make QQ-plot of the resampled ensemble:
# --------------------------------------
kernel.percentiles # quantile estimates in the ensemble

transfer.function = get.percentile.transfer.function(obs.grid.cell=Tair_Europe_ERAI.JJA.area.original, mod.grid.cell= unlist(HadRM3P.field16.JJA.means) - 273.15, n.sample=1000)



# resample to get PROBCOR adjustment:
# ------------------------------------
HadRM3P.field16.JJA.ORIG.means = c(unlist(HadRM3P.field16.JJA.means))
idx = resample.ensemble(mod.grid.cell=HadRM3P.field16.JJA.ORIG.means, transfer.function=transfer.function, sample.size= 1000)

# make quick bootstrap scheme to determine resampling scheme:
# -------------------------------------


## -----------------------------------------------
# Test figure for bias correction
## -----------------------------------------------

# Sebastian Sippel
# 28.08.2015

# source bias correction functions:
source("/Users/ssippel/code/bias_correction/ensemble_bias_correction.R")


# setwd to bias correction:
setwd("/Users/ssippel/projects/in-progress/bias_correction/data/HadRM3P_bias_corrected/univariate/Tair/HadRM3P_probabilistic_resampling_area/")
## load ERAI:
load("HadRM3P.field16.ERAI.prob_resampling_areamean_INFO.RData")

### load grid cell stuff:
setwd("/Users/ssippel/projects/in-progress/bias_correction/data/HadRM3P_bias_corrected/univariate/Tair/HadRM3P_probabilistic_resampling_cell/")
# ... files: HadRM3P.field16.BerkeleyEarth.prob_resampling_areamean[[1]]
load("HadRM3P.field16.ERAI.prob_resampling_gridcell_INFO.RData")


quality.control.ensemble <- function(transfer.function, percentile.ranges = seq(from=0, to=1, by = 0.1) ) {
  return(round(sapply(X=1:(length(percentile.ranges)-1), FUN=function(idx) transfer.function(percentile.ranges[idx+1]) - transfer.function(percentile.ranges[idx])), 3))
}


# try bar plot with ggplot2:
percentile.ranges = c(0, 0.02, seq(0.1, 0.9, 0.1), 0.98, 1)
# percentile.ranges = seq(0, 1, 0.05)
percentiles.areamean = quality.control.ensemble(transfer.function=HadRM3P.field16.ERAI.prob_resampling_areamean_INFO[[1]], percentile.ranges)

grid.cell = 156 # Jena grid cell
percentiles.gridcell = quality.control.ensemble(transfer.function=HadRM3P.field16.ERAI.prob_resampling_gridcell_INFO[[1]][[grid.cell]], percentile.ranges)

text = c("< 2", "2-10", "10-20",       
         "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-98", "> 98")

sample.size = 800 * 26



# make barplot in base R:
setwd("/Users/ssippel/projects/in-progress/bias_correction/code_new/02_methods/")

pdf(file = "Fig1c_nr-ensembles_TESTNEW.pdf", width = 9, height = 4)
par(mar=c(5,5,1,5),mfrow=c(1,2), cex=0.7)

plot(c(1,1), type='n', xlim = c(14, 24), ylim=c(14,24),
     xlab = "Observations, JJA Tair [°C]", ylab = "Model Ensemble, JJA Tair [°C]", bty='n')
lines(x=c(1, 100), y=c(1,100), col="darkblue", lwd = 2)
points(x=quantile(kernel.percentiles, seq(0.01, 0.99, 0.01)), y = quantile(HadRM3P.field16.JJA.ORIG.means[idx] - 273.15, seq(0.01, 0.99, 0.01)))
points(x=quantile(kernel.percentiles, seq(0.01, 0.99, 0.01)), y = quantile(HadRM3P.field16.ORIG_areamean_JJA, seq(0.01, 0.99, 0.01)), col="red")
legend("bottomright", c("Original ensemble", "Resampled ensemble"), pch = 1, col=c("red", "black"), bty='n')

barplot(percentiles.areamean, width=c(0.2,0.6, rep(0.8,8), 0.6, 0.2), names = text, 
        xlab = "Percentiles in observations", ylab = "Fraction of original ensemble members in each percentile bin")
axis(side=4, at=c(0, 0.05, 0.1, 0.15, 0.2), labels=c(0, 0.05, 0.1, 0.15, 0.2) * sample.size)
mtext(text="# of samples in each percentile bin", side=4, line=2.5, cex=0.75)
legend("topleft", "Area mean constraint", bty='n')

par(new=T, mar=c(0,0,0,0), oma=c(0,0,0,0), font=2, mfrow=c(1,1))
plot.new()
legend(x=-0.06,y=1.05, "c", bty='n')
legend(x=0.5,y=1.05, "d", bty='n')
dev.off()




# -----------------------------------------------------------
# Make evaluation for areamean simulation: (part of Fig. 1, bottom)
# -----------------------------------------------------------


# get seasonal means for field 16:
HadRM3P.field16.ORIG_areamean_JJA = -273.15 + get.seas.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field16.original, nr.samples.per.year=100, stat = mean, mon.idx = c(7,8,9))
HadRM3P.field16.PROBCOR_areamean_JJA = -273.15 + get.seas.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field16.ERAI.prob_resampling_areamean, nr.samples.per.year=100, stat = mean, mon.idx = c(7,8,9))
HadRM3P.field16.ISIMIP_areamean_JJA = get.seas.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field16.ISIMIP.ERAI, nr.samples.per.year=100, stat = mean, mon.idx = c(7,8,9))
ERAI.Tair_JJA_areamean = cellStats(x=get.JJA.mean(original.raster= ERAI.Tair, idx.extract=c(6,7,8)), stat="mean")


# get seasonal sums for field90:
HadRM3P.field90.ORIG_areamean_JJA = 24 * 3600 * 92 * get.seas.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field90.original, nr.samples.per.year=100, stat=mean, mon.idx = c(7,8,9))
HadRM3P.field90.PROBCOR_areamean_JJA = 24 * 3600 * 92 * get.seas.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field90.ERAI.prob_resampling_areamean, nr.samples.per.year=100, stat=mean, mon.idx = c(7,8,9))
HadRM3P.field90.ISIMIP_areamean_JJA = 3 * get.seas.mean.areamean.HadRM3P(HadRM3P.list= HadRM3P.field90.seascor_6harmonics.ERAI, nr.samples.per.year=100, stat=mean, mon.idx= c(7,8,9))
ERAI.Precip_JJA_areamean = 3 * cellStats(x= get.JJA.mean(original.raster=ERAI.Precip), stat="mean")
  

# get seasonal distribution for SWdown:
HadRM3P.field203.ORIG_areamean_JJA = get.seas.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field203.original, nr.samples.per.year=100, stat=mean, mon.idx = c(7,8,9))
HadRM3P.field203.PROBCOR_areamean_JJA = get.seas.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field203.ERAI.prob_resampling_areamean, nr.samples.per.year=100, stat=mean, mon.idx = c(7,8,9))
HadRM3P.field203.ISIMIP_areamean_JJA =  get.seas.mean.areamean.HadRM3P(HadRM3P.list= HadRM3P.field203.ERAI.ISIMIP, nr.samples.per.year=100, stat=mean, mon.idx= c(7,8,9))
ERAI.SWdown_JJA_areamean = cellStats(x= get.JJA.mean(original.raster= ERAI.SWdown), stat="mean")

# LWdown simulations:
HadRM3P.field205.ORIG_areamean_JJA = get.seas.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field205.original, nr.samples.per.year=100, stat=mean, mon.idx = c(7,8,9))
HadRM3P.field205.PROBCOR_areamean_JJA = get.seas.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field205.ERAI.prob_resampling_areamean, nr.samples.per.year=100, stat=mean, mon.idx = c(7,8,9))
HadRM3P.field205.ISIMIP_areamean_JJA =  get.seas.mean.areamean.HadRM3P(HadRM3P.list= HadRM3P.field205.ERAI.ISIMIP, nr.samples.per.year=100, stat=mean, mon.idx= c(7,8,9))
ERAI.LWdown_JJA_areamean = cellStats(x= get.JJA.mean(original.raster= ERAI.LWdown), stat="mean")




# make plot for Figure 1b.
# ---------------------------------------------
source("/Users/ssippel/code/tools/vioplot2_equalarea.R")


# Demonstration of improvement of summer precipitation:
# ------------------------------------------------------
setwd("/Users/ssippel/projects/in-progress/bias_correction/code_new/03_1_HadRM3P_evaluation/")

pdf(file="Fig3_evaluation_areamean.pdf", width=9, height=6)
par(mar=c(4,4,1,1), mfrow=c(2,2), las = 1)

# Seasonal temperature:
plot(c(1,1), type='n', xlim=c(0.5, 4.5), ylim=c(13, 26), 
     ylab="Mean temperature, JJA [°C]", xlab="", bty="n", xaxt="n")
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field16.ORIG_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 1, col="darkred")
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field16.PROBCOR_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 2, col="lightblue")
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field16.ISIMIP_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 3, col="orange")
vioplot2_equalarea(ERAI.Tair_JJA_areamean, add = T, at = 4, col="darkgray")
axis(side=1, at = c(1,2,3,4), labels=c("ORIG", "PROBCOR", "ISIMIP", "OBS"), lty=0)

# Seasonal precipitation:
plot(c(1,1), type='n', xlim=c(0.5, 4.5), ylim=c(0, 500), 
     ylab="Precipitation sum, JJA [°C]", xlab="", bty="n", xaxt="n")
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field90.ORIG_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 1, col="darkred", wex = 50)
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field90.PROBCOR_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 2, col="lightblue", wex = 50)
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field90.ISIMIP_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 3, col="orange", wex = 50)
vioplot2_equalarea(ERAI.Precip_JJA_areamean, add = T, at = 4, col="darkgray", wex = 50)
axis(side=1, at = c(1,2,3,4), labels=c("ORIG", "PROBCOR", "ISIMIP", "OBS"), lty=0)


# SW incoming radiation:
plot(c(1,1), type='n', xlim=c(0.5, 4.5), ylim=c(150, 350), 
     ylab="Mean incoming SW radiation [W m-2]", xlab="", bty="n", xaxt="n")
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field203.ORIG_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 1, col="darkred", wex = 10)
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field203.PROBCOR_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 2, col="lightblue", wex = 10)
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field203.ISIMIP_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 3, col="orange", wex = 10)
vioplot2_equalarea(ERAI.SWdown_JJA_areamean, add = T, at = 4, col="darkgray", wex = 10)
axis(side=1, at = c(1,2,3,4), labels=c("ORIG", "PROBCOR", "ISIMIP", "OBS"), lty=0)

# LW incoming radiation:
plot(c(1,1), type='n', xlim=c(0.5, 4.5), ylim=c(320, 380), 
     ylab="Mean incoming LW radiation [W m-2]", xlab="", bty="n", xaxt="n")
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field205.ORIG_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 1, col="darkred", wex = 5)
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field205.PROBCOR_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 2, col="lightblue", wex = 5)
vioplot2_equalarea(cut.kernel.density(data=HadRM3P.field205.ISIMIP_areamean_JJA, probs=c(0.01, 0.99), zero.mean = F), add = T, at = 3, col="orange", wex = 5)
vioplot2_equalarea(ERAI.LWdown_JJA_areamean, add = T, at = 4, col="darkgray", wex = 5)
axis(side=1, at = c(1,2,3,4), labels=c("ORIG", "PROBCOR", "ISIMIP", "OBS"), lty=0)

par(new=T, mar=c(0,0,0,0), oma=c(0,0,0,0), font=2, mfrow=c(1,1))
plot.new()
legend(x=-0.06,y=1.05, "a", bty='n')
legend(x=0.48,y=1.05, "b", bty='n')

legend(x=-0.06, y= 0.54, "c", bty='n')
legend(x=0.48, y= 0.54, "d", bty='n')
dev.off()







# ---------------------------------------------------------------------
# Appendix figures:
# ---------------------------------------------------------------------

# 1. Year-round monthly evaluation for all variables:
get.mon.mean.areamean.HadRM3P <- function(HadRM3P.list, nr.samples.per.year = 50, mon.idx = c(2:12, 1)) {
  
  list.out = lapply(1:length(HadRM3P.list), FUN=function(year.idx) {
    sapply(X=HadRM3P.list[[year.idx]][1:nr.samples.per.year], FUN=function(x) {
      test = cellStats(x=x, stat="mean")
      return(test)
    })
  })
  
  # swap components: years vs. months:
  mon.list.out = lapply(X=mon.idx, FUN=function(mon) {
    c(sapply(list.out, FUN=function(x) x[mon,]))
  })
  
  return(mon.list.out)
}



# OBS = ERAI.Tair
get.mon.mean.OBS = function(OBS) {
  lapply(X=1:12, FUN=function(mon) cellStats(x=subset(x=OBS, subset=seq(from=0, to=nlayers(OBS)-1, by=12) + mon), stat="mean"))
}
  

# get monthly values of CellStats:
# ----------------------------------
HadRM3P.field16.mon.ORIG = get.mon.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field16.original, nr.samples.per.year=100)
HadRM3P.field16.mon.PROBCOR = get.mon.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field16.ERAI.prob_resampling_areamean, nr.samples.per.year=100)
HadRM3P.field16.mon.ISIMIP = get.mon.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field16.ISIMIP.ERAI, nr.samples.per.year=100)
ERAI.Tair.mon = get.mon.mean.OBS(OBS=ERAI.Tair)


HadRM3P.field90.mon.ORIG = get.mon.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field90.original, nr.samples.per.year=100)
HadRM3P.field90.mon.PROBCOR = get.mon.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field90.ERAI.prob_resampling_areamean, nr.samples.per.year=100)
HadRM3P.field90.mon.ISIMIP = get.mon.mean.areamean.HadRM3P(HadRM3P.list=HadRM3P.field90.seascor_6harmonics.ERAI, nr.samples.per.year=100)
ERAI.Precip.mon = get.mon.mean.OBS(OBS=ERAI.Precip)




# make evaluation of monthly temperature as area mean:
# -----------------------------------------------------





# Tair:
# ---------------
pdf(file="SI1a_Tair_monthly_eval_areamean.pdf", width = 9, height = 4.5)
par(mfrow=c(1,1), las = 1, mar=c(3,5,0,0))
plot(c(1,1), type='n', xlim=c(0.5, 12.5), ylim=c(-5, 28), 
     ylab="Mean monthly temperature [°C]", xlab="", bty="n", xaxt="n")

sapply(X=1:12, FUN=function(mon.idx) vioplot2_equalmode(cut.kernel.density(HadRM3P.field16.mon.ORIG[[mon.idx]] -273.15, zero.mean=F), wex = 0.3, add = T, at = mon.idx - 0.3, side = "right", col="darkred"))
sapply(X=1:12, FUN=function(mon.idx) vioplot2_equalmode(cut.kernel.density(HadRM3P.field16.mon.PROBCOR[[mon.idx]] -273.15, zero.mean = F), wex = 0.3, add = T, at = mon.idx - 0.1, side = "right", col="lightblue"))
sapply(X=1:12, FUN=function(mon.idx) vioplot2_equalmode(cut.kernel.density(HadRM3P.field16.mon.ISIMIP[[mon.idx]], zero.mean = F), wex = 0.3, add = T, at = mon.idx + 0.1, side = "right", col="orange"))
sapply(X=1:12, FUN=function(mon.idx) vioplot2_equalmode(cut.kernel.density(ERAI.Tair.mon[[mon.idx]], zero.mean = F), wex = 0.3, add = T, at = mon.idx + 0.3, side = "right", col="darkgray"))

# include vertical lines to distinguish months:
sapply(X=1:11, FUN=function(x) lines(x=rep(x + 0.5, 2), y=c(-10, 30), lty=2, col="darkgrey"))
axis(side=1, at = c(1:12), labels=month.abb, lty=0)
legend("topleft", c("ORIG", "PROBCOR", "ISIMIP", "ERA-Interim"), fill = c("darkred", "lightblue", "orange", "darkgray"), bty="n")
dev.off()



# make for precipitation:
pdf(file="SI1b_Precip_monthly_eval_areamean.pdf", width = 9, height = 4.5)
par(mfrow=c(1,1), las = 1, mar=c(3,5,0,0))
plot(c(1,1), type='n', xlim=c(0.5, 12.5), ylim=c(0, 300), 
     ylab="Monthly precipitation [mm]", xlab="", bty="n", xaxt="n")

sapply(X=1:12, FUN=function(mon.idx) vioplot2_equalmode(cut.kernel.density(HadRM3P.field90.mon.ORIG[[mon.idx]] * 24 * 3600 * 30, zero.mean=F), wex = 0.3, add = T, at = mon.idx - 0.3, side = "right", col="darkred"))
sapply(X=1:12, FUN=function(mon.idx) vioplot2_equalmode(cut.kernel.density(HadRM3P.field90.mon.PROBCOR[[mon.idx]] * 24 * 3600 * 30, zero.mean = F), wex = 0.3, add = T, at = mon.idx - 0.1, side = "right", col="lightblue"))
sapply(X=1:12, FUN=function(mon.idx) vioplot2_equalmode(cut.kernel.density(HadRM3P.field90.mon.ISIMIP[[mon.idx]], zero.mean = F), wex = 0.3, add = T, at = mon.idx + 0.1, side = "right", col="orange"))
sapply(X=1:12, FUN=function(mon.idx) vioplot2_equalmode(cut.kernel.density(ERAI.Precip.mon[[mon.idx]], zero.mean = F), wex = 0.3, add = T, at = mon.idx + 0.3, side = "right", col="darkgray"))

# include vertical lines to distinguish months:
sapply(X=1:11, FUN=function(x) lines(x=rep(x + 0.5, 2), y=c(-10, 3000), lty=2, col="darkgrey"))
axis(side=1, at = c(1:12), labels=month.abb, lty=0)
legend("topleft", c("ORIG", "PROBCOR", "ISIMIP", "ERA-Interim"), fill = c("darkred", "lightblue", "orange", "darkgray"), bty="n")
dev.off()










# ---------------------------------------------------------
# Make maps of bias correction effect...(for seasonal means):
# ---------------------------------------------------------


# get seas. means and IQR:
# ---------------------------





HadRM3P.list = HadRM3P.field16.original
nr.samples.per.year = 10
stat = mean


get.seas.stat.HadRM3P <- function(HadRM3P.list, nr.samples.per.year = 10, stat=mean) {
  
  list.seas.idx = list(c(1,2,3), c(4,5,6), c(7,8,9), c(10:12))
  
  # get selection of files: 
  HadRM3P.list = unlist(lapply(X=HadRM3P.list, FUN=function(x) x[1:10]), recursive = F)
  
  
  lapply(X=list.seas.idx, FUN=function(idx) )
  
  list.out = lapply(1:length(HadRM3P.list), FUN=function(year.idx) {
    sapply(X=HadRM3P.list[[year.idx]][1:nr.samples.per.year], FUN=function(x) {
      test = mean(cellStats(x=x, stat)[mon.idx])
      return(test)
    })
  })
  return(unlist(list.out))
}








