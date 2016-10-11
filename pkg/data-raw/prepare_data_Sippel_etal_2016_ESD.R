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
hadrm3p.extent = brick("/Users/ssippel/projects/published/2016/Sippel_etal_ESD2016_bias_correction//data/CEurope_cutout/1985/hadam3p_eu_ins0_1985.nc", varname="field16") + 0
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
load("/Users/ssippel/projects/published/2016/Sippel_etal_ESD2016_bias_correction/data/HadRM3P_bias_corrected/univariate/Tair/HadRM3P_probabilistic_resampling_area/HadRM3P.field16.ERAI.prob_resampling_areamean.RData")
load("/Users/ssippel/projects/published/2016/Sippel_etal_ESD2016_bias_correction/data/HadRM3P_bias_corrected/univariate/Tair/HadRM3P_probabilistic_resampling_area/HadRM3P.field16.ERAI.prob_resampling_areamean_INFO.RData")
load("/Users/ssippel/projects/published/2016/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/precip/HadRM3P_probabilistic_resampling_area/HadRM3P.field90.ERAI.prob_resampling_areamean.RData")
load("/Users/ssippel/projects/published/2016/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/SWdown/HadRM3P.field203.ERAI.prob_resampling_areamean.RData")
load("/Users/ssippel/projects/published/2016/Sippel_etal_ESD2016_bias_correction/data//HadRM3P_bias_corrected/univariate/LWdown/HadRM3P.field205.ERAI.prob_resampling_areamean.RData")

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
ERAI.Tair =   mask(crop(convert.to.eurocentric(brick("/Users/ssippel/data/grid/0d50_monthly/ERAInterim_V2/Tair/ERAI.Tair.720.360.nc")), hadrm3p.extent), hadrm3p.mask) - 273.15
ERAI.Precip = 1000 * mask(crop(convert.to.eurocentric(brick("/Users/ssippel/data/grid/0d50_monthly/ERAInterim_V2/precip/ERAI.Precip.720.360.nc")), hadrm3p.extent), hadrm3p.mask)
ERAI.SWdown = 10^6 / 24 / 3600 * mask(crop(convert.to.eurocentric(brick("/Users/ssippel/data/grid/0d50_monthly/ERAInterim_V2/ERAI.SWdown.720.360.nc")), hadrm3p.extent), hadrm3p.mask)
ERAI.LWdown = 10^6 / 24 / 3600 * mask(crop(convert.to.eurocentric(brick("/Users/ssippel/data/grid/0d50_monthly/ERAInterim_V2/ERAI.LWdown.720.360.nc")), hadrm3p.extent), hadrm3p.mask)




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
years = 1979:2014
OBS.Tair = Tair_Europe_ERAI.JJA.area.original
OBS.Precip = Precip_Europe_ERAI.original
OBS.SWdown = SWdown_Europe_ERAI.original
OBS.LWdown = LWdown_Europe_ERAI.original

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

length(mod.Precip)
length(mod.SWdown)

# qqplot(x=OBS.LWdown, y=mod.LWdown)
mod.data = data.frame(umid[1:20798], years[1:20798], mod.Tair, mod.Precip, mod.LWdown, mod.SWdown)
names(mod.data) = c("ID", "Year", "Tair", "Precip", "SWdown", "LWdown")

mod.data$Tair = mod.data$Tair - 273.15

# Save obs/model data:
load("/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/pkg/data/HadRM3P_monthly.RData")
save(list=c("mod.data", "obs.data"), file="/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/pkg/data/HadRM3P_monthly.RData")
obs.data = obs.data$Tair + 273.15
save(list=c("mod.data", "obs.data"), file="/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/pkg/data/HadRM3P_monthly.RData")

load("/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/pkg/data/HadRM3P_monthly.RData")
ensbiascor.example1 = list(obs.data, mod.data)
names(ensbiascor.example1) = c("obs.data", "mod.data")
save(list=c("ensbiascor.example1"), file="/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/pkg/data/ensbiascor.example1.RData")

load("/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/pkg/data/ensbiascor.example1.RData")
ensbiascoR.example1 = ensbiascor.example1
save(list=c("ensbiascoR.example1"), file="/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/pkg/data/ensbiascoR.example1.RData")




