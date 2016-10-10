# ---------------------------------------------
# Prepare data for "land coupling constraint" example:
# ---------------------------------------------

# Sebastian Sippel
# 08.10.2016


library(raster)
library(ncdf4)
library("SpatialEpi")
library("rgdal")
library(RColorBrewer)
library(plotrix)
library(MASS)

# Read required functions:
source("/Volumes/BGI/people/ssippel/code/tools/vioplot2_equalarea.R")
source("/Volumes/BGI/people/ssippel/code/tools/frenchcolormap.R")
source("/Volumes/BGI/people/ssippel/code/tools/convert.to.eurocentric.R")
source("/Volumes/BGI/people/ssippel/code/tools/project_raster.R")
source("/Volumes/BGI/people/ssippel/code/tools/filled.contour2.R")
source("/Users/ssippel/projects/in-progress/land_coupling/scripts_v3/functions_plotting.R")
source("/Users/ssippel/projects/in-progress/land_coupling/scripts/gridcorts.R")



# Get SREX masks:   #/net/firebolt/data/LANDCLIM/Static/SREX_masks/2.5deg
# ------------------------
setwd("/Users/ssippel/data/grid/2d50_static/SREX_mask/")
SREX_masks = list.files(path="/Users/ssippel/data/grid/2d50_static/SREX_mask/")[-21]
SREX = sapply(X=SREX_masks, FUN=function(x) raster(x) + 0)
names(SREX) = substring(text=SREX_masks, first=1, last=3)

# Read koeppen Geiger for desert mask:
# -------------------------
koeppen = raster("/Users/ssippel/data/grid/0d50_static/Koeppen/koeppen.nc", varname="koeppen")
koeppen_nc = nc_open("/Users/ssippel/data/grid/0d50_static/Koeppen/koeppen.nc")
koeppen_legend = ncvar_get(nc=koeppen_nc, varid="legend")
nc_close(koeppen_nc)
koeppen_desert_mask = (koeppen == 6 | koeppen == 7)

# get data frame as an overview over CMIP5 models:
file.list = list.files(path="/Users/ssippel/projects/in-progress/land_coupling/data/CMIP5_all/1981/VACa/coinc70/", pattern = ".nc")



# Read Observations / VAC, correlation, Tair_metrics:
# ----------------------------------------------
setwd("/Users/ssippel/projects/in-progress/land_coupling/data_v3/R_processed/obs/")
load("00_1989_2005_OBS_all_orig_SREX.RData")
load("00_1989_2005_OBS_all.RData")
load("02_VAC_1989_2005.RData")
load("02_VAC_1989_2005_SREX.RData")
load("01_Tair_ET_PRECIP_1989_2005_metrics.RData")
load("01_Tair_ET_PRECIP_1989_2005_metrics_SREX.RData")

# Read CMIP5 model simulations +  metadata:
# --------------------------------------------
setwd("/Users/ssippel/projects/in-progress/land_coupling/data_v3/R_processed/CMIP5/")
load("00_CMIP5_monthly_metadata.RData")
load("01_CMIP5_VAC_1989_2005.RData")
load("01_CMIP5_VAC_1989_2005_SREX.RData")
load("01_CMIP5_correlation.RData")
load("01_CMIP5_metrics.RData")
load("00_CMIP5.monthly_SREX.RData")


# I. Different steps to take for package example
# a) fit kernel density estimate over all observations,
# b) get constrained ensemble,
# c) plot observations, original and constrained ensemble.
# -> DO implement:

# save model and observations:
# obs.data.VAC, mod.data.VAC, obs.data.ET, mod.data.ET
# alternative: obs.data$VACc, obs.data$T, obs.data$ET





# ----------------------------------------
# 3. plot correlation with T-mean / T-variability within CMIP5-model world.
# ----------------------------------------
# get percentiles for uncertainty assessment:
VAC70.list = list(VAC70_CRU_allET_SREX, VAC70_ERAI_allET_SREX, VAC70_CRUNCEP_allET_SREX)
names(VAC70.list) = c("CRU", "ERAI", "CRUNCEP")

VAC90.list = list(VAC90_CRU_allET_SREX, VAC90_ERAI_allET_SREX, VAC90_CRUNCEP_allET_SREX)
names(VAC90.list) = c("CRU", "ERAI", "CRUNCEP")

VAC70_allLandFluxEVAL_kdepctl = get.kernel.density.percentiles_LandFluxEVAL_3T(obs.list= VAC70.list, obs.list.idx=c(1:14, 19:22), probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
VAC90_allLandFluxEVAL_kdepctl = get.kernel.density.percentiles_LandFluxEVAL_3T(obs.list= VAC90.list, obs.list.idx=c(1:14, 19:22), probs=c(0.05, 0.25, 0.5, 0.75, 0.95))

# Get monthly CMIP5 resampling index (i.e. which models are within observations...)
CMIP5.VAC70.resampling.idx = get.CMIP5.resampling.idx(CMIP5_VAC = CMIP5_VAC_pctl70_SREX, OBS_kde = VAC70_allLandFluxEVAL_kdepctl, SREX=SREX)


## z-transform through the whole ensemble:
# ----------------------------------------

ET_all_orig_SREX.z = lapply(X=names(SREX), FUN=function(cur.SREX) {
  sapply(X=names(ET_all_orig_SREX), FUN=function(mod.name) {
    cur.ET = c(sapply(X=1:3, FUN=function(mon) {
      cur.tas = ET_all_orig_SREX[[mod.name]][[cur.SREX]][seq(from=mon, to=51, by=3)]
      return((cur.tas - mean(cur.tas)) / sd(cur.tas)) }))
  })
})
names(ET_all_orig_SREX.z) <- names(SREX)

Tair_all_orig_SREX.z = lapply(X=names(SREX), FUN=function(cur.SREX) {
  sapply(X=names(Tair_all_orig_SREX), FUN=function(mod.name) {
    cur.ET = c(sapply(X=1:3, FUN=function(mon) {
      cur.tas = Tair_all_orig_SREX[[mod.name]][[cur.SREX]][seq(from=mon, to=51, by=3)]
      return((cur.tas - mean(cur.tas)) / sd(cur.tas)) }))
  })[,1:3]
})
names(Tair_all_orig_SREX.z) <- names(SREX)


# Plot for CMIP5 ensemble:
CMIP5.year.idx = match(x=1989:2005, table=1870:2100)

CMIP5.tas.z = lapply(X=names(SREX), FUN=function(cur.SREX) {
  CMIP5.tas.z = sapply(X=1:301, FUN=function(mod.idx) {
    cur.ET = c(sapply(X=6:8, FUN=function(mon) {
      cur.tas = CMIP5.tas.monthly.arr[[cur.SREX]][mon,CMIP5.year.idx,mod.idx]
      return((cur.tas - mean(cur.tas)) / sd(cur.tas)) }))
    return(cur.ET) })
  return(CMIP5.tas.z)
})
names(CMIP5.tas.z) <- names(SREX)

CMIP5.ET.z = lapply(X=names(SREX), FUN=function(cur.SREX) {
  CMIP5.ET.z = sapply(X=1:301, FUN=function(mod.idx) {
    cur.ET = c(sapply(X=6:8, FUN=function(mon) {
      cur.tas = CMIP5.hfls.monthly.arr[[cur.SREX]][mon, CMIP5.year.idx, mod.idx]
      return((cur.tas - mean(cur.tas)) / sd(cur.tas)) }))
    return(cur.ET) })
  return(CMIP5.ET.z)
})
names(CMIP5.ET.z) <- names(SREX)





# i. Plot differences in VACc between observations and models
# ------------------------------------------------------------- 
{
  cur.SREX = "CEU"
  cur.VAC = "VACc"

    par(mfrow=c(1,1), mar=c(5, 4, 1, 1))
    plot(c(1,1), xlim = c(0.5, 6.5), ylim = c(0,1), type='n', xaxt="n", bty='n', 
         ylab = "VACc", xlab = "")
    legend("topleft", c("CEU"), bty='n')
    
    axis(side=1, at = 1:6, labels=c("MAM-CMIP5", "MAM-OBS", "JJA-CMIP5", "JJA-OBS", "SON-CMIP5", "SON-OBS"), las=2, cex.axis=0.7)
    
    
    # MAM:    
    vioplot2_equalarea(x= CMIP5_VAC_pctl70_SREX$MAM[[cur.SREX]][[cur.VAC]][CMIP5.unique.mod.idx], at=1, add=T, wex=0.05, col="darkred")
    # All-ET
    vioplot2_equalarea(x= c(unlist(VAC70.list$CRU$MAM[[cur.SREX]][[cur.VAC]])[c(1:14, 19:22)],
                            unlist(VAC70.list$ERAI$MAM[[cur.SREX]][[cur.VAC]])[c(1:14, 19:22)],
                            unlist(VAC70.list$CRUNCEP$MAM[[cur.SREX]][[cur.VAC]])[c(1:14, 19:22)]), 
                       at=1.8, wex = 0.02, add=T, col="lightgrey")
    
    lines(x=rep(2, 2), y= VAC70_allLandFluxEVAL_kdepctl$MAM[[cur.SREX]][[cur.VAC]][c(1,5)], col="darkgray")
    lines(x=rep(2, 2), y= VAC70_allLandFluxEVAL_kdepctl$MAM[[cur.SREX]][[cur.VAC]][c(2,4)], col="darkgray", lwd = 2)
    
    # ERAI-ET diagnostic:
    points(x=rep(2, 5), y= unlist(VAC70.list$ERAI$MAM[[cur.SREX]][[cur.VAC]])[1:5], pch = 16, col="darkgray", cex = 0.7)
    points(x=rep(2, 5), y= unlist(VAC70.list$ERAI$MAM[[cur.SREX]][[cur.VAC]])[6:10], pch = 18, col="darkgray", cex = 0.7)
    points(x=rep(2, 4), y= unlist(VAC70.list$ERAI$MAM[[cur.SREX]][[cur.VAC]])[11:14], pch = 17, col="darkgray", cex=0.7)
    
    # median of LandFluxEval-ORIG:
    points(2.2, y = VAC70.list$ERAI$MAM[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.all.median, pch = 15)
    points(2.2, y = VAC70.list$ERAI$MAM[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.diagnostic.median, pch = 16, cex = 0.7)
    points(2.2, y = VAC70.list$ERAI$MAM[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.lsm.median, pch = 18, cex = 0.7)
    points(2.2, y = VAC70.list$ERAI$MAM[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.reanalyses.median, pch = 17, cex = 0.7)
        
    # JJA:
    vioplot2_equalarea(x= CMIP5_VAC_pctl70_SREX$JJA[[cur.SREX]][[cur.VAC]][CMIP5.unique.mod.idx], at=3, add=T, wex=0.05, col="darkred")    
    # All-ET
    vioplot2_equalarea(x= c(unlist(VAC70.list$CRU$JJA[[cur.SREX]][[cur.VAC]])[c(1:14, 19:22)],
                            unlist(VAC70.list$ERAI$JJA[[cur.SREX]][[cur.VAC]])[c(1:14, 19:22)],
                            unlist(VAC70.list$CRUNCEP$JJA[[cur.SREX]][[cur.VAC]])[c(1:14, 19:22)]), 
                       at=3.8, wex = 0.02, add=T, col="lightgrey")
    
    lines(x=rep(4, 2), y= VAC70_allLandFluxEVAL_kdepctl$JJA[[cur.SREX]][[cur.VAC]][c(1,5)], col="darkgray")
    lines(x=rep(4, 2), y= VAC70_allLandFluxEVAL_kdepctl$JJA[[cur.SREX]][[cur.VAC]][c(2,4)], col="darkgray", lwd = 2)
  
    points(x=rep(4, 5), y= unlist(VAC70.list$ERAI$JJA[[cur.SREX]][[cur.VAC]])[1:5], pch = 16, col="darkgray", cex = 0.7)
    points(x=rep(4, 5), y= unlist(VAC70.list$ERAI$JJA[[cur.SREX]][[cur.VAC]])[6:10], pch = 18, col="darkgray", cex = 0.7)
    points(x=rep(4, 4), y= unlist(VAC70.list$ERAI$JJA[[cur.SREX]][[cur.VAC]])[11:14], pch = 17, col="darkgray", cex=0.7)
    
    # median of LandFluxEval-ORIG:
    points(4.2, y = VAC70.list$ERAI$JJA[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.all.median, pch = 15)
    points(4.2, y = VAC70.list$ERAI$JJA[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.diagnostic.median, pch = 16, cex = 0.7)
    points(4.2, y = VAC70.list$ERAI$JJA[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.lsm.median, pch = 18, cex = 0.7)
    points(4.2, y = VAC70.list$ERAI$JJA[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.reanalyses.median, pch = 17, cex = 0.7)
    
    # SON:    
    vioplot2_equalarea(x= CMIP5_VAC_pctl70_SREX$SON[[cur.SREX]][[cur.VAC]][CMIP5.unique.mod.idx], at=5, add=T, wex=0.05, col="darkred")
    
    vioplot2_equalarea(x= c(unlist(VAC70.list$CRU$SON[[cur.SREX]][[cur.VAC]])[c(1:14, 19:22)],
                            unlist(VAC70.list$ERAI$SON[[cur.SREX]][[cur.VAC]])[c(1:14, 19:22)],
                            unlist(VAC70.list$CRUNCEP$SON[[cur.SREX]][[cur.VAC]])[c(1:14, 19:22)]), 
                       at=5.8, wex = 0.02, add=T, col="lightgrey")
    
    lines(x=rep(6, 2), y= VAC70_allLandFluxEVAL_kdepctl$SON[[cur.SREX]][[cur.VAC]][c(1,5)], col="darkgray")
    lines(x=rep(6, 2), y= VAC70_allLandFluxEVAL_kdepctl$SON[[cur.SREX]][[cur.VAC]][c(2,4)], col="darkgray", lwd = 2)
    
    # ERAI-ET diagnostic:
    points(x=rep(6, 5), y= unlist(VAC70.list$ERAI$SON[[cur.SREX]][[cur.VAC]])[1:5], pch = 16, col="darkgray", cex = 0.7)
    points(x=rep(6, 5), y= unlist(VAC70.list$ERAI$SON[[cur.SREX]][[cur.VAC]])[6:10], pch = 18, col="darkgray", cex = 0.7)
    points(x=rep(6, 4), y= unlist(VAC70.list$ERAI$SON[[cur.SREX]][[cur.VAC]])[11:14], pch = 17, col="darkgray", cex=0.7)
    
    # median of LandFluxEval-ORIG:
    points(6.2, y = VAC70.list$ERAI$SON[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.all.median, pch = 15)
    points(6.2, y = VAC70.list$ERAI$SON[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.diagnostic.median, pch = 16, cex = 0.7)
    points(6.2, y = VAC70.list$ERAI$SON[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.lsm.median, pch = 18, cex = 0.7)
    points(6.2, y = VAC70.list$ERAI$SON[[cur.SREX]][[cur.VAC]]$LandFluxEVAL.reanalyses.median, pch = 17, cex = 0.7)
    
    # draw legend to plot:
    # --------------------------------------------------------------------------
    legend("topright", c("Diagnostic", "LSM", "Reanalyses", 
                         "Median Diagnostic (LFE)", "Median LSM (LFE)", "Median Reanalyses (LFE)",
                         "LandFluxEVAL-Median"), bty='n', pch = c(16, 18, 17, 16, 18, 17, 15), cex = 0.7,
           col=c("lightgrey", "lightgrey", "lightgrey", "black", "black", "black", "black"))
}




# ii. Plot different 2-d structures:
# -----------------------------------
cur.SREX = "CEU"

  # 2-d density in CMIP5 & prepare vectors for 2-D kernel:
  CMIP5.unique.mod.good.idx = CMIP5.unique.mod.idx[which(CMIP5.unique.mod.idx %in% CMIP5.VAC70.resampling.idx$JJA[[cur.SREX]]$VACc)]
  y = c(rep(c(Tair_all_orig_SREX.z[[cur.SREX]][,1]), times=22), rep(c(Tair_all_orig_SREX.z[[cur.SREX]][,2]), times=22), rep(c(Tair_all_orig_SREX.z[[cur.SREX]][,3]), times=22))
  x = rep(c(ET_all_orig_SREX.z[[cur.SREX]][,1:22]), 3)
  
  OBS_kde = kde2d(x=x, y = y, lims=c(-3,3,-3,3), n=30)
  CMIP5_kde = kde2d(y=c(CMIP5.tas.z[[cur.SREX]][,CMIP5.unique.mod.idx]), x = c(CMIP5.ET.z[[cur.SREX]][,CMIP5.unique.mod.idx]), lims=c(-3,3,-3,3), n=30)
  CMIP5_constraint_kde = kde2d(y=c(CMIP5.tas.z[[cur.SREX]][,CMIP5.unique.mod.good.idx]), x = c(CMIP5.ET.z[[cur.SREX]][,CMIP5.unique.mod.good.idx]), lims=c(-3,3,-3,3), n=30)
  
  # get VACc and VACb idx for observations and CMIP5 models:
  CRU.VACc.idx = lapply(X=1:22, FUN=function(idx) which(c(Tair_all_orig_SREX.z[[cur.SREX]][,1]) > quantile(x=c(Tair_all_orig_SREX.z[[cur.SREX]][,1]), probs=0.7, na.rm=T) &
                                                          c(ET_all_orig_SREX.z[[cur.SREX]][,idx]) < quantile(x=c(ET_all_orig_SREX.z[[cur.SREX]][,idx]), probs=0.3, na.rm=T)))
  CRU.VACb.idx = lapply(X=1:22, FUN=function(idx) which(c(Tair_all_orig_SREX.z[[cur.SREX]][,1]) > quantile(x=c(Tair_all_orig_SREX.z[[cur.SREX]][,1]), probs=0.7, na.rm=T) &
                                                          c(ET_all_orig_SREX.z[[cur.SREX]][,idx]) > quantile(x=c(ET_all_orig_SREX.z[[cur.SREX]][,idx]), probs=0.7, na.rm=T)))
  ERAI.VACc.idx = lapply(X=1:22, FUN=function(idx) which(c(Tair_all_orig_SREX.z[[cur.SREX]][,2]) > quantile(x=c(Tair_all_orig_SREX.z[[cur.SREX]][,2]), probs=0.7, na.rm=T) &
                                                           c(ET_all_orig_SREX.z[[cur.SREX]][,idx]) < quantile(x=c(ET_all_orig_SREX.z[[cur.SREX]][,idx]), probs=0.3, na.rm=T)))
  ERAI.VACb.idx = lapply(X=1:22, FUN=function(idx) which(c(Tair_all_orig_SREX.z[[cur.SREX]][,2]) > quantile(x=c(Tair_all_orig_SREX.z[[cur.SREX]][,2]), probs=0.7, na.rm=T) &
                                                           c(ET_all_orig_SREX.z[[cur.SREX]][,idx]) > quantile(x=c(ET_all_orig_SREX.z[[cur.SREX]][,idx]), probs=0.7, na.rm=T)))
  CRUNCEP.VACc.idx = lapply(X=1:22, FUN=function(idx) which(c(Tair_all_orig_SREX.z[[cur.SREX]][,3]) > quantile(x=c(Tair_all_orig_SREX.z[[cur.SREX]][,3]), probs=0.7, na.rm=T) &
                                                              c(ET_all_orig_SREX.z[[cur.SREX]][,idx]) < quantile(x=c(ET_all_orig_SREX.z[[cur.SREX]][,idx]), probs=0.3, na.rm=T)))
  CRUNCEP.VACb.idx = lapply(X=1:22, FUN=function(idx) which(c(Tair_all_orig_SREX.z[[cur.SREX]][,3]) > quantile(x=c(Tair_all_orig_SREX.z[[cur.SREX]][,3]), probs=0.7, na.rm=T) &
                                                              c(ET_all_orig_SREX.z[[cur.SREX]][,idx]) > quantile(x=c(ET_all_orig_SREX.z[[cur.SREX]][,idx]), probs=0.7, na.rm=T)))
  
  # get idx for CMIP5:
  VACc.idx.CMIP5 = lapply(X=1:301, FUN=function(idx) which(c(CMIP5.tas.z[[cur.SREX]][,idx]) > quantile(x=c(CMIP5.tas.z[[cur.SREX]][,idx]), probs=0.7, na.rm=T) &
                                                             c(CMIP5.ET.z[[cur.SREX]][,idx]) < quantile(x=c(CMIP5.ET.z[[cur.SREX]][,idx]), probs=0.3, na.rm=T)))
  VACb.idx.CMIP5 = lapply(X=1:301, FUN=function(idx) which(c(CMIP5.tas.z[[cur.SREX]][,idx]) > quantile(x=c(CMIP5.tas.z[[cur.SREX]][,idx]), probs=0.7, na.rm=T) &
                                                             c(CMIP5.ET.z[[cur.SREX]][,idx]) > quantile(x=c(CMIP5.ET.z[[cur.SREX]][,idx]), probs=0.7, na.rm=T)))
  
  # get correlations and coincidence rates:
  OBS_kde_cor = cor(x=x, y = y)
  OBS_VACc_rate = length(c(unlist(CRU.VACc.idx), unlist(ERAI.VACc.idx), unlist(CRUNCEP.VACc.idx))) / (length(x) / 100 * 30)
  OBS_VACb_rate = length(c(unlist(CRU.VACb.idx), unlist(ERAI.VACb.idx), unlist(CRUNCEP.VACb.idx))) / (length(x) / 100 * 30)
  
  CMIP5_kde_cor = cor(x=c(CMIP5.tas.z[[cur.SREX]][,CMIP5.unique.mod.idx]), y = c(CMIP5.ET.z[[cur.SREX]][,CMIP5.unique.mod.idx]))
  CMIP5_kde_VACc_rate = length(unlist(VACc.idx.CMIP5[CMIP5.unique.mod.idx])) / (length(CMIP5.unique.mod.idx) * 51 / 100 * 30)
  CMIP5_kde_VACb_rate = length(unlist(VACb.idx.CMIP5[CMIP5.unique.mod.idx])) / (length(CMIP5.unique.mod.idx) * 51 / 100 * 30)
  
  CMIP5_constraint_kde_cor = cor(x=c(CMIP5.tas.z[[cur.SREX]][,CMIP5.unique.mod.good.idx]), y = c(CMIP5.ET.z[[cur.SREX]][,CMIP5.unique.mod.good.idx]))
  CMIP5_constraint_kde_VACc_rate = length(unlist(VACc.idx.CMIP5[CMIP5.unique.mod.good.idx])) / (length(CMIP5.unique.mod.good.idx) * 51 / 100 * 30)
  CMIP5_constraint_kde_VACb_rate = length(unlist(VACb.idx.CMIP5[CMIP5.unique.mod.good.idx])) / (length(CMIP5.unique.mod.good.idx) * 51 / 100 * 30)
  
    
  # -----------------------------------------------
  # Make plot for each SREX region:
  # -----------------------------------------------
#  setwd("/Users/ssippel/projects/in-progress/land_coupling/figures_v3/03_constraint_evaluation/02_bivariate_structure_Tanom-ETanom/")
  frenchcolormap.vegetation <- function(n = 1024) {
    return(colorRampPalette(c("ghostwhite", "orange", "red", "brown"))(n))
  }
  col = frenchcolormap.vegetation(n=20)
  cex.points = 0.3
  
  # ii. a) Observations:
  par(mar=c(4,4,0.5, 0.5))
  filled.contour2(OBS_kde, xlim = c(-3, 3), ylim = c(-3, 3), zlim=c(0, 0.3), col=col, nlevels=20, bty='n', 
                  xlab ="ET standardized anomalies (1989-2005)", ylab="Temperature, standardized anomalies (1989-2005)")
  points(x= x, y= y, col = make.transparent.color("darkgray", alpha=250), pch = 16, cex = 0.2)
  
  # CRU:
  lapply(X=1:22, FUN=function(idx) points(y= c(Tair_all_orig_SREX.z[[cur.SREX]][CRU.VACc.idx[[idx]],1]), x= c(ET_all_orig_SREX.z[[cur.SREX]][CRU.VACc.idx[[idx]],idx]), col = "darkred", pch = 16, cex = cex.points))
  lapply(X=1:22, FUN=function(idx) points(y= c(Tair_all_orig_SREX.z[[cur.SREX]][CRU.VACb.idx[[idx]],1]), x= c(ET_all_orig_SREX.z[[cur.SREX]][CRU.VACb.idx[[idx]],idx]), col = "darkorange", pch = 16, cex = cex.points))
  # ERAI:
  lapply(X=1:22, FUN=function(idx) points(y= c(Tair_all_orig_SREX.z[[cur.SREX]][ERAI.VACc.idx[[idx]],2]), x= c(ET_all_orig_SREX.z[[cur.SREX]][ERAI.VACc.idx[[idx]],idx]), col = "darkred", pch = 16, cex = cex.points))
  lapply(X=1:22, FUN=function(idx) points(y= c(Tair_all_orig_SREX.z[[cur.SREX]][ERAI.VACb.idx[[idx]],2]), x= c(ET_all_orig_SREX.z[[cur.SREX]][ERAI.VACb.idx[[idx]],idx]), col = "darkorange", pch = 16, cex = cex.points))
  # CRUNCEP:
  lapply(X=1:22, FUN=function(idx) points(y= c(Tair_all_orig_SREX.z[[cur.SREX]][CRUNCEP.VACc.idx[[idx]],3]), x= c(ET_all_orig_SREX.z[[cur.SREX]][CRUNCEP.VACc.idx[[idx]],idx]), col = "darkred", pch = 16, cex = cex.points))
  lapply(X=1:22, FUN=function(idx) points(y= c(Tair_all_orig_SREX.z[[cur.SREX]][CRUNCEP.VACb.idx[[idx]],3]), x= c(ET_all_orig_SREX.z[[cur.SREX]][CRUNCEP.VACb.idx[[idx]],idx]), col = "darkorange", pch = 16, cex = cex.points))
  legend("topleft", paste(c(cur.SREX), "-OBS", sep=""), bty='n')
  legend("bottomleft", c(paste("r = ", round(OBS_kde_cor, 2), sep=""),
                         paste("VACb-rate = ", round(OBS_VACb_rate*100, 1), "%", sep=""),
                         paste("VACc-rate = ", round(OBS_VACc_rate*100, 1), "%", sep="")), bty='n')
  legend("bottomright", c("VACb-indiv.", "VACc-indiv."), pch=16, col=c("darkorange", "darkred"), bty='n')
  
  
  # ii. b) CMIP5-constrained:
  par(mar=c(4,4,0.5, 0.5))
  filled.contour2(CMIP5_constraint_kde, xlim = c(-3, 3), ylim = c(-3, 3), zlim=c(0, 0.3), col=col, nlevels=20, bty='n', 
                  xlab ="ET standardized anomalies (1989-2005)", ylab="Temperature, standardized anomalies (1989-2005)")
  points(y= c(CMIP5.tas.z[[cur.SREX]][,CMIP5.unique.mod.good.idx]), x= c(CMIP5.ET.z[[cur.SREX]][,CMIP5.unique.mod.good.idx]), col = "darkgray", pch = 16, cex = 0.2)
  lapply(X= CMIP5.unique.mod.good.idx, FUN=function(idx) points(y= c(CMIP5.tas.z[[cur.SREX]][VACc.idx.CMIP5[[idx]],idx]), x= c(CMIP5.ET.z[[cur.SREX]][VACc.idx.CMIP5[[idx]],idx]), col = "darkred", pch = 16, cex = cex.points))
  lapply(X= CMIP5.unique.mod.good.idx, FUN=function(idx) points(y= c(CMIP5.tas.z[[cur.SREX]][VACb.idx.CMIP5[[idx]],idx]), x= c(CMIP5.ET.z[[cur.SREX]][VACb.idx.CMIP5[[idx]],idx]), col = "darkorange", pch = 16, cex = cex.points))  
  legend("topleft", paste(c(cur.SREX), "-CMIP5-constraint", sep=""), bty='n')
  legend("bottomleft", c(paste("r = ", round(CMIP5_constraint_kde_cor, 2), sep=""),
                         paste("VACb-rate = ", round(CMIP5_constraint_kde_VACb_rate*100, 1), "%", sep=""),
                         paste("VACc-rate = ", round(CMIP5_constraint_kde_VACc_rate*100, 1), "%", sep="")), bty='n')
  legend("bottomright", c("VACb-indiv.", "VACc-indiv."), pch=16, col=c("darkorange", "darkred"), bty='n')
  
  
  # ii. c) CMIP5-original:
  par(mar=c(4,4,0.5, 0.5))
  filled.contour2(CMIP5_kde, xlim = c(-3, 3), ylim = c(-3, 3), zlim=c(0, 0.3), col=col, nlevels=20, bty='n', 
                  xlab ="ET standardized anomalies (1989-2005)", ylab="Temperature, standardized anomalies (1989-2005)")
  points(y= c(CMIP5.tas.z[[cur.SREX]][,CMIP5.unique.mod.idx]), x= c(CMIP5.ET.z[[cur.SREX]][,CMIP5.unique.mod.idx]), col = "darkgray", pch = 16, cex = 0.2)
  lapply(X= CMIP5.unique.mod.idx, FUN=function(idx) points(y= c(CMIP5.tas.z[[cur.SREX]][VACc.idx.CMIP5[[idx]],idx]), x= c(CMIP5.ET.z[[cur.SREX]][VACc.idx.CMIP5[[idx]],idx]), col = "darkred", pch = 16, cex = cex.points))
  lapply(X= CMIP5.unique.mod.idx, FUN=function(idx) points(y= c(CMIP5.tas.z[[cur.SREX]][VACb.idx.CMIP5[[idx]],idx]), x= c(CMIP5.ET.z[[cur.SREX]][VACb.idx.CMIP5[[idx]],idx]), col = "darkorange", pch = 16, cex = cex.points))  
  legend("topleft", paste(c(cur.SREX), "-CMIP5-original", sep=""), bty='n')
  legend("bottomleft", c(paste("r = ", round(CMIP5_kde_cor, 2), sep=""),
                         paste("VACb-rate = ", round(CMIP5_kde_VACb_rate*100, 1), "%", sep=""),
                         paste("VACc-rate = ", round(CMIP5_kde_VACc_rate*100, 1), "%", sep="")), bty='n')
  legend("bottomright", c("VACb-indiv.", "VACc-indiv."), pch=16, col=c("darkorange", "darkred"), bty='n')



