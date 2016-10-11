# ---------------------------------------
# Script to evaluate bias correction strategy
# ---------------------------------------

# Sebastian Sippel
# 17.03.2016

library(ncdf4)
library(raster)
library(zoo)


# source code for calculation of wbtemp-etc.:
source("/Users/ssippel/code/wetbulbtemp/wetbulbvap.R")
source("/Users/ssippel/projects/in-progress/attribution_CEU2015/scripts/functions.R")


# ---------------------------------------
# load station data:
# HadRM3P - Oxpewwes:
load("/Users/ssippel/projects/in-progress/attribution_CEU2015/data/HadRM3P/03_fullDF/Oxpewwes_BC.RData")

# EOBS - station observations:
load("/Users/ssippel/projects/in-progress/attribution_CEU2015/data/stations/EOBS_stations_GEV.RData")

# load year data
Oxpewwes_list = list.files(path="/Volumes/BGI/projects/LPJEns/HadRM3P4LPJml/ncdf/", pattern=".nc", recursive=T)
Oxpewwes_ORIG = data.frame(as.numeric(substring(text = Oxpewwes_list, first=1, last=4)), substring(text = Oxpewwes_list, first=17, last=20))
names(Oxpewwes_ORIG) <- c("Year", "umid")
str(Oxpewwes_ORIG)


# Match into simple data.frame structure:
# --------------------------------------------------------

# 1. AT-WIE
obs.data = data.frame(AT_Wien_JJA$Year, 
                      AT_Wien_JJA$TG3_max, AT_Wien_JJA$TG7_max, AT_Wien_JJA$TG15_max, AT_Wien_JJA$TG21_max, AT_Wien_JJA$TG29max, 
                      AT_Wien_JJA$TX3_max, 
                      AT_Wien_JJA$WBGT3_max, AT_Wien_JJA$wbtemp3_max, AT_Wien_JJA$TG3, AT_Wien_JJA$RR3)
names(obs.data) <- c("Year", "Tair_3d_seasmax", "Tair_7d_seasmax", "Tair_15d_seasmax", "Tair_21d_seasmax", "Tair_29d_seasmax",
                     "TX_3d_seasmax", "WBGT_3d_seasmax", "WBT_3d_seasmax", "Tair_JJA", "Precip_JJA")

AT_WIE_Oxpewwes_years = Oxpewwes_ORIG$Year[match(x=AT_WIE_Oxpewwes_BC21d$ORIG$umid, table=Oxpewwes_ORIG$umid)]
mod.data = data.frame(AT_WIE_Oxpewwes_BC21d$ORIG$umid, AT_WIE_Oxpewwes_years,
                      AT_WIE_Oxpewwes_BC21d$ORIG$tair_3dmax, AT_WIE_Oxpewwes_BC21d$ORIG$tair_7dmax, AT_WIE_Oxpewwes_BC21d$ORIG$tair_15dmax, AT_WIE_Oxpewwes_BC21d$ORIG$tair_21dmax, AT_WIE_Oxpewwes_BC21d$ORIG$tair_29dmax,
                      AT_WIE_Oxpewwes_BC21d$ORIG$tx_3dmax,
                      AT_WIE_Oxpewwes_BC21d$ORIG$WBGT_3dmax, AT_WIE_Oxpewwes_BC21d$ORIG$wbtemp_3dmax, AT_WIE_Oxpewwes_BC21d$ORIG$tair_JJA, AT_WIE_Oxpewwes_BC21d$ORIG$rr_JJA)
names(mod.data) <- c("umid", "Year", "Tair_3d_seasmax", "Tair_7d_seasmax", "Tair_15d_seasmax", "Tair_21d_seasmax", "Tair_29d_seasmax",
                     "TX_3d_seasmax", "WBGT_3d_seasmax", "WBT_3d_seasmax", "Tair_JJA", "Precip_JJA")
mod.data = mod.data[-unique(which(is.na(mod.data), arr.ind=T)[,1]),]

AT_WIE = list(obs.data, mod.data)
names(AT_WIE) = c("obs.data", "mod.data")


# 2. DE-JEN
obs.data = data.frame(DE_Jena_JJA$Year, 
                      DE_Jena_JJA$TG3_max, DE_Jena_JJA$TG7_max, DE_Jena_JJA$TG15_max, DE_Jena_JJA$TG21_max, DE_Jena_JJA$TG29max, 
                      DE_Jena_JJA$TX3_max, 
                      DE_Jena_JJA$WBGT3_max, DE_Jena_JJA$wbtemp3_max, DE_Jena_JJA$TG3, DE_Jena_JJA$RR3)
names(obs.data) <- c("Year", "Tair_3d_seasmax", "Tair_7d_seasmax", "Tair_15d_seasmax", "Tair_21d_seasmax", "Tair_29d_seasmax",
                     "TX_3d_seasmax", "WBGT_3d_seasmax", "WBT_3d_seasmax", "Tair_JJA", "Precip_JJA")

DE_JEN_Oxpewwes_years = Oxpewwes_ORIG$Year[match(x=DE_JEN_Oxpewwes_BC21d$ORIG$umid, table=Oxpewwes_ORIG$umid)]
mod.data = data.frame(DE_JEN_Oxpewwes_BC21d$ORIG$umid, DE_JEN_Oxpewwes_years,
                      DE_JEN_Oxpewwes_BC21d$ORIG$tair_3dmax, DE_JEN_Oxpewwes_BC21d$ORIG$tair_7dmax, DE_JEN_Oxpewwes_BC21d$ORIG$tair_15dmax, DE_JEN_Oxpewwes_BC21d$ORIG$tair_21dmax, DE_JEN_Oxpewwes_BC21d$ORIG$tair_29dmax,
                      DE_JEN_Oxpewwes_BC21d$ORIG$tx_3dmax,
                      DE_JEN_Oxpewwes_BC21d$ORIG$WBGT_3dmax, DE_JEN_Oxpewwes_BC21d$ORIG$wbtemp_3dmax, DE_JEN_Oxpewwes_BC21d$ORIG$tair_JJA, DE_JEN_Oxpewwes_BC21d$ORIG$rr_JJA)
names(mod.data) <- c("umid", "Year", "Tair_3d_seasmax", "Tair_7d_seasmax", "Tair_15d_seasmax", "Tair_21d_seasmax", "Tair_29d_seasmax",
                     "TX_3d_seasmax", "WBGT_3d_seasmax", "WBT_3d_seasmax", "Tair_JJA", "Precip_JJA")
mod.data = mod.data[-unique(which(is.na(mod.data), arr.ind=T)[,1]),]

DE_JEN = list(obs.data, mod.data)
names(DE_JEN) = c("obs.data", "mod.data")

# 3. NL-DEB
obs.data = data.frame(NL_DeBilt_JJA$Year, 
                      NL_DeBilt_JJA$TG3_max, NL_DeBilt_JJA$TG7_max, NL_DeBilt_JJA$TG15_max, NL_DeBilt_JJA$TG21_max, NL_DeBilt_JJA$TG29max, 
                      NL_DeBilt_JJA$TX3_max, 
                      NL_DeBilt_JJA$WBGT3_max, NL_DeBilt_JJA$wbtemp3_max, NL_DeBilt_JJA$TG3, NL_DeBilt_JJA$RR3)
names(obs.data) <- c("Year", "Tair_3d_seasmax", "Tair_7d_seasmax", "Tair_15d_seasmax", "Tair_21d_seasmax", "Tair_29d_seasmax",
                     "TX_3d_seasmax", "WBGT_3d_seasmax", "WBT_3d_seasmax", "Tair_JJA", "Precip_JJA")

NL_DEB_Oxpewwes_years = Oxpewwes_ORIG$Year[match(x=NL_DEB_Oxpewwes_BC21d$ORIG$umid, table=Oxpewwes_ORIG$umid)]
mod.data = data.frame(NL_DEB_Oxpewwes_BC21d$ORIG$umid, NL_DEB_Oxpewwes_years,
                      NL_DEB_Oxpewwes_BC21d$ORIG$tair_3dmax, NL_DEB_Oxpewwes_BC21d$ORIG$tair_7dmax, NL_DEB_Oxpewwes_BC21d$ORIG$tair_15dmax, NL_DEB_Oxpewwes_BC21d$ORIG$tair_21dmax, NL_DEB_Oxpewwes_BC21d$ORIG$tair_29dmax,
                      NL_DEB_Oxpewwes_BC21d$ORIG$tx_3dmax,
                      NL_DEB_Oxpewwes_BC21d$ORIG$WBGT_3dmax, NL_DEB_Oxpewwes_BC21d$ORIG$wbtemp_3dmax, NL_DEB_Oxpewwes_BC21d$ORIG$tair_JJA, NL_DEB_Oxpewwes_BC21d$ORIG$rr_JJA)
names(mod.data) <- c("umid", "Year", "Tair_3d_seasmax", "Tair_7d_seasmax", "Tair_15d_seasmax", "Tair_21d_seasmax", "Tair_29d_seasmax",
                     "TX_3d_seasmax", "WBGT_3d_seasmax", "WBT_3d_seasmax", "Tair_JJA", "Precip_JJA")
mod.data = mod.data[-unique(which(is.na(mod.data), arr.ind=T)[,1]),]

NL_DEB = list(obs.data, mod.data)
names(NL_DEB) = c("obs.data", "mod.data")



# fit together to a list structure:
ensbiascoR.example2 = list(AT_WIE, DE_JEN, NL_DEB)
names(ensbiascoR.example2) <- c("AT_WIE", "DE_JEN", "NL_DEB")

setwd("/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/pkg/data/")
save(list=c("ensbiascoR.example2"), file="ensbiascoR.example2.RData")


