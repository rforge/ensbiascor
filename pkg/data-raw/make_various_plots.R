# ----------------------------------------------
# Make various plots:
# ----------------------------------------------

# Sebastian Sippel
# 27.09.2016

load("/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/pkg/data/ensbiascor.example1.RData")
# source("/Users/ssippel/code/tools/frenchcolormap.R")


# Run "resampling sequence" according to ESD paper:
# -----------------------------------------------------
setwd("/Users/ssippel/code/R_packages/ensbiascoR/ensbiascoR/www/")

# 1. Define an observable meteorological metric as constraint: Tair
# 2. and 3. Estimate CDF's using e.g. a Gaussian kernel with reliable band-width estimation (Sheather and Jones, 1991):
obs.kernel = get.kernel.percentiles(data=ensbiascor.example1$obs.data$Tair[which(ensbiascor.example1$obs.data$Year %in% 1985:2010)])
mod.kernel = get.kernel.percentiles(data=ensbiascor.example1$mod.data$Tair)


# i. CDF plot:
png(filename="_example1_fig1.png", width=500, height=400)
par(mar=c(4,4,1,1))
plot(ecdf(obs.data$Tair[which(obs.data$Year %in% 1985:2010)]), bty='n', 
     xlab = "Observational metric", ylab="Cumulative Density Function", xlim=c(14,25), 
     main="", yaxt="n")
axis(side=2, las=1)
kernel.quantiles = seq(0, 1, 0.001)
lines(x = obs.kernel, y=kernel.quantiles, col="blue", lwd=2)
lines(x = mod.kernel, y=kernel.quantiles, col="red", lwd = 2)
dev.off()

# 4. Derive percentile transfer function between CDF's:
cur.fun = get.percentile.transfer.function(obs.grid.cell= obs.data$Tair, mod.grid.cell=mod.data$Tair)

# ii. percent-percent plot and resampling function:
png(filename="_example1_fig2.png", width=500, height=400)
par(mar=c(4,4,1,1))
plot(x = kernel.quantiles, y = cur.fun(kernel.quantiles), 
     xlab="Observations - CDF", ylab="Model Ensemble - CDF", bty="n",
     xlim = c(0,1), ylim=c(0, 1), type="n")
lines(x = kernel.quantiles, y = cur.fun(kernel.quantiles), col="red")
lines(x=c(-10^6, 10^6), y=c(-10^6,10^6), col="darkgray", lwd = 2)
# plot lines for translation:
sapply(X= seq(0, 1, by=0.1), FUN=function(x) lines(x = c(x, x), y = c(0, cur.fun(x)), lty=2, col="black"))
sapply(X= seq(0, 1, by=0.1), FUN=function(x) lines(x = c(0, x), y = c(cur.fun(x), cur.fun(x)), lty = 2, col="black"))
legend("topleft", c("Hermite Splines"), bty='n', lty=1, col="red")
dev.off()

# 5. Resample ensemble: Check quantile-quantile plot:

# resample ensemble WITH replacement:
resample.idx_replacement = resample.ensemble(mod.grid.cell = mod.data$Tair, transfer.function= cur.fun, 
                                 replacement = T, sample.size=10000)
# resampled ensemble WITHOUT replacement:
resample.idx_no_replacement = resample.ensemble(mod.grid.cell = mod.data$Tair, transfer.function= cur.fun, 
                                 replacement = F, sample.size=8000, search.radius = 0.5)



# iii. quantile-quantile plot:

# get plot range:
lim = c(floor(min(c(quantile(obs.kernel, probs=0.01), quantile(mod.kernel, probs=0.01)))), ceiling(max(c(quantile(obs.kernel, probs=0.99), quantile(mod.kernel, probs=0.99)))))

png(filename="_example1_fig3.png", width=500, height=400)
par(mar=c(4,4,1,1))
plot(c(1,1), type='n', xlim = lim, ylim= lim,
     xlab = "Observations", ylab = "Model Ensemble", bty='n')
lines(x=c(1, 100), y=c(1,100), col="darkblue", lwd = 2)
# plot.idx = sample(x=1:length(obs.kernel), size=100)
points(x = quantile(obs.kernel, seq(0.01, 0.99, 0.01)), y= quantile(mod.kernel, seq(0.01, 0.99, 0.01)), col="red", pch = 16, cex = 0.8)
points(x = quantile(obs.kernel, seq(0.01, 0.99, 0.01)), y= quantile(mod.data$Tair[resample.idx_replacement], seq(0.01, 0.99, 0.01)), col="black", pch = 16, cex = 0.8)
# points(x = quantile(obs.kernel, seq(0.01, 0.99, 0.01)), y= quantile(mod.data$Tair[resample.idx_no_replacement], seq(0.01, 0.99, 0.01), na.rm=T), col="brown", pch = 16, cex = 0.8)
legend("bottomright", c("Original ensemble", "Resampled ensemble"), col=c("red", "black"), bty='n', pch = 16)
dev.off()


# 6. Quality-control of resampled ensemble:
# percentile.ranges = c(0, 0.02, seq(0.1, 0.9, 0.1), 0.98, 1)
# text = c("< 2", "2-10", "10-20",       
#         "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-98", "> 98")

percentile.ranges = seq(0, 1, 0.1)
percentiles.areamean = quality.control.ensemble(transfer.function=cur.fun, percentile.ranges)

png(filename="_example1_fig4a.png", width=500, height=400)
par(mar=c(4,4,1,1))
text = c("< 10", "10-20",       
         "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "> 90")
ylim = c(0, 0.2)
barplot(percentiles.areamean, names = text, ylim = ylim,
        xlab = "Observations - Percentiles", ylab = "Fraction of original ensemble members in each percentile bin", 
        cex.lab = 0.7, cex.axis=0.7, cex.names = 0.7)
dev.off()

# make normalized version:
png(filename="_example1_fig4b.png", width=500, height=400)
par(mar=c(4,4,1,1))
ylim = c(0, 2)
barplot(percentiles.areamean / diff(percentile.ranges), names = text, ylim = ylim,
        xlab = "Observations - Percentiles", ylab = "Relative change in sample size")
dev.off()


# 7. Analyze resampled ensemble:
xlim = c(floor(min(c(quantile(obs.kernel, probs=0.01), quantile(mod.kernel, probs=0.01)))), ceiling(max(c(quantile(obs.kernel, probs=0.99), quantile(mod.kernel, probs=0.99)))))
year.idx = which(obs.data$Year %in% 1985:2010)
ylim = c(min(c(obs.data$Precip[which(obs.data$Year %in% 1985:2010)], quantile(mod.data$Precip, 0.01))), max(c(obs.data$Precip[which(obs.data$Year %in% 1985:2010)], quantile(mod.data$Precip, 0.999))))


png(filename="_example1_fig5.png", width=500, height=400)
par(mar=c(4,4,1,1))
plot(c(1,1), xlim = xlim, ylim=ylim,
     xlab = "Temperature [°C]", ylab = "Precipitation [mm]", bty="n")
# model data: 
points(x=mod.data$Tair, y=mod.data$Precip, col = make.transparent.color("darkred", alpha=25), pch = 16)
points(x=mod.data$Tair[resample.idx_no_replacement], y=mod.data$Precip[resample.idx_no_replacement], col = make.transparent.color("darkblue", alpha=25), pch=16)

# observations data:
points(x=obs.data$Tair[year.idx], y=obs.data$Precip[year.idx], col = "black", pch = 16)
legend("topright", c("Original Ensemble", "Resampled Ensemble", "Observations"), pch = 16, col=c("darkred", "darkblue", "black"), bty='n')
legend("bottomleft", c(paste("R_ens-orig = ", round(cor(mod.data$Tair, mod.data$Precip), 2), sep=""),
                       paste("R_ens-cor = ", round(cor(mod.data$Tair[resample.idx_replacement], mod.data$Precip[resample.idx_replacement]), 2), sep=""),
                       paste("R_obs = ", round(cor(obs.data$Tair[year.idx], obs.data$Precip[year.idx]), 2), sep="")), text.col=c("darkred", "darkblue", "black"), bty="n")
dev.off()



# Which plots are essential?
# ---------------------------

# 1. check biases: CDF comparison
# 2. check resampling function: percent-percent plot
# 3. check resampled ensemble: qq plot of original ensemble
# 4. check sample size of resampled ensemble
# 5. plot constraint in original and resampled ensemble against obs and mod!



