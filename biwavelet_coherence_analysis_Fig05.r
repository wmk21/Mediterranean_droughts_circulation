#-----------------------
#   Wavelet coherence analysis
#   Wavelet plot  until 1850. 
#-----------------------

fval <- 100    # Moving window values: Nevents/fval; if the loaded time series is the number of occurrences.

rm=(list=ls())

library(ncdf4)
library(biwavelet)

maindir <- ""

"""
   Read time series data for each or all model simulations, for each of the regions.
   Time series of SOIL or of Number of occurrence
"""

freq_east_data <- nc_open(paste(maindir,"model.frequency.east.",fval,".0850-2005.nc",sep=""))

freq_west_data <- nc_open(paste(maindir,"model.frequency.west.",fval,".0850-2005.nc",sep=""))

cesm_freq_east_data <- nc_open(paste(maindir,"cesm.frequency.east.",fval,".0850-2005.nc",sep=""))

cesm_freq_west_data <- nc_open(paste(maindir,"cesm.frequency.west.",fval,".0850-2005.nc",sep=""))

model_name <- ncvar_get(freq_east_data, "model")


freq_east <- ncvar_get(freq_east_data, "__xarray_dataarray_variable__")
freq_west <- ncvar_get(freq_west_data, "__xarray_dataarray_variable__")

cesm_freq_east <- ncvar_get(cesm_freq_east_data, "__xarray_dataarray_variable__")
cesm_freq_west <- ncvar_get(cesm_freq_west_data, "__xarray_dataarray_variable__")


#-----------------------
#   Standardizing the time series
#-----------------------

giss_soil_e <- (freq_east[1, 1,] - mean(freq_east[1, 1,]))/sd(freq_east[1, 1,])
giss_soil_w <- (freq_west[1, 1, ] - mean(freq_west[1, 1,]))/sd(freq_west[1, 1,])

ccsm_soil_e <- (freq_east[2, 1,] - mean(freq_east[2, 1,]))/sd(freq_east[2, 1,])
ccsm_soil_w <- (freq_west[2, 1,] - mean(freq_west[2, 1,]))/sd(freq_west[2, 1, ])

bcc_soil_e <- (freq_east[3, 1,] - mean(freq_east[3, 1,]))/sd(freq_east[3, 1,])
bcc_soil_w <- (freq_west[3, 1,] - mean(freq_west[3, 1,]))/sd(freq_west[3, 1, ])

miroc_soil_e <- (freq_east[4, 1,] - mean(freq_east[4, 1,]))/sd(freq_east[4, 1,])
miroc_soil_w <- (freq_west[4, 1,] - mean(freq_west[4, 1,]))/sd(freq_west[4, 1, ])

cesm_soil_e <- (cesm_freq_east[1,] - mean(cesm_freq_east[1,]))/sd(cesm_freq_east[1,])

cesm_soil_w <- (cesm_freq_west[1,] - mean(cesm_freq_west[1,]))/sd(cesm_freq_west[1, ])


years <- seq(850, 2005)



bcc_soil_ef <- cbind(years, bcc_soil_e)
bcc_soil_wf <- cbind(years, bcc_soil_w)

cesm_soil_ef <- cbind(years, cesm_soil_e)
cesm_soil_wf <- cbind(years, cesm_soil_w)

ccsm_soil_ef <- cbind(years, ccsm_soil_e)
ccsm_soil_wf <- cbind(years, ccsm_soil_w)

miroc_soil_ef <- cbind(years, miroc_soil_e)
miroc_soil_wf <- cbind(years, miroc_soil_w)

giss_soil_ef <- cbind(years, giss_soil_e)
giss_soil_wf <- cbind(years, giss_soil_w)

nrands = 500
m=128

#--- cesm
wtc.AB = wtc(cesm_soil_wf, cesm_soil_ef, max.scale =m, nrands = nrands)

#--- giss
wtc.CD = wtc(giss_soil_wf, giss_soil_ef, max.scale =m, nrands = nrands)

#--- ccsm
wtc.EF = wtc(ccsm_soil_wf, ccsm_soil_ef, max.scale =m, nrands = nrands)

#--- bcc
wtc.HI = wtc(bcc_soil_wf, bcc_soil_ef, max.scale =m, nrands = nrands)

#--- miroc
wtc.JK = wtc(miroc_soil_wf, miroc_soil_ef, max.scale =m, nrands = nrands)

# CESM
dev.new(width=8.5, height=4, unit="in")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 7) + 0.1)

plot(wtc.AB, ncol=32,plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, lwd.sig = 2, arrow.lwd = 0.05, arrow.len = 0.1, ylab = "Period", xlab = "", plot.cb = TRUE, main = "CESM west - east",xaxt="n")

axis(side = 1, at = c(seq(900, 2000, 100)),labels=c(seq(900, 2000, 100)))

dev.copy2pdf(file = paste(maindir,"wtc_freq",fval,"_cesm_0850-2005.pdf", sep=""))

# GISS

dev.new(width=8.5, height=4, unit="in")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 7) + 0.1)

plot(wtc.CD, ncol=32,plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, lwd.sig = 2, arrow.lwd = 0.05, arrow.len = 0.1, ylab = "Period", xlab = "", plot.cb = TRUE, main = "GISS west - east",xaxt="n")

axis(side = 1, at = c(seq(900, 2000, 100)),labels=c(seq(900, 2000, 100)))

dev.copy2pdf(file = paste(maindir,"wtc_freq",fval,"_giss_0850-2005.pdf", sep=""))


# CCSM

dev.new(width=8.5, height=4, unit="in")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 7) + 0.1)

plot(wtc.EF, ncol=32,plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, lwd.sig = 2, arrow.lwd = 0.05, arrow.len = 0.1, ylab = "Period", xlab = "", plot.cb = TRUE, main = "CCSM west - east",xaxt="n")

axis(side = 1, at = c(seq(900, 2000, 100)),labels=c(seq(900, 1850, 100)))

dev.copy2pdf(file = paste(maindir,"wtc_freq",fval,"_ccsm_0850-2005.pdf", sep=""))


# BCC

dev.new(width=8.5, height=4, unit="in")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 7) + 0.1)

plot(wtc.HI, ncol=32,plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, lwd.sig = 2, arrow.lwd = 0.05, arrow.len = 0.1, ylab = "Period", xlab = "", plot.cb = TRUE, main = "BCC west - east",xaxt="n")

axis(side = 1, at = c(seq(900, 2000, 100)),labels=c(seq(900, 2000, 100)))

dev.copy2pdf(file = paste(maindir,"wtc_freq",fval,"_bcc_0850-2005.pdf", sep=""))


# MIROC

dev.new(width=8.5, height=4, unit="in")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 7) + 0.1)

plot(wtc.JK, ncol=32,plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, lwd.sig = 2, arrow.lwd = 0.05, arrow.len = 0.1, ylab = "Period", xlab = "", plot.cb = TRUE, main = "MIROC west - east",xaxt="n")

axis(side = 1, at = c(seq(900, 2000, 100)),labels=c(seq(900, 2000, 100)))

dev.copy2pdf(file = paste(maindir,"wtc_freq",fval,"miroc_0850-2005.pdf", sep=""))

