# slv_tapwater_source_model.R

# Code description, and rough timeline of updates:
# ================================================
# 30nov16 - added multivariate estimation using d18O, d2H, and 4 
# mixtures of the two largest source volumes.
#
#
#
#
#
#
#
#

# remove old variables before running script, good practice for debugging
rm(list=ls())

###########################################################################
###########################################################################
# LOAD NECESSARY PACKAGES
###########################################################################
############################################################################

# packages for multivariate normal distribution generation
library(mvtnorm)

# packages for spatial data analysis
library(rgdal)
library(maptools)
library(ggplot2)
library(plyr)
library(rgeos)

###########################################################################
###########################################################################
# LOAD DATA, RESTRUCTURE FOR ANALYSIS
###########################################################################
###########################################################################

# Yusuf's computer
#setwd("C:/Users/u0817369/Dropbox/Utah_water_isotopes/JVWCD")
#setwd("C:/Users/yusuf/Dropbox/Utah_water_isotopes/JVWCD")
# Rich's computer
setwd("~/Documents/SALT-LAKE-TAP-WATER-MODEL/")
data <- read.csv("May_5_nomix.csv",header=TRUE)

# break into source and supplyline data frames
source.data <- subset(data,data$Site_Description=="Source")
supply.data <- subset(data,data$Site_Description=="Supplyline")

# turn off warnings quickly while setting up the data (data.frame issues
# a few relatively unhelpful warnings about the length of the output) -
# they can safely be disregarded.
options(warn=-1)

#-----------------------------------------------------
# structure source data frame first - some sites were 
# sampled more than once, we've chosen here to 
# average multiple measurements of individual sites.
#-----------------------------------------------------

# find unique site IDs
source.siteIDs <- unique(source.data$Site_ID)

# initialize data vectors
source.mean_d18O <- vector()
source.mean_d2H <- vector()
source.mean_DEX <- vector() 
source.stdev_d18O <- vector()
source.stdev_d2H <- vector()
source.stdev_DEX <- vector()
source.month <- vector()
source.year <- vector()
source.site_ID <- vector()
source.lat <- vector()
source.lon <- vector()
source.uniqueID <-vector()
source.address <-vector()

# loop across unique site IDs
for (i in 1:length(source.siteIDs)) {
	# calculate mean isotope values
	source.mean_d18O[i] <- mean(subset(source.data$d18O,source.data$Site_ID==source.siteIDs[i]))
	source.mean_d2H[i] <- mean(subset(source.data$d2H,source.data$Site_ID==source.siteIDs[i]))
	source.mean_DEX[i] <- mean(subset(source.data$DEX,source.data$Site_ID==source.siteIDs[i]))
	# propagate uncertainties in new means - CHECK THE FORMULA HERE
	source.stdev_d18O[i] <- sqrt(sum(subset(source.data$d18O_sd,source.data$Site_ID==source.siteIDs[i])^2))
	source.stdev_d2H[i] <- sqrt(sum(subset(source.data$d2H_sd,source.data$Site_ID==source.siteIDs[i])^2))
	source.stdev_DEX[i] <- sqrt(sum(subset(source.data$DEX_sd,source.data$Site_ID==source.siteIDs[i])^2))
	# pull out coordinate variables, metadata
	source.month[i] <- subset(source.data$MONTH,source.data$Site_ID==source.siteIDs[i])
	source.year[i] <- subset(source.data$YEAR,source.data$Site_ID==source.siteIDs[i])[1]
	source.address[i] <- subset(source.data$Address,source.data$Site_ID==source.siteIDs[i])[1]
	source.lat[i] <- subset(source.data$Latitude,source.data$Site_ID==source.siteIDs[i])[1]
	source.lon[i] <- subset(source.data$Longitude,source.data$Site_ID==source.siteIDs[i])[1]
}
	# assign a unique number to site IDs - this line makes a sequence starting at 1, and goes to the number of source locations.
	source.uniqueID <- seq(1,length(source.lon),1)

	# volumes of different treatment plants - this will inform our Bayesian priors!
	source.volumes <- c(150,300,300,3000,1200)

# combine vectors back into single, reduced dataframe - STDEV DEX missing!
source.reduced <- data.frame("mean_d18O"=source.mean_d18O,"mean_d2H"=source.mean_d2H,"mean_DEX"=source.mean_DEX,
	"stdev_d18O"=source.stdev_d18O,"stdev_d2H"=source.stdev_d2H,"stdev_DEX"=source.stdev_DEX,"month"=source.month,"year"=source.year,
	"lat"=source.lat,"lon"=source.lon, "SITE_ID" = source.siteIDs,"uniqueIDS" = source.uniqueID,"volumes"=source.volumes)

#-----------------------------------------------------
# add 4 mixtures of sources 4 and 5.
#-----------------------------------------------------

# mixture isotope means
source.mixtures.d18O = c(20,40,60,80)*source.reduced$mean_d18O[4]/100 + c(80,60,40,20)*source.reduced$mean_d18O[5]/100
source.mixtures.d2H = c(20,40,60,80)*source.reduced$mean_d2H[4]/100 + c(80,60,40,20)*source.reduced$mean_d2H[5]/100
source.mixtures.DEX = 8*source.mixtures.d18O - source.mixtures.d2H

# mixture isotope uncertainties
source.mixtures.d18O_sd = rep(0.2,4)
source.mixtures.d2H_sd = rep(0.6,4)
source.mixtures.DEX_sd = rep(NA,4)

# mixture metadata
source.mixtures.month = rep(source.month[1],4)
source.mixtures.year = rep(source.year[1],4)
source.mixtures.address = rep(NA,4)
source.mixtures.lat = rep(NA,4)
source.mixtures.lon = rep(NA,4)

# mixture unique IDs
source.mixtures.uniqueID <- seq(length(source.lon)+1,length(source.lon)+4,1)

# mixture volumes (doesn't really make physical sense)
source.mixtures.volumes <- rep(1000,4)

# site IDS
source.mixtures.siteIDs <- c("4.80M5.20","4.60M5.40","4.40M5.60","4.20M5.80")

# combine vectors back into single, reduced dataframe - STDEV DEX missing!
source.mixtures <- data.frame("mean_d18O"=source.mixtures.d18O,"mean_d2H"=source.mixtures.d2H,"mean_DEX"=source.mixtures.DEX,
	"stdev_d18O"=source.mixtures.d18O_sd,"stdev_d2H"=source.mixtures.d2H_sd,"stdev_DEX"=source.mixtures.DEX_sd,"month"=source.mixtures.month,"year"=source.mixtures.year,
	"lat"=source.mixtures.lat,"lon"=source.mixtures.lon, "SITE_ID" = source.mixtures.siteIDs,"uniqueIDS" = source.mixtures.uniqueID,"volumes"=source.mixtures.volumes)

# combine mixture dataframe with "raw" dataframe
source.reduced <- rbind(source.reduced,source.mixtures)

#-----------------------------------------------------
# now structure the sampling location data frame
#-----------------------------------------------------

# find unique site IDs
supply.siteIDs <- unique(supply.data$Site_ID)

# initialize data vectors
supply.mean_d18O <- vector()
supply.mean_d2H <- vector()
supply.mean_DEX <- vector() 
supply.stdev_d18O <- vector()
supply.stdev_d2H <- vector()
supply.stdev_DEX <- vector()
supply.month <- vector()
supply.year <- vector()
supply.address <- vector()
supply.lat <- vector()
supply.lon <- vector()

# loop across unique site IDs
for (i in 1:length(supply.siteIDs)) {
	# calculate mean isotope values
	supply.mean_d18O[i] <- mean(subset(supply.data$d18O,supply.data$Site_ID==supply.siteIDs[i]))
	supply.mean_d2H[i] <- mean(subset(supply.data$d2H,supply.data$Site_ID==supply.siteIDs[i]))
	supply.mean_DEX[i] <- mean(subset(supply.data$DEX,supply.data$Site_ID==supply.siteIDs[i]))
	# propagate uncertainties in new means - CHECK THE FORMULA HERE
	supply.stdev_d18O[i] <- sqrt(sum(subset(supply.data$d18O_sd,supply.data$Site_ID==supply.siteIDs[i])^2))
	supply.stdev_d2H[i] <- sqrt(sum(subset(supply.data$d2H_sd,supply.data$Site_ID==supply.siteIDs[i])^2))
	supply.stdev_DEX[i] <- sqrt(sum(subset(supply.data$DEX_sd,supply.data$Site_ID==supply.siteIDs[i])^2))
	# pull out coordinate variables, metadata
	supply.month[i] <- subset(supply.data$MONTH,supply.data$Site_ID==supply.siteIDs[i])[1]
	supply.year[i] <- subset(supply.data$YEAR,supply.data$Site_ID==supply.siteIDs[i])[1]
	#supply.[i] <- subset(supply.data$Address,supply.data$Site_ID==supply.siteIDs[i])[1]
	supply.lat[i] <- subset(supply.data$Latitude,supply.data$Site_ID==supply.siteIDs[i])[1]
	supply.lon[i] <- subset(supply.data$Longitude,supply.data$Site_ID==supply.siteIDs[i])[1]
}

# combine vectors back into single, reduced dataframe - STDEV DEX missing!
supply.reduced <- data.frame("mean_d18O"=supply.mean_d18O,"mean_d2H"=supply.mean_d2H,"mean_DEX"=supply.mean_DEX,
	"stdev_d18O"=supply.stdev_d18O,"stdev_d2H"=supply.stdev_d2H,"stdev_DEX"=supply.stdev_DEX,"month"=supply.month,"year"=supply.year,
	"lat"=supply.lat,"lon"=supply.lon, "SITE_ID" = supply.siteIDs)

# reinstate warnings beyond this point.
options(warn=0)

###########################################################################
###########################################################################
# INITIAL DIAGNOSTIC PLOTS - d18O vs d2H, d2H vs dex, d18O vs dex
###########################################################################
###########################################################################

quartz()
par(mfrow=c(1,3))

# identify min/max values for all three parameters, use these for plot windows
minO = min(c(min(supply.mean_d18O),min(source.mean_d18O)))
maxO = max(c(max(supply.mean_d18O),max(source.mean_d18O)))

minH = min(c(min(supply.mean_d2H),min(source.mean_d2H)))
maxH = max(c(max(supply.mean_d2H),max(source.mean_d2H)))

mind = min(c(min(supply.mean_DEX),min(source.mean_DEX)))
maxd = max(c(max(supply.mean_DEX),max(source.mean_DEX)))

# d18O vs d2H
#-------------
# plot source data
plot(source.mean_d18O,source.mean_d2H,ylim=c(minH-0.1*(maxH-minH),maxH+0.1*(maxH-minH)),
	xlim=c(minO-0.1*(maxO-minO),maxO+0.1*(maxO-minO)),pch=15,col="red")
# plot supply data
points(supply.mean_d18O,supply.mean_d2H,pch=18,col="blue")
# plot source errors
arrows(source.mean_d18O,source.mean_d2H+source.stdev_d2H,source.mean_d18O,source.mean_d2H-source.stdev_d2H,
	code=3,angle=90,length=0,col="red")
arrows(source.mean_d18O-source.stdev_d18O,source.mean_d2H,source.mean_d18O+source.stdev_d18O,source.mean_d2H,
	code=3,angle=90,length=0,col="red")
# plot supply errors
arrows(supply.mean_d18O,supply.mean_d2H+supply.stdev_d2H,supply.mean_d18O,supply.mean_d2H-supply.stdev_d2H,
	code=3,angle=90,length=0,col="blue")
arrows(supply.mean_d18O-supply.stdev_d18O,supply.mean_d2H,supply.mean_d18O+supply.stdev_d18O,supply.mean_d2H,
	code=3,angle=90,length=0,col="blue")

# d2H vs DEX
#-------------
# plot source data
plot(source.mean_d2H,source.mean_DEX,ylim=c(mind-0.1*(maxd-mind),maxd+0.1*(maxd-mind)),
	xlim=c(minH-0.1*(maxH-minH),maxH+0.1*(maxH-minH)),pch=15,col="red")
# plot supply errors
points(supply.mean_d2H,supply.mean_DEX,pch=18,col="blue")
# plot source errors
arrows(source.mean_d2H,source.mean_DEX+source.stdev_DEX,source.mean_d2H,source.mean_DEX-source.stdev_DEX,
	code=3,angle=90,length=0,col="red")
arrows(source.mean_d2H-source.stdev_d2H,source.mean_DEX,source.mean_d2H+source.stdev_d2H,source.mean_DEX,
	code=3,angle=90,length=0,col="red")
# plot supply errors
arrows(supply.mean_d2H,supply.mean_DEX+supply.stdev_DEX,supply.mean_d2H,supply.mean_DEX-supply.stdev_DEX,
	code=3,angle=90,length=0,col="blue")
arrows(supply.mean_d2H-supply.stdev_d2H,supply.mean_DEX,supply.mean_d2H+supply.stdev_d2H,supply.mean_DEX,
	code=3,angle=90,length=0,col="blue")

# d18O vs DEX
#-------------
# plot source data
plot(source.mean_d18O,source.mean_DEX,ylim=c(mind-0.1*(maxd-mind),maxd+0.1*(maxd-mind)),
	xlim=c(minO-0.1*(maxO-minO),maxO+0.1*(maxO-minO)),pch=15,col="red")
# plot supply data
points(supply.mean_d18O,supply.mean_DEX,pch=18,col="blue")
# plot source errors
arrows(source.mean_d18O,source.mean_DEX+source.stdev_DEX,source.mean_d18O,source.mean_DEX-source.stdev_DEX,
	code=3,angle=90,length=0,col="red")
arrows(source.mean_d18O-source.stdev_d18O,source.mean_DEX,source.mean_d18O+source.stdev_d18O,source.mean_DEX,
	code=3,angle=90,length=0,col="red")
# plot supply errors
arrows(supply.mean_d18O,supply.mean_DEX+supply.stdev_DEX,supply.mean_d18O,supply.mean_DEX-supply.stdev_DEX,
	code=3,angle=90,length=0,col="blue")
arrows(supply.mean_d18O-supply.stdev_d18O,supply.mean_DEX,supply.mean_d18O+supply.stdev_d18O,supply.mean_DEX,
	code=3,angle=90,length=0,col="blue")

###########################################################################
###########################################################################
# PERFORM BAYESIAN ANALYSIS - this section needs a lot more notes about what we've done!
###########################################################################
###########################################################################

# bivariate using d18O and d2H
lvals_H <- matrix(nrow=nrow(supply.reduced),ncol=nrow(source.reduced))

sigma_d18O <- 0.2
sigma_d2H <- 0.5
rho <- cor(data$d18O,data$d2H) # hardcode to the actual dataset being analyzed.

corr_matrix <- matrix(c(sigma_d18O^2, sigma_d18O*sigma_d2H*rho, sigma_d18O*sigma_d2H*rho, sigma_d2H^2),2)

upper <- matrix(nrow=nrow(supply.reduced),ncol=2)
lower <- matrix(nrow=nrow(supply.reduced),ncol=2)
mean <- matrix(nrow=nrow(source.reduced),ncol=2)
for (i in 1:nrow(supply.reduced)) {
  for (j in 1:nrow(source.reduced)) {
    # assign likelihood for given i,j
    upper[i,] = c((supply.reduced$mean_d2H[i] + 3*supply.reduced$stdev_d2H[i]), (supply.reduced$mean_d18O[i] + 3*supply.reduced$stdev_d18O[i]))
    lower[i,] = c((supply.reduced$mean_d2H[i] - 3*supply.reduced$stdev_d2H[i]), (supply.reduced$mean_d18O[i] - 3*supply.reduced$stdev_d18O[i]))
    mean[j,] = c((source.reduced$mean_d2H[j]),(source.reduced$mean_d18O[j]))
    
    lvals_H[i,j] <- pmvnorm(upper = upper[i,], lower = lower[i,], mean = mean[j,], sigma = corr_matrix)
  }
}

# due to monte carlo simulation and stochastic nature of pmvnorm, some likelihoods are 
# negative. set a floor here to require values less than zero to be zero
lvals_H[lvals_H<0] <- 0.0

colnames(lvals_H) <- c("solena_way","JVWTP","SEWTP","SWWTP","WELL_1300","4.80M5.20","4.60M5.40","4.40M5.60","4.20M5.80")
lH<-data.frame(lvals_H)

post_H_uniform <- matrix(nrow=nrow(lvals_H), ncol=ncol(lvals_H))
post_H_volwgt <- matrix(nrow=nrow(lvals_H), ncol=ncol(lvals_H))

i <- 1:nrow(lvals_H)
j <- 1:ncol(lvals_H)
# assign posterior for given i,j supply.
for (j in 1:nrow(source.reduced)) {
	post_H_uniform[i,j] <- (lvals_H[i,j]*(1/ncol(lvals_H)))/((1/ncol(lvals_H))*rowSums(lvals_H))
	post_H_volwgt[i,j] <- (lvals_H[i,j]*(source.reduced$volume[j]/sum(source.reduced$volume,na.rm=T)))/
		((source.reduced$volume[j]/sum(source.reduced$volume,na.rm=T))*rowSums(lvals_H))
}

# plotting the POSTERIOR values

quartz()
pp_H <- plot(x= seq(1,9,1),post_H_volwgt[1,], col ="red", type = "l", xlim = c(0.8,9.1), ylim = c(0,1) )
title(main = "d2H_POSTERIOR")
p = sample(rainbow(63))
for (i in 2:64){
  lines(x= seq(1,9,1),post_H_volwgt[i,], type = "l", col = p[i])
}

pp_export <- data.frame(post_H_volwgt, supply.reduced$lat, supply.reduced$lon, supply.reduced$SITE_ID)

# write.csv(pp_export, file ="post_may2015_H.csv")

###########################################################################
###########################################################################
# PLOT DATA AND LIKELIHOODS ON SPATIAL MAPS - this section still needs some cleaning.
###########################################################################
###########################################################################

# load the salt lake valley shapefile.
shapefile <- readOGR(dsn=path.expand("~/Documents/SALT-LAKE-TAP-WATER-MODEL/slc_tap_shape_file") ,layer = "Export_Output")
shapefile@data$id <- rownames(shapefile@data)
shapefile1 <- fortify(shapefile, region = "id")
#shapefileDF <- merge(shapefile1, shapefile@data, by = "id")

#y = ggplot(shapefileDF, aes(long,lat, group = group, fill = "WR_ID")) + geom_polygon() + geom_path(color = "black") + coord_equal() + scale_fill_manual(values=c("white"), guide="none")# + coord_map("albers", lat0=40, lat1=40) #+ (color = 'gray', fill = 'blue', size = 2) 
#y

#coordinates(pp_export)<-~supply.reduced.lon+supply.reduced.lat
#proj4string(pp_export)
#proj4string(pp_export)<-CRS("+proj=longlat +datum=NAD83")
#pp_export1<-spTransform(pp_export, CRS(proj4string(shapefile)))
#identical(proj4string(pp_export1),proj4string(shapefile))
#pp_export1<-data.frame(pp_export1)

# source.reduced1<-subset(source.reduced[1:5,])
# #coordinates(source.reduced1)<-~lon+lat
# proj4string(source.reduced1)
# proj4string(source.reduced1)<-CRS("+proj=longlat +datum=NAD83")
# source.reduced1<-spTransform(source.reduced1, CRS(proj4string(shapefile)))
# #identical(proj4string(pp_export1),proj4string(shapefile))
# source.reduced1<-data.frame(source.reduced1)

quartz()
x =  ggplot() + geom_point(data = pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X1, color = X1), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced$lon[1], y= source.reduced$lat[1]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
#x = x + geom_point(aes(x=source.reduced1$lon[4], y= source.reduced1$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
print(x)

quartz()
x =  ggplot() + geom_point(data = pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X2, color = X2), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced$lon[2], y= source.reduced$lat[2]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
#x = x + geom_point(aes(x=source.reduced1$lon[4], y= source.reduced1$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
print(x)

quartz()
x =  ggplot() + geom_point(data = pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X3, color = X3), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced$lon[3], y= source.reduced$lat[3]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
#x = x + geom_point(aes(x=source.reduced1$lon[4], y= source.reduced1$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
print(x)

quartz()
x =  ggplot() + geom_point(data = pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X4, color = X4), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced$lon[4], y= source.reduced$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
#x = x + geom_point(aes(x=source.reduced1$lon[4], y= source.reduced1$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
print(x)

quartz()
x =  ggplot() + geom_point(data = pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X5, color = X5), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced$lon[5], y= source.reduced$lat[5]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
#x = x + geom_point(aes(x=source.reduced1$lon[4], y= source.reduced1$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
print(x)

quartz()
x =  ggplot() + geom_point(data = pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X6, color = X6), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced$lon[5], y= source.reduced$lat[5]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
x = x + geom_point(aes(x=source.reduced$lon[4], y= source.reduced$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
print(x)

quartz()
x =  ggplot() + geom_point(data = pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X7, color = X7), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced$lon[5], y= source.reduced$lat[5]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
x = x + geom_point(aes(x=source.reduced$lon[4], y= source.reduced$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
print(x)

quartz()
x =  ggplot() + geom_point(data = pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X8, color = X8), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced$lon[5], y= source.reduced$lat[5]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
x = x + geom_point(aes(x=source.reduced$lon[4], y= source.reduced$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
print(x)

quartz()
x =  ggplot() + geom_point(data = pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X9, color = X9), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced$lon[5], y= source.reduced$lat[5]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
x = x + geom_point(aes(x=source.reduced$lon[4], y= source.reduced$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
print(x)

# rel_likelihood_H <- matrix(nrow=nrow(post),ncol=ncol(post))

# for (i in 1:nrow(post_H)) {
#   rel_likelihood_H[i,] <- post_H[i,]/max(post_H[i,]) 
# }


# ##############################################
# # calculation of likelihood values for DEX
# ##############################################

# lvals_DEX <- matrix(nrow=nrow(supply.reduced),ncol=nrow(source.reduced))
# for (i in 1:nrow(supply.reduced)) {
#   for (j in 1:nrow(source.reduced)) {
#     # assign likelihood for given i,j
#     lvals_DEX[i,j] <- pnorm((supply.reduced$mean_DEX[i] + supply.reduced$stdev_DEX[i]), 
#                         mean = source.reduced$mean_DEX[j], sd = source.reduced$stdev_DEX[j], lower.tail = TRUE) - 
#       pnorm((supply.reduced$mean_DEX[i] - supply.reduced$stdev_DEX[i]), 
#             mean = source.reduced$mean_DEX[j], sd = source.reduced$stdev_DEX[j], lower.tail = TRUE)
#   }
# }
# colnames(lvals_DEX) <- c("solena_way", "JVWTP", "SEWTP", "SWWTP" , "WELL_1300" ,  "M1", "M2" )
# l_DEX<-data.frame(lvals_DEX)

# ###########################
# # PLOTTING THE LIKELIHOOD FOR dex
# ################
# lvals_DEX_df <- data.frame("SITE_ID" = supply.reduced$SITE_ID, "solena_way" = l_DEX$solena_way, "JVWTP" = l_DEX$JVWTP,
#                            "SEWTP" =l_DEX$SEWTP, "SWWTP" =l_DEX$SWWTP,"WELL_1300" = l_DEX$WELL_1300, "M1" = l_DEX$M1, "M2" = l_DEX$M2  )
# ll_DEX <- plot(x= c(1,2,3,4,5,6,7),lvals_DEX_df[1,2:8], col ="red", type = "b", xlim = c(0.8,7.1), ylim = c(0,1) )
# title(main = "D-EXCESS_LIKELIHOOD")
# p = sample(rainbow(63))
# for (i in 2:64){
#   lines(x= c(1,2,3,4,5,6,7),lvals_DEX_df[i,2:8], type = "b", col = p[i])
# }

# #calculation of posterior distributions
# post_DEX <- matrix(nrow=nrow(lvals_DEX), ncol=ncol(lvals_DEX))
# i <- 1:nrow(lvals_DEX)
# j <- 1:ncol(lvals_DEX)
# # assign posterior for given i,jsupply.
# post_DEX[i,j] <- (lvals_DEX[i,j]*ncol(lvals_DEX))/(ncol(lvals_DEX)* rowSums(lvals_DEX))
# #rowSums(post)
# #colSums(post)
# #########################################
# ##plotting the POSTERIOR values
# pp_DEX <- plot(x= c(1,2,3,4,5,6,7),post_DEX[1,1:7], col ="red", type = "l", xlim = c(0.8,7.1), ylim = c(0,1) )
# title(main = "DEX_POSTERIOR")
# p = sample(rainbow(63))
# for (i in 2:64){
#   lines(x= c(1,2,3,4,5,6,7),post_DEX[i,1:7], type = "l", col = p[i])
# }

# pp_DEX_exp <- data.frame(post_DEX, supply.reduced$lat, supply.reduced$lon, supply.reduced$SITE_ID)
# #write.csv(pp_export, file ="post_may2015_H.csv")


# library(ggplot2)
# x = ggplot(pp_export, aes(x=supply.reduced.lon, y=supply.reduced.lat, color=X3)) + geom_point() 
# x = x + ggtitle("d2H") + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
# x

# ###############################################################
# ###############################################################
# ##############3xporting data with supply and corresponding sources d18O
# ####################################################################
# export_source_O <- matrix(nrow=nrow(post)+4, ncol=ncol(post)+4)
# i<-5:nrow(export_source_O)
# a<- 1:nrow(rel_likelihood)
# b<- 1:ncol(rel_likelihood)
# export_source_O[i,1] <- supply.reduced[,10]
# export_source_O[i,2] <- supply.reduced[,1]
# export_source_O[i,3] <- supply.reduced[,8]
# export_source_O[i,4] <- supply.reduced[,9]
# j<-5:ncol(export_source_O)
# export_source_O[i,j] <-rel_likelihood[a,b]
# export_source_O[1,j] <- source.reduced[,1]
# export_source_O[2,j] <- source.reduced[,8]
# export_source_O[3,j] <- source.reduced[,9]
# export_source_O[4,j] <-j-4
# export_source_O
# ############################
# ##########################################
# ###############source and supply matrix combined
# export_2 <- matrix(nrow=64,ncol=6)
# i <- 1:nrow(export_2)
# export_2[i,1] <- export_source_O[i+4,2]
# export_2[i,2] <- export_source_O[i+4,3]
# export_2[i,3] <- export_source_O[i+4,4]
# q <-vector()
# lat<-vector()
# long<-vector()
# for (b in 5:68) {
#   for (a in 5:23) {
#     if (export_source_O[b,a] == 1) {
#      q[b-4] <- export_source_O[1,a]
#      lat[b-4] <-export_source_O[2,a]
#      long[b-4] <-export_source_O[3,a]
     
#     } 
#  }
# }
# export_2[i,4] <- q[i]
# export_2[i,5] <-lat[i]
# export_2[i,6] <-long[i]
# export_2
# ss_month <-data.frame ("meansupplyd18O" = export_2[i,1], "meansupplylat" = export_2[i,2], "meansupplylong" = export_2[i,3],
#                        "meansourced18O" = export_2[i,4], "meansourcelat" = export_2[i,5], "meansourcelong" = export_2[i,6])

# ss_month2 <- transform(ss_month, id=match(meansourced18O, unique(meansourced18O)))
# all.equal(ss_month,ss_month2)
# ss_month2
# write.csv(ss_month2, file ="supply_source_may2015_O.csv")
# ####################################
# ###############################################################
# ##############3xporting data with supply and corresponding sources d2H
# ####################################################################
# export_source_H <- matrix(nrow=nrow(post)+4, ncol=ncol(post)+4)
# i<-5:nrow(export_source_H)
# a<- 1:nrow(rel_likelihood_H)
# b<- 1:ncol(rel_likelihood_H)
# export_source_H[i,1] <- supply.reduced[,10]
# export_source_H[i,2] <- supply.reduced[,2]
# export_source_H[i,3] <- supply.reduced[,8]
# export_source_H[i,4] <- supply.reduced[,9]
# j<-5:ncol(export_source_H)
# export_source_H[i,j] <-rel_likelihood_H[a,b]
# export_source_H[1,j] <- source.reduced[,2]
# export_source_H[2,j] <- source.reduced[,8]
# export_source_H[3,j] <- source.reduced[,9]
# export_source_H[4,j] <-j-4
# export_source_H
# ############################
# ##########################################
# ###############source and supply matrix combined
# export_2 <- matrix(nrow=64,ncol=6)
# i <- 1:nrow(export_2)
# export_2[i,1] <- export_source_H[i+4,2]
# export_2[i,2] <- export_source_H[i+4,3]
# export_2[i,3] <- export_source_H[i+4,4]
# q <-vector()
# lat<-vector()
# long<-vector()
# for (b in 5:68) {
#   for (a in 5:23) {
#     if (export_source_H[b,a] == 1) {
#       q[b-4] <- export_source_H[1,a]
#       lat[b-4] <-export_source_H[2,a]
#       long[b-4] <-export_source_H[3,a]
      
#     } 
#   }
# }
# export_2[i,4] <- q[i]
# export_2[i,5] <-lat[i]
# export_2[i,6] <-long[i]
# export_2
# ss_month <-data.frame ("meansupplyd2H" = export_2[i,1], "meansupplylat" = export_2[i,2], "meansupplylong" = export_2[i,3],
#                        "meansourced2H" = export_2[i,4], "meansourcelat" = export_2[i,5], "meansourcelong" = export_2[i,6])

# ss_month2 <- transform(ss_month, id=match(meansourced2H, unique(meansourced2H)))
# all.equal(ss_month,ss_month2)
# ss_month2
# write.csv(ss_month2, file ="supply_source_may2015_H.csv")

# # ###################################
# # ## Now a 2-part mixture.

# # mixture.d18O <- matrix(nrow=5*5,ncol=21)

# # for (a in 1:nrow(source.reduced)) {
# # 	for (b in 1:nrow(source.reduced)) {
# # 		for (i in 1:21) {
# # 			mixture.d18O[(a-1)*5+b,i] <- ((i-1)/20)*source.reduced$mean_d18O[a] + (1-(i-1)/20)*source.reduced$mean_d18O[b] 
# # 		}
# # 	}
# # }


