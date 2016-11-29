# slv_tapwater_source_model.R


# load data file into dataframe

# Yusuf's computer
setwd("C:/Users/yusu8/Dropbox/Utah_water_isotopes/JVWCD")
#setwd("C:/Users/u0817369/Dropbox/Utah_water_isotopes/JVWCD")
#setwd("C:/Users/yusuf/Dropbox/Utah_water_isotopes/JVWCD")
# Rich's computer
#setwd("~/Downloads")
data <- read.csv("May_5.csv",header=TRUE)

# break into source and supplyline data frames
source.data <- subset(data,data$Site_Description=="Source")
supply.data <- subset(data,data$Site_Description=="Supplyline")

# average multiple measurements of individual sites.
#----------------------------------------------------

# processing source dataframe

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
	#######SITE ID
	source.uniqueID[i] <- subset(source.data$Site_ID,source.data$Site_ID==source.siteIDs)
}

# combine vectors back into single, reduced dataframe - STDEV DEX missing!
source.reduced <- data.frame("mean_d18O"=source.mean_d18O,"mean_d2H"=source.mean_d2H,"mean_DEX"=source.mean_DEX,
	"stdev_d18O"=source.stdev_d18O,"stdev_d2H"=source.stdev_d2H,"stdev_DEX"=source.stdev_DEX, "month"=source.month,"year"=source.year,
	"lat"=source.lat,"lon"=source.lon, "SITE_ID" = source.siteIDs,"uniqueIDS" = source.uniqueID)

#------------------------------
# processing supply dataframe

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

#==========================================================================

##############################################
# calculation of likelihood values for d2H and DEX combined
##############################################

lvals_H <- matrix(nrow=nrow(supply.reduced),ncol=nrow(source.reduced))
n<-2
s1 <- 0.1
s2 <- 0.3
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2)
library(MASS)
library(mvtnorm)
upper <- matrix(nrow=nrow(supply.reduced),ncol=2)
lower <- matrix(nrow=nrow(supply.reduced),ncol=2)
mean <- matrix(nrow=nrow(source.reduced),ncol=2)
for (i in 1:nrow(supply.reduced)) {
  for (j in 1:nrow(source.reduced)) {
    # assign likelihood for given i,j
    upper[i,] = c((supply.reduced$mean_d2H[i] + supply.reduced$stdev_d2H[i]), (supply.reduced$mean_DEX[i] + supply.reduced$stdev_DEX[i]))
    lower[i,] = c((supply.reduced$mean_d2H[i] - supply.reduced$stdev_d2H[i]), (supply.reduced$mean_DEX[i] - supply.reduced$stdev_DEX[i]))
    mean[j,] = c((source.reduced$mean_d2H[j]),(source.reduced$mean_DEX[j]))
    
    lvals_H[i,j] <-pmvnorm(upper = upper[i,], lower = lower[i,], mean = mean[,j], sigma = sigma)
  }
}

colnames(lvals_H) <- c("solena_way", "JVWTP", "SEWTP", "SWWTP" , "WELL_1300",  "M1", "M2")
lH<-data.frame(lvals_H)

###########################
# dataframe for likelihood functions
################
lvals_dH <- data.frame("SITE_ID" = supply.reduced$SITE_ID,
                       "solena_way" = lH$solena_way, "JVWTP" = lH$JVWTP, "SEWTP" =lH$SEWTP, 
                       "SWWTP" =lH$SWWTP,"WELL_1300" = lH$WELL_1300, "M1" = lH$M1, "M2" = lH$M2  )
ll_d2H <- plot(x= c(1,2,3,4,5,6,7),lvals_dH[1,2:8], col ="red", type = "l", xlim = c(0.8,7.1), ylim = c(0,1) )
title(main = "d2H_LIKELIHOOD")
p = sample(rainbow(63))
for (i in 2:64){
  lines(x= c(1,2,3,4,5,6,7),lvals_dH[i,2:8], type = "l", col = p[i])
}

#################
#calculation of posterior distributions
post_H <- matrix(nrow=nrow(lvals_H), ncol=ncol(lvals_H))
i <- 1:nrow(lvals_H)
j <- 1:ncol(lvals_H)
# assign posterior for given i,jsupply.
post_H[i,j] <- (lvals_H[i,j]*ncol(lvals_H))/(ncol(lvals_H)* rowSums(lvals_H))
#rowSums(post)
#colSums(post)
#########################################
##plotting the POSTERIOR values
pp_H <- plot(x= c(1,2,3,4,5,6,7),post_H[1,1:7], col ="red", type = "l", xlim = c(0.8,7.1), ylim = c(0,1) )
title(main = "d2H_POSTERIOR")
p = sample(rainbow(63))
for (i in 2:64){
  lines(x= c(1,2,3,4,5,6,7),post_H[i,1:7], type = "l", col = p[i])
}

pp_export <- data.frame(post_H, supply.reduced$lat, supply.reduced$lon, supply.reduced$SITE_ID)
#write.csv(pp_export, file ="post_may2015_H.csv")

require("rgdal") # requires sp, will use proj.4 if installed
#install.packages("gpclib")
library("maptools")
require(ggplot2)
gpclibPermit()
require("plyr")
require(rgeos)


#shapefile <- readOGR(dsn = "C:/Users/u0817369/Dropbox/Utah_water_isotopes/slc_county_shape_file _1" ,layer = "counties_clip")
shapefile <- readOGR(dsn = "C:/Users/u0817369/Dropbox/Utah_water_isotopes/slc_county_shape_file _1" ,layer = "Export_Output")
proj4string(shapefileDF) <- CRS("+proj=longlat")
#shapefile <- readShapePoly( "Export_Output.shp")
shapefile@data$id <- rownames(shapefile@data)
shapefile1 <- fortify(shapefile, region = "id")
shapefileDF <- merge(shapefile1, shapefile@data, by = "id")

y = ggplot(shapefileDF, aes(long,lat, group = group, fill = "WR_ID")) + geom_polygon() + geom_path(color = "black") + coord_equal() + scale_fill_manual(values=c("white"), guide="none")# + coord_map("albers", lat0=40, lat1=40) #+ (color = 'gray', fill = 'blue', size = 2) 
y

coordinates(pp_export)<-~supply.reduced.lon+supply.reduced.lat
proj4string(pp_export)
proj4string(pp_export)<-CRS("+proj=longlat +datum=NAD83")
pp_export1<-spTransform(pp_export, CRS(proj4string(shapefile)))
identical(proj4string(pp_export1),proj4string(shapefile))
pp_export1<-data.frame(pp_export1)

source.reduced1<-subset(source.reduced[1:5,])
coordinates(source.reduced1)<-~lon+lat
proj4string(source.reduced1)
proj4string(source.reduced1)<-CRS("+proj=longlat +datum=NAD83")
source.reduced1<-spTransform(source.reduced1, CRS(proj4string(shapefile)))
#identical(proj4string(pp_export1),proj4string(shapefile))
source.reduced1<-data.frame(source.reduced1)




require(ggplot2)
#for (i in 1:4) {
#i <- 1:4
x =  y +  geom_point(data = pp_export1, aes(x=supply.reduced.lon, y=supply.reduced.lat, group = X7, color = X7), size = 4) + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x = x + geom_point(aes(x=source.reduced1$lon[5], y= source.reduced1$lat[5]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
x = x + geom_point(aes(x=source.reduced1$lon[4], y= source.reduced1$lat[4]), size =5, shape = 23, color ="black", fill = "black") + ggtitle("d2H")
x
pdf("source7_H_may.ai")
print(x)
dev.off()
#}

rel_likelihood_H <- matrix(nrow=nrow(post),ncol=ncol(post))

for (i in 1:nrow(post_H)) {
  rel_likelihood_H[i,] <- post_H[i,]/max(post_H[i,])
  
}


##############################################
# calculation of likelihood values for DEX
##############################################

lvals_DEX <- matrix(nrow=nrow(supply.reduced),ncol=nrow(source.reduced))
for (i in 1:nrow(supply.reduced)) {
  for (j in 1:nrow(source.reduced)) {
    # assign likelihood for given i,j
    lvals_DEX[i,j] <- pnorm((supply.reduced$mean_DEX[i] + supply.reduced$stdev_DEX[i]), 
                        mean = source.reduced$mean_DEX[j], sd = source.reduced$stdev_DEX[j], lower.tail = TRUE) - 
      pnorm((supply.reduced$mean_DEX[i] - supply.reduced$stdev_DEX[i]), 
            mean = source.reduced$mean_DEX[j], sd = source.reduced$stdev_DEX[j], lower.tail = TRUE)
  }
}
colnames(lvals_DEX) <- c("solena_way", "JVWTP", "SEWTP", "SWWTP" , "WELL_1300" ,  "M1", "M2" )
l_DEX<-data.frame(lvals_DEX)

###########################
# PLOTTING THE LIKELIHOOD FOR dex
################
lvals_DEX_df <- data.frame("SITE_ID" = supply.reduced$SITE_ID, "solena_way" = l_DEX$solena_way, "JVWTP" = l_DEX$JVWTP,
                           "SEWTP" =l_DEX$SEWTP, "SWWTP" =l_DEX$SWWTP,"WELL_1300" = l_DEX$WELL_1300, "M1" = l_DEX$M1, "M2" = l_DEX$M2  )
ll_DEX <- plot(x= c(1,2,3,4,5,6,7),lvals_DEX_df[1,2:8], col ="red", type = "b", xlim = c(0.8,7.1), ylim = c(0,1) )
title(main = "D-EXCESS_LIKELIHOOD")
p = sample(rainbow(63))
for (i in 2:64){
  lines(x= c(1,2,3,4,5,6,7),lvals_DEX_df[i,2:8], type = "b", col = p[i])
}

#calculation of posterior distributions
post_DEX <- matrix(nrow=nrow(lvals_DEX), ncol=ncol(lvals_DEX))
i <- 1:nrow(lvals_DEX)
j <- 1:ncol(lvals_DEX)
# assign posterior for given i,jsupply.
post_DEX[i,j] <- (lvals_DEX[i,j]*ncol(lvals_DEX))/(ncol(lvals_DEX)* rowSums(lvals_DEX))
#rowSums(post)
#colSums(post)
#########################################
##plotting the POSTERIOR values
pp_DEX <- plot(x= c(1,2,3,4,5,6,7),post_DEX[1,1:7], col ="red", type = "l", xlim = c(0.8,7.1), ylim = c(0,1) )
title(main = "DEX_POSTERIOR")
p = sample(rainbow(63))
for (i in 2:64){
  lines(x= c(1,2,3,4,5,6,7),post_DEX[i,1:7], type = "l", col = p[i])
}

pp_DEX_exp <- data.frame(post_DEX, supply.reduced$lat, supply.reduced$lon, supply.reduced$SITE_ID)
#write.csv(pp_export, file ="post_may2015_H.csv")


require(ggplot2)
x = ggplot(pp_DEX_exp, aes(x=supply.reduced.lon, y=supply.reduced.lat, color=X5)) + geom_point(size=8) 
x = x + ggtitle("DEX") + scale_colour_gradient(limits=c(0, 1), low="grey", high="red")
x





###############################################################
###############################################################
##############3xporting data with supply and corresponding sources d18O
####################################################################
export_source_O <- matrix(nrow=nrow(post)+4, ncol=ncol(post)+4)
i<-5:nrow(export_source_O)
a<- 1:nrow(rel_likelihood)
b<- 1:ncol(rel_likelihood)
export_source_O[i,1] <- supply.reduced[,10]
export_source_O[i,2] <- supply.reduced[,1]
export_source_O[i,3] <- supply.reduced[,8]
export_source_O[i,4] <- supply.reduced[,9]
j<-5:ncol(export_source_O)
export_source_O[i,j] <-rel_likelihood[a,b]
export_source_O[1,j] <- source.reduced[,1]
export_source_O[2,j] <- source.reduced[,8]
export_source_O[3,j] <- source.reduced[,9]
export_source_O[4,j] <-j-4
export_source_O
############################
##########################################
###############source and supply matrix combined
export_2 <- matrix(nrow=64,ncol=6)
i <- 1:nrow(export_2)
export_2[i,1] <- export_source_O[i+4,2]
export_2[i,2] <- export_source_O[i+4,3]
export_2[i,3] <- export_source_O[i+4,4]
q <-vector()
lat<-vector()
long<-vector()
for (b in 5:68) {
  for (a in 5:23) {
    if (export_source_O[b,a] == 1) {
     q[b-4] <- export_source_O[1,a]
     lat[b-4] <-export_source_O[2,a]
     long[b-4] <-export_source_O[3,a]
     
    } 
 }
}
export_2[i,4] <- q[i]
export_2[i,5] <-lat[i]
export_2[i,6] <-long[i]
export_2
ss_month <-data.frame ("meansupplyd18O" = export_2[i,1], "meansupplylat" = export_2[i,2], "meansupplylong" = export_2[i,3],
                       "meansourced18O" = export_2[i,4], "meansourcelat" = export_2[i,5], "meansourcelong" = export_2[i,6])

ss_month2 <- transform(ss_month, id=match(meansourced18O, unique(meansourced18O)))
all.equal(ss_month,ss_month2)
ss_month2
write.csv(ss_month2, file ="supply_source_may2015_O.csv")
####################################
###############################################################
##############3xporting data with supply and corresponding sources d2H
####################################################################
export_source_H <- matrix(nrow=nrow(post)+4, ncol=ncol(post)+4)
i<-5:nrow(export_source_H)
a<- 1:nrow(rel_likelihood_H)
b<- 1:ncol(rel_likelihood_H)
export_source_H[i,1] <- supply.reduced[,10]
export_source_H[i,2] <- supply.reduced[,2]
export_source_H[i,3] <- supply.reduced[,8]
export_source_H[i,4] <- supply.reduced[,9]
j<-5:ncol(export_source_H)
export_source_H[i,j] <-rel_likelihood_H[a,b]
export_source_H[1,j] <- source.reduced[,2]
export_source_H[2,j] <- source.reduced[,8]
export_source_H[3,j] <- source.reduced[,9]
export_source_H[4,j] <-j-4
export_source_H
############################
##########################################
###############source and supply matrix combined
export_2 <- matrix(nrow=64,ncol=6)
i <- 1:nrow(export_2)
export_2[i,1] <- export_source_H[i+4,2]
export_2[i,2] <- export_source_H[i+4,3]
export_2[i,3] <- export_source_H[i+4,4]
q <-vector()
lat<-vector()
long<-vector()
for (b in 5:68) {
  for (a in 5:23) {
    if (export_source_H[b,a] == 1) {
      q[b-4] <- export_source_H[1,a]
      lat[b-4] <-export_source_H[2,a]
      long[b-4] <-export_source_H[3,a]
      
    } 
  }
}
export_2[i,4] <- q[i]
export_2[i,5] <-lat[i]
export_2[i,6] <-long[i]
export_2
ss_month <-data.frame ("meansupplyd2H" = export_2[i,1], "meansupplylat" = export_2[i,2], "meansupplylong" = export_2[i,3],
                       "meansourced2H" = export_2[i,4], "meansourcelat" = export_2[i,5], "meansourcelong" = export_2[i,6])

ss_month2 <- transform(ss_month, id=match(meansourced2H, unique(meansourced2H)))
all.equal(ss_month,ss_month2)
ss_month2
write.csv(ss_month2, file ="supply_source_may2015_H.csv")

###################################
## Now a 2-part mixture.

mixture.d18O <- matrix(nrow=5*5,ncol=21)

for (a in 1:nrow(source.reduced)) {
	for (b in 1:nrow(source.reduced)) {
		for (i in 1:21) {
			mixture.d18O[(a-1)*5+b,i] <- ((i-1)/20)*source.reduced$mean_d18O[a] + (1-(i-1)/20)*source.reduced$mean_d18O[b] 
		}
	}
}


