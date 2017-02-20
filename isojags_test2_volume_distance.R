## Taken from v11 of original tests
## This uses a multivariate normal for the sources
## i: obs
## j: isotope
## k: source
##
## The likelihood then becomes:
## y[,j] ~ dmnorm( inprod( p*s_mean[1:K]), tau )
## tau should have Wishart prior
setwd("C:/Users/u0817369/Dropbox/Isotope Models/newscripts")
#setwd("C:/Users/yusu8/Dropbox/Isotope Models/newscripts")
library(rjags)
library(geosphere)
library(mcmcplots)
library(ggplot2)

jagsfile = "isojags_v11.txt"
rfile = paste0(substr(jagsfile,0,10),".Rdata")

#isotope = read.csv("C:/Users/u0817369/Dropbox/Isotope Models/isotope_data_may_2105.csv")
isotope = read.csv("C:/Users/yusu8/Dropbox/Isotope Models/isotope_data_may_2105.csv")


###############to extract the isotope data of the sources and the supply site
isite = subset(isotope, type=="supply_site")
isource = subset(isotope, type=="source_water")
#isource = isource[c(2,5),]


#############to calculate the distance between the different supplysites and the sources
isite.sp = isite
coordinates(isite.sp) <- ~lon+lat
proj4string(isite.sp) <- CRS("+init=epsg:4326")

isource.sp = isource
coordinates(isource.sp) <- ~lon+lat
proj4string(isource.sp) <- CRS("+init=epsg:4326")

distM = cbind( # Distance sites to source
  distGeo(p1=isite.sp, p2=isource.sp[1,]),
  distGeo(p1=isite.sp, p2=isource.sp[2,]),
  distGeo(p1=isite.sp, p2=isource.sp[3,]),
  distGeo(p1=isite.sp, p2=isource.sp[4,]),
  distGeo(p1=isite.sp, p2=isource.sp[5,])
)



# s_mean = cbind(isource$mean_d18O,isource$mean_d2H,isource$mean_DEX)
# y = cbind(isite$mean_d18O,isite$mean_d2H,isite$mean_DEX)
s_mean = cbind(isource$mean_d18O,isource$mean_d2H)
y = cbind(isite$mean_d18O,isite$mean_d2H)
# y = matrix(y[5,], ncol=2)
# y = s_mean + rnorm(dim(s_mean)[1] * dim(s_mean)[2])
# y = matrix(jitter(s_mean[2,], factor=0.01), ncol=3)

##  isotope Plots of the source and the suuply sites
plot.df = rbind(isite,isource)
ggplot(plot.df, aes(x=mean_d18O, y=mean_d2H, shape=type, col=type)) + geom_point(size=2)

# Dimensions of number of sources , supply sites and isotope ( we use 2 isotopes) matrix
Ntotal = dim(y)[1]
Jtotal = dim(y)[2]
Ktotal = dim(s_mean)[1]

# Alpha is base matrix for Dirichelet.We should have n alphas (alpha 1 to n, where n is the number of sources) 
#Here we have defined alpha as volume  and distance weighted for each source. Aditionally alpha is multiplied 
#by a large number to prevent the JAGS from getting stuck at a given value

alpha = matrix(1, nrow=Ntotal, ncol=Ktotal)
alpha = rep(1, Ktotal)
alpha = isource$volumes/sum(isource$volumes)*(distM^-1)*100000
alpha

## S3 is the identity matrix for the Wishart
S3 = diag(Jtotal)/4 

dataList = list(
  s_mean = s_mean,
  S3 = S3,
  y = y,
  alpha = alpha,
  K = Ktotal,
  N = Ntotal,
  J = dim(y)[2]
)

#using JAGS we are looking at the posterior distribution of p's and mu's. P's are the proportion value from different sources
# each p should be >0 and the sum of all p's = 1. mu is the value obtained using by multiplying the p'
parList = c("p", "mu" )

modelString = "
  model {
    # Likelihood
    # Now p must become p[i,1:k] row per obs, col for each source
    for (i in 1:N) {
      # Generate mu
      for (j in 1:J) {
        mu[i,j] = inprod(p[i,1:K],s_mean[1:K,j]) ### calculating mu as the product of source isotope value and its respective contribution
      }
      y[i,1:J] ~ dmnorm( mu[i,1:J], tau ) ##this is the likelihood function
    }
    # See cable eta l 2011 for details of the mathematics. The next lines are setting up the priors.
    # Dirichlet prior on p
    # Wrap this in N loop
    for (i in 1:N) {
      # First set proportions using Gamma prior
      for(k in 1:K) {
        f[i,k] ~ dgamma(alpha[i,k], 1)
      }
      # Now estimate p from sum of alpha
      for(k in 1:K) {
        p[i,k] <- f[i,k] / sum(f[i,1:K])
      }
    }
    # Tau - precision on sample variance
    # Wishart prior: J d.f. (can probably remove indices)
    tau[1:J,1:J] ~  dwish(S3[,],J)
    sigma[1:J,1:J] <- inverse(tau[,])
  }
"
cat( modelString, file="jags_isotopes_v11.txt")

jagsModel = jags.model( file = "jags_isotopes_v11.txt",
                        data = dataList,
                        n.chains = 3)

## Burnin
update(jagsModel, n.iter=15000)

## Run
codaSamples = coda.samples(jagsModel, variable.names=parList,
                           n.iter=300000)
gelman.diag(codaSamples, multivariate = FALSE)

x <- summary(codaSamples)
x$quantiles[c("p[1,1]","p[1,2]","p[1,3]","p[1,4]","p[1,5]"),3]
sum(x$quantiles[c("p[1,1]","p[1,2]","p[1,3]","p[1,4]","p[1,5]"),3])

save(isite, isource, jagsModel, codaSamples, file=rfile)

stop()
## test plots to analyze the results
traplot(codaSamples, c("p[1,1]","p[1,2]","p[1,3]","p[1,4]","p[1,5]"))
denplot(codaSamples, c("p[64,1]","p[64,2]","p[64,3]","p[64,4]","p[64,5]"))

isite.sp$pSRC1 <- x$quantiles[paste0("p[",seq(1,Ntotal),",1]"),][,3]
isite.sp$pSRC2 <- x$quantiles[paste0("p[",seq(1,Ntotal),",2]"),][,3]
isite.sp$pSRC3 <- x$quantiles[paste0("p[",seq(1,Ntotal),",3]"),][,3]
isite.sp$pSRC4 <- x$quantiles[paste0("p[",seq(1,Ntotal),",4]"),][,3]
isite.sp$pSRC5 <- x$quantiles[paste0("p[",seq(1,Ntotal),",5]"),][,3]
spplot(isite.sp, c("pSRC1"), 
       sp.layout=list(isource.sp, pch=c("1","2","3","4","5"), cex=1.2))
spplot(isite.sp, c("pSRC2"), 
       sp.layout=list(isource.sp, pch=c("1","2","3","4","5"), cex=1.2))
spplot(isite.sp, c("pSRC3"), 
       sp.layout=list(isource.sp, pch=c("1","2","3","4","5"), cex=1.2))
spplot(isite.sp, c("pSRC4"), 
       sp.layout=list(isource.sp, pch=c("1","2","3","4","5"), cex=1.2))
spplot(isite.sp, c("pSRC5"), 
       sp.layout=list(isource.sp, pch=c("1","2","3","4","5"), cex=1.2))

isite$pSRC1 <- x$quantiles[paste0("p[",seq(1,Ntotal),",1]"),][,3]
ggplot(isite, aes(x=mean_d18O, y=mean_d2H, col=pSRC1)) + geom_point(size=2)
isite$pSRC2 <- x$quantiles[paste0("p[",seq(1,Ntotal),",2]"),][,3]
ggplot(isite, aes(x=mean_d18O, y=mean_d2H, col=pSRC2)) + geom_point(size=2)
isite$pSRC3 <- x$quantiles[paste0("p[",seq(1,Ntotal),",3]"),][,3]
ggplot(isite, aes(x=mean_d18O, y=mean_d2H, col=pSRC3)) + geom_point(size=2)
isite$pSRC4 <- x$quantiles[paste0("p[",seq(1,Ntotal),",4]"),][,3]
ggplot(isite, aes(x=mean_d18O, y=mean_d2H, col=pSRC4)) + geom_point(size=2)
isite$pSRC5 <- x$quantiles[paste0("p[",seq(1,Ntotal),",5]"),][,3]
ggplot(isite, aes(x=mean_d18O, y=mean_d2H, col=pSRC5)) + geom_point(size=2)

## Comparison: d180
pred.d180 = x$statistics[paste0("mu[",seq(1,Ntotal),",1]"),"Mean"]
obs.d180 = isite$mean_d18O
## RMSE
sqrt(mean((pred.d180 - obs.d180)^2))
max(abs((pred.d180 - obs.d180)))

## Comparison: d2H
pred.d2H = x$statistics[paste0("mu[",seq(1,Ntotal),",2]"),"Mean"]
obs.d2H = isite$mean_d2H
## RMSE
sqrt(mean((pred.d2H - obs.d2H)^2))
max(abs((pred.d2H - obs.d2H)))
