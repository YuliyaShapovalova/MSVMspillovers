library("xtable")
library("coda")
library("stats")
#library('latex2exp')
Mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
setpar<-function(...) par(mar=c(5.1,4.1,4.1,2.1),mgp=c(1.4,0.2,0),
                          cex=0.7,bty="l",las=1,lwd=0.5,tck=0.01,...)

setwd("/Users/Yuliya/Dropbox/2dimSVspillovers_newmac/figures")

results_name = paste("results_00_","seed1235_","15042021",".RData", sep="")
load(paste("/Users/Yuliya/Dropbox/2dimSVspillovers_newmac/",results_name, sep = ""))
save_name = paste("_seed1235_","15042021", sep="")

BurnIn=5000
M=10000
thetamh=results$thetamh
HPD<-HPDinterval(mcmc(thetamh[BurnIn:M,,1]))

tb<-cbind(apply(thetamh[BurnIn:M,,1],2,mean),apply(thetamh[BurnIn:M,,1],2,median),
          apply(thetamh[BurnIn:M,,1],2,Mode),HPD[,1],HPD[,2])

xtable(tb,digits=4)
source("/Users/Yuliya/Dropbox/2dimSVspillovers_newmac/trivariate_apf.R")

thetaname=c("phi", "mu[2]", "phi[11]","phi[21]",
            "phi[12]","phi[22]", "sigma[11]","sigma[22]",
            "rho")

thetaname_mu = c("mu[1]")

n=2
mu_m=c(0.025,0.025)
phi_m=matrix(c(0.0),n,n)
diag(phi_m)=c(0.8,0.8)
phi_m[1,2]=0.0
phi_m[2,1]=0.0
seta_m<-diag(c(0.1,0.1),n)
Rhotrue=0.6
seps_m<-matrix(Rhotrue,n,n)
diag(seps_m)=1.0
thetatrue=c(mu_m,phi_m,diag(seta_m), Rhotrue)

##################
# Plot mu
##################
setEPS()
postscript(width=10, height = 5,paste("results_mu", save_name, ".eps", sep=""))
setpar(mfrow=c(2,3),mar=c(3.5,3.2,3.0,2.0),  cex.lab = 1.2)
high=c(max(results$thetamh[,1,1]),max(results$thetamh[,2,1]))
low=c(min(results$thetamh[,1,1]),min(results$thetamh[,2,1]))

plot(results$thetamh[,1,1],type="l",ylim=c(low[1],high[1]),main=TeX('Trace plot for $\\mu_1$'), xlab='', ylab=parse(text=paste(thetaname[1])))
abline(h=thetatrue[1], col="darkred", lwd = 3)
hist(results$thetamh[BurnIn:M,1,1], xlab='', ylab="", breaks=20, main=TeX('Histogram for $\\mu_1$'))
abline(v =thetatrue[1], col = "darkred", lwd = 3)
acf(results$thetamh[BurnIn:M,1,1], lag=1000, main=TeX('ACF for $\\mu_1$'))

plot(results$thetamh[,2,1],type="l",ylim=c(low[2],high[2]),main=TeX('Trace plot for $\\mu_2$'),xlab='',ylab=parse(text=paste(thetaname[2])))
abline(h=thetatrue[2],col="darkred", lwd = 3)
hist(results$thetamh[BurnIn:M,2,1], xlab='', ylab="", breaks=20, main=TeX('Histogram for $\\mu_2$'))
abline(v =thetatrue[2], col = "darkred", lwd = 3)
acf(results$thetamh[BurnIn:M,2,1], lag=1000, main=TeX('ACF for $\\mu_2$'))
dev.off()


##################
# Plot phi
##################
setEPS()
postscript(width=10, height = 10,paste("results_phi", save_name, ".eps", sep=""))
setpar(mfrow=c(4,3),mar=c(2.0,3.2,3.0,2.0), cex.lab = 1.2,mgp=c(2.0,0.2,0))
high=c(max(results$thetamh[,3,1]),max(results$thetamh[,4,1]),max(results$thetamh[,5,1]),
       max(results$thetamh[,6,1]))
low=c(min(results$thetamh[,3,1]),min(results$thetamh[,4,1]),min(results$thetamh[,5,1]),
      min(results$thetamh[,6,1]))

plot(results$thetamh[,3,1],type="l",ylim=c(low[1],high[1]),ylab=parse(text=paste(thetaname[3])), xlab='',main=TeX('Trace plot for $\\phi_{11}$'))
abline(h=thetatrue[3],col = "darkred", lwd = 3)
hist(results$thetamh[BurnIn:M,3,1], xlab=parse(text=paste(thetaname[3])), ylab="", breaks=20, main=TeX('Histogram for $\\phi_{11}$'))
abline(v =thetatrue[3], col = "darkred", lwd = 3)
acf(results$thetamh[BurnIn:M,3,1], lag=1000, main=TeX('ACF for $\\phi_{11}$'))

plot(results$thetamh[,4,1],type="l",ylim=c(low[2],high[2]),ylab=parse(text=paste(thetaname[4])),xlab='', main=TeX('Trace plot for $\\phi_{12}$'))
abline(h=thetatrue[4],col = "darkred", lwd = 3)
hist(results$thetamh[BurnIn:M,4,1], xlab=parse(text=paste(thetaname[4])), ylab="", breaks=20, main=TeX('Histogram for $\\phi_{12}$'))
abline(v =thetatrue[4], col = "darkred", lwd = 3)
acf(results$thetamh[BurnIn:M,4,1],  lag=1000, main=TeX('ACF for $\\phi_{12}$'))

plot(results$thetamh[,5,1],type="l",ylim=c(low[3],high[3]),ylab=parse(text=paste(thetaname[5])),xlab='', main=TeX('Trace plot for $\\phi_{21}$'))
abline(h=thetatrue[5],col = "darkred", lwd = 3)
hist(results$thetamh[BurnIn:M,5,1], xlab=parse(text=paste(thetaname[5])), ylab="", breaks=20, main=TeX('Histogram for $\\phi_{21}$'))
abline(v =thetatrue[5], col = "darkred", lwd = 3)
acf(results$thetamh[BurnIn:M,5,1], lag=1000, main=TeX('ACF for $\\phi_{21}$'))

plot(results$thetamh[,6,1],type="l",ylim=c(low[4],high[4]),ylab=parse(text=paste(thetaname[6])),xlab='', main=TeX('Trace plot for $\\phi_{22}$'))
abline(h=thetatrue[6],col = "darkred", lwd = 3)
hist(results$thetamh[BurnIn:M,6,1], xlab=parse(text=paste(thetaname[6])), ylab="", breaks=20, main=TeX('Histogram for $\\phi_{22}$'))
abline(v =thetatrue[6], col = "darkred", lwd = 3)
acf(results$thetamh[BurnIn:M,6,1], lag=1000, main=TeX('ACF for $\\phi_{22}$'))
dev.off()


######################
# Plot sigmas and rho
######################
setEPS()
postscript(width=10, height = 7.5, paste("results_sigma", save_name, ".eps", sep=""))
setpar(mfrow=c(3,3),mar=c(3.0,3.5,3.0,2.0),cex.lab = 1.2,mgp=c(2.0,0.2,0))
high=c(max(results$thetamh[,7,1]),max(results$thetamh[,8,1]),max(results$thetamh[,9,1]))
low=c(min(results$thetamh[,7,1]),min(results$thetamh[,8,1]),min(results$thetamh[,9,1]))

plot(results$thetamh[,7,1],type="l",ylim=c(low[1],high[1]),ylab=parse(text=paste(thetaname[7])), xlab='',main=TeX('Trace plot for $\\sigma_{1}$'))
abline(h=thetatrue[7],col = "darkred", lwd = 3)
hist(results$thetamh[BurnIn:M,7,1], xlab='', ylab="", breaks=20, main=TeX('Histogram for $\\sigma_{1}$'))
abline(v =thetatrue[7], col = "darkred", lwd = 3)
acf(results$thetamh[BurnIn:M,7,1],main=TeX('ACF for $\\sigma_{1}$'), lag=1000)

plot(results$thetamh[,8,1],type="l",ylim=c(low[2],high[2]),ylab=parse(text=paste(thetaname[8])),xlab='', main=TeX('Trace plot for $\\sigma_{2}$'))
abline(h=thetatrue[8],col = "darkred", lwd = 3)
hist(results$thetamh[BurnIn:M,8,1], xlab='', ylab="", breaks=20, main=TeX('Histogram for $\\sigma_{2}$'))
abline(v =thetatrue[8], col = "darkred", lwd = 3)
acf(results$thetamh[BurnIn:M,8,1],main=TeX('ACF for $\\sigma_{2}$'), lag=1000)

plot(results$thetamh[,9,1],type="l",ylim=c(low[3],high[3]),ylab=parse(text=paste(thetaname[9])),xlab='', main=TeX('Trace plot for $\\rho$'))
abline(h=thetatrue[9],col = "darkred", lwd = 3)
hist(results$thetamh[BurnIn:M,9,1], xlab='', ylab="", breaks=20, main=TeX('Histogram for $\\rho$'))
abline(v =thetatrue[9], col = "darkred", lwd = 3)
acf(results$thetamh[BurnIn:M,9,1],main=TeX('ACF for $\\rho'),lag=1000)

dev.off()



