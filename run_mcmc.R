setwd("/Users/Yuliya/MSVMspillovers")
source("utils/utils.R")
set.seed(1235)

n=2      # dimension of the models
T=1000   # lingth of time series   
N=10000  # number of particles
M=10000   # number of MH steps

BurnIn=1000 # only relevant for plotting/posterior inference

# parameters of the model for data generation
mutrue<-c(0.025,0.025)
phitrue<-matrix(c(0.0),n,n)
diag(Phitrue)=c(0.8,0.8)
phitrue[1,2]=0.4
phitrue[2,1]=0.0
eigen(phitrue)$values
setatrue<-diag(c(0.1,0.1),n)
rhotrue=0.6
sepstrue<-matrix(rhotrue,n,n)
diag(sepstrue)=1.0
######################################################################
# Generate data 
######################################################################
Data=SimulateStochVolModel(T,mutrue,phitrue,setatrue,sepstrue)
Y=Data$Y
n = dim(Y)[1]
T = dim(Y)[2]

######################################################################
# Proposal steps
######################################################################
constcold=0.5
vucold=constcold*c(rep(0.1,n),0.01,rep(0.09,n),0.01,rep(0.05,n),0.03) 
consthot=2.0
vuhot=consthot*c(rep(0.1,n),0.01,rep(0.09,n),0.01,rep(0.05,n),0.03) 

vuphi1hot=consthot*c(0.001,0.01)*1.8
vuphi1cold=constcold*c(0.001,0.01)*2.8
vuphi2hot=consthot*c(0.01,0.001)*1.5
vuphi2cold=constcold*c(0.01,0.001)*1.5  

#Mu
numhmucold=rmultnorm(M,rep(0,n),diag((vucold[1:n])),Cholesky=TRUE)
numhmuhot=rmultnorm(M,rep(0,n),diag((vuhot[1:n])),Cholesky=TRUE)
#Phi
vuphi1hot=consthot*c(0.01,0.1)*0.8
vuphi1cold=constcold*c(0.01,0.1)*0.8
vuphi2hot=consthot*c(0.1,0.01)*1.2
vuphi2cold=constcold*c(0.1,0.01)*1.2

vuphihot=consthot*c(0.001,0.1,0.1,0.001)
vuphicold=constcold*c(0.001,0.1,0.1,0.001)  

numhphicold=rmultnorm(M,rep(0,n*n),diag(vuphicold),Cholesky=TRUE)
numhphihot=rmultnorm(M,rep(0,n*n),diag(vuphihot),Cholesky=TRUE)

numhphi1cold=rmultnorm(M,rep(0,n),diag(vuphi1cold),Cholesky=TRUE)
numhphi1hot=rmultnorm(M,rep(0,n),diag(vuphi1hot),Cholesky=TRUE)
numhphi2cold=rmultnorm(M,rep(0,n),diag(vuphi2cold),Cholesky=TRUE)
numhphi2hot=rmultnorm(M,rep(0,n),diag(vuphi2hot),Cholesky=TRUE)
#Seta
numhsighot=consthot*rmultnorm(M,rep(0,n),diag(vuhot[(n^2+n+1):(n^2+2*n)]),Cholesky=TRUE)
numhsigcold=constcold*rmultnorm(M,rep(0,n),diag(vucold[(n^2+n+1):(n^2+2*n)]),Cholesky=TRUE)
#Rho
numhrhocold=rnorm(M,0,(vucold[n^2+2*n+1]))
numhrhohot=rnorm(M,0,(vuhot[n^2+2*n+1]))

######################################################################
# Initial values for the sampler 
######################################################################
mustart=matrix(rep(0.01,n),n,1)
phistart=diag(c(0.8,0.8))
setastart=diag(c(0.01,0.01))
rhostart=0.1
sepsstart=matrix(rep(rhostart,n^2), n, n)
diag(sepsstart)=1.0
thetamh=array(0,c(M,(n+n*n+n+1),1))
thetamh[1,,1]=c(mustart,c(phistart),diag(setastart),rhostart)

Likmh=rep(0,M)
LikFin=likfort2020(n,T,N,Y,mustart,phistart,setastart,sepsstart)
proc.time() - ptm
Likmh[1]=LikFin
LikC=rep(0,M)
#Initial value of Lik is taken from particle filter algorithm    
q=matrix(0,8,M)
thetac=matrix(0,1,(n^2+2*n+1+n*n))
Likratio=rep(0,M)
Priorratio=rep(0,M)
Priorratiomh=rep(0,M)
PriorratioC=rep(0,M)

######################################################################
# Indicators for proposal parameters
######################################################################
Imu=c(1:n)
Iphi1=(n+1):(2*n)
Iphi2=(2*n+1):(3*n)
Iseta=(n*n+n+1):(n*n+2*n)
Irho=(n*n+2*n+1)

U=matrix(runif(8*M),8,M)
S=matrix(runif(8*M),8,M)
s=0.33

for (j in 2:M){
  thetamh[j,,1]=thetamh[j-1,,1]  
  #######################
  # proposal for means 
  #######################
  thetac=thetamh[j,,1]
  Likmh[j]=Likmh[j-1]
  if(S[1,j]<s){
    thetac[Imu]=thetamh[j,Imu,1]+t(numhmuhot[,j])
  }else{
    thetac[Imu]=thetamh[j,Imu,1]+t(numhmucold[,j])    
  }
  Mumh=thetac[Imu]
  Phimh=matrix(thetac[Iphi],n,n)
  SigmaEtamh=diag(thetac[Iseta])    
  SigmaEpsmh=matrix(thetac[Irho],n,n)
  diag(SigmaEpsmh)=1.0
  #Calculate likelihood with thetac as a parameter vector...  
  #record likelihood of a candidate on step j     
  LikC[j]=likfort2020(n,T,N,Y,Mumh,Phimh,SigmaEtamh,SigmaEpsmh)
  #LikC[j]=LogLikl(T,n,N,Mumh,Phimh,SigmaEtamh,SigmaEpsmh,Y)
  if(is.nan(LikC[j])){
    LikC[j]=-Inf
  }
  #calculate the acceptance probablity
  Likratio[j] = LikC[j]-Likmh[j-1]   
  Priorratio[j] = logPriorValue(thetac) - logPriorValue(thetamh[j-1,,1])
  q[1,j]=min(1,exp(Likratio[j]+Priorratio[j]))    
  if(is.nan(q[1,j])){
    q[1,j]=-Inf
  }
  #q[1,j]=min(0,(Likratio[j]+Priorratio[j]))      
  # Check https://umbertopicchini.wordpress.com/2017/12/18/tips-for-coding-a-metropolis-hastings-sampler/
  
  #accept or reject candidate
  u=runif(1,0,1)
  if(u<q[1,j]){
    thetamh[j,,1]=thetac
    Likmh[j]=LikC[j]
  }
  
  #############################
  #proposal for phi, First row
  #############################
  thetac=thetamh[j,,1]
  if(S[2,j]<s){
    thetac[Iphi1]=thetamh[j,Iphi1,1]+t(numhphi1hot[,j])
  }else{
    thetac[Iphi1]=thetamh[j,Iphi1,1]+t(numhphi1cold[,j])
  }
  
  Mumh=matrix(thetac[Imu],n,1)
  Phimh=matrix(c(thetac[Iphi1],thetamh[j,Iphi2,1]),n,n)
  SigmaEtamh=diag(thetac[Iseta])
  SigmaEpsmh=matrix(thetac[Irho],n,n)
  diag(SigmaEpsmh)=1.0
  Phiev=eigen(Phimh)$values
  if (any(abs(Phiev)>0.98)) {
    #q[2,j]=0
    q[2,j]=-10000000.0
    thetamh[j,,1]=thetamh[j,,1]
    Likmh[j]=Likmh[j]
  } else {
    #Calculate likelihood with thetac as a parameter vector...
    #record likelihood of a candidate on step j
    LikC[j]=likfort2020(n,T,N,Y,Mumh,Phimh,SigmaEtamh,SigmaEpsmh)
    #calculate the acceptance probablity
    Likratio[j]= LikC[j]-Likmh[j]
    Priorratio[j]=logPriorValue(thetac)-logPriorValue(thetamh[j,,1])
    q[2,j]= min(1,exp(Likratio[j]+Priorratio[j]))
    
  }
  if(is.nan(q[2,j])){
    q[2,j]=-Inf
  }
  #accept or reject candidate
  u=runif(1,0,1)
  if (u<q[2,j]) {
    thetamh[j,,1]=thetac
    Likmh[j]=LikC[j]
  }
  
  # #############################
  # #proposal for phi, Second row
  # #############################
  thetac=thetamh[j,,1]
  if(S[3,j]<s){
    thetac[Iphi2]=thetamh[j,Iphi2,1]+t(numhphi2hot[,j])
  }else{
    thetac[Iphi2]=thetamh[j,Iphi2,1]+t(numhphi2cold[,j])
  }
  
  Mumh=matrix(thetac[Imu],n,1)
  Phimh=matrix(c(thetamh[j,Iphi1,1],thetac[Iphi2]),n,n)
  SigmaEtamh=diag(thetac[Iseta])
  SigmaEpsmh=matrix(thetac[Irho],n,n)
  diag(SigmaEpsmh)=1.0
  Phiev=eigen(Phimh)$values
  if (any(abs(Phiev)>0.98)) {
    #q[3,j]=0
    q[3,j]=-10000000.0
    thetamh[j,,1]=thetamh[j,,1]
    Likmh[j]=Likmh[j]
  } else {
    #Calculate likelihood with thetac as a parameter vector...
    #record likelihood of a candidate on step j
    LikC[j]=likfort2020(n,T,N,Y,Mumh,Phimh,SigmaEtamh,SigmaEpsmh)
    
    #calculate the acceptance probablity
    Likratio[j]=LikC[j]-Likmh[j]
    Priorratio[j]=logPriorValue(thetac)-logPriorValue(thetamh[j,,1])
    q[3,j]=min(1,exp(Likratio[j] + Priorratio[j]))
    if(is.nan(q[3,j])){
      q[3,j]=-Inf
    }
  }
  #accept or reject candidate
  u=runif(1,0,1)
  if (u<q[3,j]) {
    thetamh[j,,1]=thetac
    Likmh[j]=LikC[j]
  }
  
  #######################
  # proposal for sigmas 
  #######################
  thetac=thetamh[j,,1]
  
  if(S[4,j]<s){
    thetac[Iseta]=thetamh[j,Iseta,1]+numhsighot[,j]
  }else{
    thetac[Iseta]=thetamh[j,Iseta,1]+numhsigcold[,j]
  }
  
  Mumh=matrix(thetac[Imu],n,1)
  Phimh=matrix(thetac[Iphi],n,n)
  SigmaEtamh=diag(thetac[Iseta])    
  SigmaEpsmh=matrix(thetac[Irho],n,n)
  diag(SigmaEpsmh)=1.0
  #posdef=ifelse(isPositiveDefinite(SigmaEtamh),1,0)
  if (any(thetac[Iseta]<0.00001)) {
    q[4,j]=-1000000.0
  }else{
    #Calculate likelihood with thetac as a parameter vector...       
    #record likelihood of a candidate on step j     
    LikC[j]=likfort2020(n,T,N,Y,Mumh,Phimh,SigmaEtamh,SigmaEpsmh)
    #calculate the acceptance probablity
    Likratio[j]=LikC[j]-Likmh[j]     
    Priorratio[j]=logPriorValue(thetac) - logPriorValue(thetamh[j,,1]) 
    q[4,j]= min(1,exp(Likratio[j]+Priorratio[j]))
    if(is.nan(q[4,j])){
      q[4,j]=-Inf
    }
  } 
  #accept or reject candidate
  u=runif(1,0,1)
  if (u<q[4,j]) {
    thetamh[j,,1]=thetac
    Likmh[j]=LikC[j]
  }
  
  #######################
  #proposal for rhos 
  #######################
  thetac=thetamh[j,,1]
  
  if(S[5,j]<s){
    thetac[Irho]=thetac[Irho]+numhrhohot[j]
  }else{
    thetac[Irho]=thetac[Irho]+numhrhohot[j]
  }
  
  Mumh=matrix(thetac[Imu],n,1)
  Phimh=matrix(thetac[Iphi],n,n)
  SigmaEtamh=diag(thetac[Iseta])    
  SigmaEpsmh=matrix(thetac[Irho],n,n)
  diag(SigmaEpsmh)=1.0
  
  Sigmaev=eigen(SigmaEpsmh)$values
  
  if(abs(thetac[Irho])>0.99 || any(Sigmaev)<0.0) {
    q[5,j]=-1000000.0
    #thetamh[j,,1]=thetamh[j,,1]
    #Likmh[j]=Likmh[j]
  } else {
    #Calculate likelihood with thetac as a parameter vector...        
    #record likelihood of a candidate on step j     
    LikC[j]=likfort2020(n,T,N,Y,Mumh,Phimh,SigmaEtamh,SigmaEpsmh)
    #calculate the acceptance probablity
    Likratio[j]=LikC[j]-Likmh[j]  
    Priorratio[j]=logPriorValue(thetac)-logPriorValue(thetamh[j,,1]) 
    q[5,j]=min(1,exp(Likratio[j]+Priorratio[j]))  
    if(is.nan(q[5,j])){
      q[5,j]=-Inf
    }
  } 
  #accept or reject candidate
  u=runif(1,0,1)
  if (u<q[5,j]) {
    thetamh[j,,1]=thetac
    Likmh[j]=LikC[j]
  }
  
  cat("MC iteration j=",j," : ",thetamh[j,,1],"\n")
  
}
time_spent = proc.time() - ptm   

for(i in 1:9){
  print(accrate(thetamh[1:M,i,1]))
}

results = list("Y"=Y,"thetamh"=thetamh, "vucold"=vucold, "vuhot"=vuhot, "vuphi1hot"=vuphi1hot, "vuphi1cold"=vuphi1cold, "vuphi2hot"=vuphi2hot, "vuphi2cold"=vuphi2cold, "constcold"=constcold, "consthot"=consthot, "seed"=1235, "mu_true"= mutrue, "phi_true" = phitrue, "seta_true" = diag(setatrue), "rho_true" = rhotrue)
save(results, file="results/results_seed1235_19042021.RData")
save(time_spent, file="results/time_spent_seed1235_19042021.RData")
