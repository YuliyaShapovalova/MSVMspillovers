setpar<-function(...) par(mar=c(2.5,2.5,0.5,0.5),mgp=c(1.4,0.2,0),
                          cex=0.7,bty="l",las=1,lwd=0.5,tck=0.01,...)
library('invgamma')

SimulateStochVolModel=function(T,MuTrue,PhiTrue,SigEtaTrue,SigEpsTrue,Seed) {
  if (!missing(Seed)) set.seed(Seed)
  n=length(MuTrue)
  epstrue=rmultnorm(T,rep(0,n),SigEpsTrue)
  etatrue=rmultnorm(T,rep(0,n),SigEtaTrue)
  htrue=matrix(0,n,T)
  Y=matrix(0,n,T)
  htrue[,1]=MuTrue+etatrue[,1]
  Y[,1]=exp(htrue[,1]/2)*epstrue[,1]
  
  for (t in 2:T) { 
    htrue[,t]=MuTrue+PhiTrue%*%(htrue[,t-1]-MuTrue)+etatrue[,t]
    Y[,t]=exp(htrue[,t]/2)*epstrue[,t]
  }
  return(list(Y=Y,h=htrue))
}

resampleSystematic<-function(w,N){
        # [ indx ] = resampleSystematic[ w, N]
        # Systematic resampling method for particle filtering. 
        # Author: Tiancheng Li,Ref:
        #T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
        # submit to IEEE Signal Processing Magazine, August 2013
        
        # Input:
        #w    the input weight sequence 
        #       N    the desired length of the output sequence[i.e. the desired number of resampled particles]
        # Output:
        #indx the resampled index according to the weight sequence
        
        #if (nargin ==1){
        #N = length(w)
        #}
        M = length(w)
        w = w / sum(w)
        Q = cumsum(w)
        indx = matrix(c(0),1, N)
        T = seq(0,1-1/N,length=N) + runif(1)/N#
        
        i = 1
        j = 1
        while(i<=N & j<=M){
                while (Q[j] < T[i]){
                        j = j + 1
                }
                indx[i] = j
                i = i + 1
        }
        return(indx)
}
rmultnorm<-function(N,mu,S,Cholesky=FALSE){
  n=length(mu)
  Z=matrix(rnorm(n*N),n,N)
  if (Cholesky) {
    M=c(mu)+S%*%Z
  } else {
    M=c(mu)+t(chol(S))%*%Z
  }
  return(M)
}

rmultnormtest<-function(Cholesky=FALSE){
      
        if (Cholesky) {
                print(1)
        } else {
                print(2)
        }
}

#acceptance rate function

accrate<-function(Likelihoodvector){
  N=length(Likelihoodvector)
  x=sum(ifelse(Likelihoodvector[-1]!=Likelihoodvector[-N],1,0))
  return(x/N)
}

digamma<-function(x,a,b) {
  return(b^a/gamma(a)*x^(-a-1)*exp(-b/x))
}

logdigamma<-function(x,a,b) {
  return(a*log(b)-log(gamma(a))+(-a-1)*log(x))
}

PriorValue=function(theta) {
  diag_seq=seq(from=1,to=n*n,by=n+1)
  nondiag_seq=seq(from=1,to=n*n,by=1)
  nondiag_seq=nondiag_seq[-diag_seq]
  # prior mu
  prior_mu = prod(dnorm(theta[1:n],c(Mumean),c(Musd)))
  # prior phi
  prior_phi = prod(dnorm(theta[Iphi[nondiag_seq]],c(0,0),c(1,1)))*
    prod(dbeta((1+theta[Iphi[diag_seq]])/2, 20, 1.5))
  # prior seta
  prior_seta = prod(digamma(theta[(n+n*n+1):(2*n+n*n)],SigmaEtaAlpha,SigmaEtaBeta))
  # prior seps
  prior_rho = prod(dunif(theta[(2*n+n*n+1):(ncol(combn(n,2)))]))
  
  prior = prior_mu*prior_phi*prior_seta*prior_rho
  return(prior)
}

logPriorValue=function(theta) {
  diag_seq=seq(from=1,to=n*n,by=n+1)
  nondiag_seq=seq(from=1,to=n*n,by=1)
  nondiag_seq=nondiag_seq[-diag_seq]
  # prior mu
  prior_mu = sum(dnorm(theta[1:n],c(Mumean),c(Musd), log=TRUE))
  # prior phi
  prior_phi = sum(dnorm(theta[Iphi[nondiag_seq]],c(0,0),c(1,1), log=TRUE))+
    #sum(dbeta((1+theta[Iphi[diag_seq]])/2, 20, 1.5, log=TRUE))
    sum(dnorm((theta[Iphi[diag_seq]]), 0.9, 0.05, log=TRUE))
  # prior seta
  #prior_seta = sum(logdigamma(theta[(n+n*n+1):(2*n+n*n)],SigmaEtaAlpha,SigmaEtaBeta))
  prior_seta = sum(dinvgamma(theta[(n+n*n+1):(2*n+n*n)],SigmaEtaAlpha,SigmaEtaBeta, log=TRUE))
  # prior seps
  prior_rho = dunif(theta[9], -1, 1, log=TRUE)
  
  prior = prior_mu+prior_phi+prior_seta+prior_rho
  return(prior)
}

LogLikl=function(T,n,N,mu,Phi,Seta,Seps,Y){
  ESS=rep(0,T)
  #bpf implementation
  InvSeps=solve(Seps)
  DetSeps=2*pi*sqrt(det(Seps))
  RSeta=t(chol(Seta))
  muvec<-matrix(c(t(mu)),n,N)
  # sample particles
  hp=muvec+rmultnorm(N,rep(0,n),RSeta,Cholesky=TRUE) #dimension n*N, mean at time t
  omega=(exp(hp/2))
  Yh=Y[,1]/omega
  wh=1/(DetSeps*apply(omega,2,prod))*exp(-0.5*apply((InvSeps%*%Yh)*Yh,2,sum)) #1*N
  wh1=wh/sum(wh)
  ESS[1]=1/sum(wh1^2)
  if(ESS[1]<N/2){
    k=resampleSystematic(wh1,N)
    hp=hp[,k]
    }
  Lik=rep(0,T)
  Lik[1]=log(mean(wh))
  
  for(t in 2:T){
    # predict
    hp=muvec+Phi%*%(hp-muvec)+rmultnorm(N,rep(0,n),RSeta,Cholesky=TRUE)
    omega=(exp(hp/2))
    Yh=Y[,t]/omega
    # update
    wh=1/(DetSeps*apply(omega,2,prod))*exp(-0.5*apply((InvSeps%*%Yh)*Yh,2,sum)) 
    wh1=wh/sum(wh)
    ESS[t]=1/sum(wh1^2)
    # resample if degeneracy criterion is satisfied
    if(ESS[t]<N/2){
      #k=resampleSystematic(wh1,N)
      k=resampleSystematic(wh1,N)
      hp=hp[,k]
      }
    Lik[t]=Lik[t-1]+log(mean(wh))
  }
  return(Lik[T])
}

# Non-blind implementation


LogLikl_apf=function(T,n,N,mu,Phi,Seta,Seps,Y){
  ESS=rep(0,T)
  #apf implementation
  InvSeps=solve(Seps)
  DetSeps=2*pi*sqrt(det(Seps))
  RSeta=t(chol(Seta))
  muvec<-matrix(c(t(mu)),n,N)
  # sample particles
  h=muvec+rmultnorm(N,rep(0,n),RSeta,Cholesky=TRUE) #dimension n*N, mean at time t
  omega=(exp(h/2))
  Yh=Y[,1]/omega
  wh=1/(DetSeps*apply(omega,2,prod))*exp(-0.5*apply((InvSeps%*%Yh)*Yh,2,sum)) #1*N
  wnorm=wh/sum(wh)
  #ESS[1]=1/sum(wnorm^2)
  k=resampleSystematic(wh,N)
  #hp=h[k]+seta^2/2*((Y[1]^2/beta^2)*exp(-h[k])-1)+rnorm(N,0,seta) #dimension n*R
  hp=h[,k]+(Seta^2/2)%*%((Y[,1]^2)*exp(-h[,k])-1)+rmultnorm(N,rep(0,n),RSeta,Cholesky=TRUE) #dimension n*R
  omega=(exp(hp/2))
  #calculate second-stage weights, we do not need to resample second time
  wn=(-Y[,1]^2/2)*(exp(-hp)-exp(-h[,k])*(1-(hp-h[,k])))
  
  Lik=rep(0,T)
  Lik[1]=mean(wn)+log(mean(wh))
  
  for(t in 2:T){
    # predict
    hp=muvec+Phi%*%(hp-muvec)
    omega=exp(hp/2)
    Yh=Y[,t]/omega
    # update
    wh=1/(DetSeps*apply(omega,2,prod))*exp(-0.5*apply((InvSeps%*%Yh)*Yh,2,sum)) 
    wnorm=wh/sum(wh)
    #ESS[t]=1/sum(whnorm^2)
    # resample if degeneracy criterion is satisfied
    k=resampleSystematic(wh,N)
    hp=h[,k]+(Seta^2/2)%*%((Y[,t]^2)*exp(-h[,k])-1)+rmultnorm(N,rep(0,n),RSeta,Cholesky=TRUE)
    omega=(exp(hp/2))
    wn=(-Y[,t]^2/2)*(exp(-hp)-exp(-h[,k])*(1-(hp-h[,k])))
    w=wn
    hp=hp[,k]
    Lik[t]=Lik[t-1]+mean(w)+log(mean(wh))
  }
  return(Lik[T])
}

likfort2020<-function(n,t,l,y,mu,phi,seta,seps){
  dyn.load('utils/likfort2020.so')
  #print(a)
  #likfin=0.0
  retdata <- .Fortran("likfort2020",
                      as.integer(n),
                      as.integer(t),
                      as.integer(l),
                      as.double(y),
                      as.double(mu),
                      as.double(phi),
                      as.double(seta),
                      as.double(seps),
                      #hp=as.double(n*l))$hp
                      likfin=as.double(1))$likfin
  return(retdata)
}

############################################################################################
# Functions for Calculating Spillover Measures for MSVM
#
#############################################################################################
Mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#Function to calculate variance shares
varshar<-function(M, Phi, Seta){
  #here M is "M-steps ahead"
  n=ncol(Phi)
  Psi=array(0,c(n,n,M))
  for(i in 1:M){
    Psi[,,i]=Phi^i
  }
  theta=matrix(c(0),n,n)
  theta=array(0,c(n,n,M))
  e_i=matrix(c(0),1,n)
  e_j=matrix(c(0),1,n)
  for(m in 1:M){
    for(i in 1:n){
      e_i[1,i]=1
      for(j in 1:n){
        e_j[1,j]=1        
        theta[i,j,m]=(sqrt(Seta[i,i])^-1)*(e_i%*%Psi[,,m]%*%Seta%*%t(e_j))^2/
          (e_i%*%Psi[,,m]%*%Seta%*%t(Psi[,,m])%*%t(e_i))
        e_j[1,j]=0
      }
      e_i[1,i]=0
    }
  }
  return(theta)
}

vsh_sum<-function(vsh){
  result<-matrix(0, nrow=n, ncol=n)
  for (i in seq(n)){
    for (j in seq(n)) {
      result[i,j]<-sum(vsh[i,j,])
    }
  }
  return(result)
}

#measure of total spillover
total<-function(result_normalized){
  totalsum=sum(result_normalized)
  diag(result_normalized)=0
  spill=sum(result_normalized)/totalsum
  return(spill)
}

#measure of directional spillover
directional<-function(result_normalized){
  result_normalized_nondiag=result_normalized
  diag(result_normalized_nondiag)<-0
  n=ncol(result_normalized)
  #sum over rows
  to_i=apply(result_normalized_nondiag,1,sum)/apply(result_normalized,1,sum)      
  #sum over cols
  from_i=apply(result_normalized_nondiag,2,sum)/apply(result_normalized,2,sum)            
  return(list(from_i=from_i,to_i=to_i))
}

post_summaries<-function(x){
  a=round(mean(x),4)
  b=round(median(x),4)
  c=round(Mode(x), 4)
  d_l=round(HPDinterval(mcmc(x))[1],4)
  d_h=round(HPDinterval(mcmc(x))[2],4)
  return(list(mean=a,median=b,mode=c,hpd_l=d_l,hpd_h=d_h))
}