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