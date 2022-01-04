############################################################################################
# Calculate Spillover Measures for MSVM
#
#############################################################################################
source("utils/utils.R")
library("coda")
library(xtable)

#load posterior of the parameters
results_name = paste("results_00_","seed1235_","15042021",".RData", sep="")
load(paste("/Users/Yuliya/Dropbox/2dimSVspillovers_newmac/",results_name, sep = ""))

results_name = paste("results_015015_","seed1238_","15042021",".RData", sep="")
load(paste("/Users/Yuliya/Dropbox/2dimSVspillovers_newmac/",results_name, sep = ""))

results_name = paste("results_","seed1236_","15042021",".RData", sep="")
load(paste("/Users/Yuliya/Dropbox/2dimSVspillovers_newmac/",results_name, sep = ""))

n=2 #dimension of the model
M=10000 #length of the MCMC chain
K=5#number of lags
BurnIn=1000

#result_normalized<-t(t(result)/apply(result,2,sum))
var_shar<-array(0,c(n,n,M))
poststats_total<-matrix(c(0),K,5)
poststats_dir_to<-array(0,c(n,5,K))
poststats_dir_from<-array(0,c(n,5,K))

for(j in 1:K){
        for(i in BurnIn:M){
                #specify Phi and Seta
                Phi=matrix(thetamh[i,(n+1):(n+n*n),1],n,n)
                Seta=diag(thetamh[i,(n+n*n+1):(2*n+n*n),1])
                var_shar[,,i]=vsh_sum(varshar(j, Phi, Seta))
        }
        result_normalized<-array(0,c(n,n,M))
        for(s in BurnIn:M){
                result_normalized[,,s]<-t(t(var_shar[,,s])/apply(var_shar[,,s],2,sum))
        }
        
        #total spillover
        total_spill=rep(0,M)
        for(s in BurnIn:M){
                total_spill[s]=total(result_normalized[,,s])
        }
        p_total<-post_summaries(total_spill[BurnIn:M])
        poststats_total[j,]<-c(p_total$mean,p_total$median,p_total$mode,p_total$hpd_l,p_total$hpd_h)
        #measure of directional spillover effect
        dir_spill_from_i<-matrix(c(0),n,M)
        dir_spill_to_i<-matrix(c(0),n,M)
        for(s in BurnIn:M){
                a<-directional(result_normalized[,,s])
                dir_spill_from_i[,s]=a$from_i
                dir_spill_to_i[,s]=a$to_i
        }
        
        p_dir_to_1<-post_summaries(dir_spill_to_i[1,BurnIn:M])
        p_dir_to_2<-post_summaries(dir_spill_to_i[2,BurnIn:M])
        
        poststats_dir_to[1,,j]<-c(p_dir_to_1$mean,p_dir_to_1$median,p_dir_to_1$mode,p_dir_to_1$hpd_l,
                                  p_dir_to_1$hpd_h)
        poststats_dir_to[2,,j]<-c(p_dir_to_2$mean,p_dir_to_2$median,p_dir_to_2$mode,p_dir_to_2$hpd_l,
                                  p_dir_to_2$hpd_h)
        
        p_dir_from_1<-post_summaries(dir_spill_from_i[1,BurnIn:M])
        p_dir_from_2<-post_summaries(dir_spill_from_i[2,BurnIn:M])
        
        poststats_dir_from[1,,j]<-c(p_dir_from_1$mean,p_dir_from_1$median,p_dir_from_1$mode,p_dir_from_1$hpd_l,
                                    p_dir_to_1$hpd_h)
        poststats_dir_from[2,,j]<-c(p_dir_from_2$mean,p_dir_from_2$median,p_dir_from_2$mode,p_dir_from_2$hpd_l,
                                    p_dir_from_2$hpd_h)
        
        
        
}

xtable(poststats_total,digits=4)
xtable(apply(poststats_dir_to,2,c),digits=4)
xtable(apply(poststats_dir_from,2,c),digits=4)
