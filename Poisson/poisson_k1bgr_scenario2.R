

library(foreach)
library(doParallel)
library(iterators)
library(MASS)
library(kernlab)
library(REBayes)
library(ggplot2)
library(gridExtra)
library(truncdist)

source('ksd_lib_dev.R')


#------------------- Parameters ----------------------------
n.grid<- 5000
reps<- 10
K.cv<-10
p<- 0.9
#------------------ Initializations --------------------------

h.bgr.grid<- seq(0,5,length.out=5)

h.bgr.est<- matrix(0,nrow=length(n.grid),ncol = reps)
risk.bgr<- matrix(0,nrow=length(n.grid),ncol = reps)

cl <- makeCluster(8) # create a cluster with 8 cores
registerDoParallel(cl) # register the cluster

for (N in 1:length(n.grid)){
  n<- n.grid[N]
 
  for (r in 1:reps){
    
    set.seed(10*r)
    q<-runif(n,0,1)
    set.seed(10*r)
    lam<- (q<=0.75)*rgamma(n,shape=5,scale=1)+(q>0.75)*rgamma(n,shape=10,scale=1)
    set.seed(r)
    X<- rpois(n,lam)
    
    # Begin cv for bgr
    ##########################################
    cvscore.bgr<- matrix(0,nrow=length(h.bgr.grid),ncol=1)
    pp<- p/(1-p)
    for(i in 1:length(h.bgr.grid)){
       h.bgr<-h.bgr.grid[i]
      result<-foreach(k = 1:K.cv,.packages="CVXR")%dopar%{
        set.seed(k)
        u<-rbinom(n,X,p)
        v<-pp*(X-u)
        out.bgr<- bgr(u,h.bgr)$del.3

        return(list("bgr"=sum((((out.bgr-v[order(u)])^2)))))
      }
      cvscore.bgr[i]<-mean((sapply(1:K.cv,function(i) result[[i]]$bgr)))
    }
    idx.bgr<- which(cvscore.bgr == min(cvscore.bgr))
    h.bgr.est[N,r]<-h.bgr.grid[idx.bgr]
    # ##################################################
    #--------------- BGR -----------------------
    out.bgr<- bgr(X,h.bgr.est[N,r])
    risk.bgr[N,r]<-sum(((out.bgr$del.3-lam[order(X)])^2)/lam[order(X)])
    print(r)
  }
  
}
stopCluster(cl)
registerDoSEQ()

risk.est<- rowMeans(risk.bgr)

