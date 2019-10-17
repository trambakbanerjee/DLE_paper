
library(foreach)
library(doParallel)
library(iterators)
library(readr)
library(REBayes)
library(ggplot2)
library(gridExtra)
suppressMessages(suppressWarnings(library(CVXR)))

source('ksd_lib_dev.R')

Crime <- read_csv("Book3.csv")

crime_dat<-setNames(data.frame(matrix(ncol = 3, nrow = 2803)), c("id", "X", "Y"))

crime_dat$id<-Crime$uniqueid
crime_dat$X<-Crime$X
crime_dat$Y<-Crime$lam

#-------------- Poisson count data ----------------
X<- as.integer(crime_dat$X)
lam<- crime_dat$Y
n<-length(X)
#------------------ Initializations --------------------------
reps=1
h.ksd.grid<- exp(seq(log(1*10^(-4)),log(0.01),length.out=10))
h.bgr.grid<- seq(0,5,length.out=10)
h.tf.grid<- exp(seq(log(1*10^(-3)),log(5),length.out=10))
p<-0.9
K.cv<-25
# Begin cv for ksd and bgr
##########################################
cl <- makeCluster(8) # create a cluster with 10 cores
registerDoParallel(cl) # register the cluster

cvscore.bgr<- matrix(0,nrow=length(h.bgr.grid),ncol=1)
pp<- p/(1-p)
for(i in 1:length(h.ksd.grid)){
 
  h.bgr<-h.bgr.grid[i]
  result<-foreach(k = 1:K.cv,.packages="CVXR")%dopar%{
    set.seed(k)
    u<-rbinom(n,X,p)
    v<-pp*(X-u)
    v1<-v
    v[v==0]<-0.1
    
    out.bgr<- bgr(u,h.bgr)$del.3

    return(list("bgr"=sum((((out.bgr-v1[order(u)])^2)))))
  }

  cvscore.bgr[i]<-mean((sapply(1:K.cv,function(i) result[[i]]$bgr)))
}
stopCluster(cl)
registerDoSEQ()

idx.bgr<- which(cvscore.bgr == min(cvscore.bgr))

h.bgr.est<-h.bgr.grid[idx.bgr]
##################################################

# Begin estimation of h est and h orc
##########################################
ure.delta<- matrix(0,nrow=length(h.ksd.grid),ncol=1)
loss.delta<- matrix(0,nrow=length(h.ksd.grid),ncol=1)
for(i in 1:length(h.ksd.grid)){
  h.ksd<-h.ksd.grid[i]
  out<-tryCatch(ksd_poiss_k1_new(X,h.ksd,n),
                error=function(e) NULL)
  step1.p<- out$delt.p
  Xp<-out$xp
  step1<- out$delt
  out.ure<-ure_poiss_new(X,Xp,step1,step1.p,h.ksd,1)
  ure.delta[i]<-out.ure$ure*(!is.null(out))+10^12*(is.null(out))
  loss.delta[i]<-sum(((step1-lam)^2)/lam)*(!is.null(step1))+10^12*(is.null(step1))
}
idx.ksd.or<- which(loss.delta == min(loss.delta))
idx.ksd<- which(ure.delta == min(ure.delta))
h.ksd.est<-h.ksd.grid[idx.ksd]
h.ksd.or<- h.ksd.grid[idx.ksd.or]
####################################################

#--------------- NEB h or -----------------------
out.ksd.or<- ksd_poiss_k1_new(X,h.ksd.or,n)
w.hat.or<-out.ksd.or$w.hat
delta.ksd.or<-out.ksd.or$delt
risk.ksd.or<-sum(((delta.ksd.or-lam)^2)/lam)/n

#--------------- NEB h est -----------------------
out.ksd<- ksd_poiss_k1_new(X,h.ksd.est,n)
w.hat<-out.ksd$w.hat
delta.ksd<-out.ksd$delt
risk.ksd<-sum(((delta.ksd-lam)^2)/lam)/n

#--------------- BGR -----------------------
out.bgr<- bgr(X,h.bgr.est)
risk.bgr<-sum(((out.bgr$del.3-lam[order(X)])^2)/lam[order(X)])/n

#------------- KM ----------------
out.km<- Pmix(X,v=300,exposure=rep(1,n))
km.dy<- predict(out.km,X,Loss=2,newexposure=rep(1,n))
risk.km<- sum(((km.dy-lam)^2)/lam)/n

#--------------- Normal approximation -----------------------
Y<- 2*sqrt(X+0.25)
sigma<- rep(1,n)
out.gauss.km<- GLmix(Y,sigma=sigma)
delta.gauss.km<- 0.25*(out.gauss.km$dy^2)
risk.gauss.km<-sum(((delta.gauss.km-lam)^2)/lam)/n

#--------------- Naive Estimator -----------------------
risk.naive<-sum(((X-lam)^2)/lam)/n

# #--------------- TF -----------------------
loss.tfor<- matrix(0,nrow=length(h.tf.grid),ncol=1)
for(i in 1:length(h.tf.grid)){
  h.tf<-h.tf.grid[i]
  step1<- log(tf_poiss(X,h.tf,n))
  loss.tfor[i]<-sum((step1-log(lam))^2)
}
idx.tf<- which(loss.tfor == min(loss.tfor))
h.tf.est<-mean(h.tf.grid[idx.tf])
delta.tf<- tf_poiss(X,h.tf.est,n)
risk.tf<-sum(((delta.tf-lam)^2)/lam)/n