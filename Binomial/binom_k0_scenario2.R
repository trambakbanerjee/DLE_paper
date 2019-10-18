

library(foreach)
library(doParallel)
library(iterators)
library(MASS)
library(kernlab)
library(REBayes)
library(ggplot2)
library(gridExtra)
library(truncdist)
suppressMessages(suppressWarnings(library(CVXR)))

source('ksd_lib_dev.R')
  
#------------------- Parameters ----------------------------

n.grid<- seq(500,5000,by=500)
reps<- 40
#------------------ Initializations --------------------------

h.ksd.grid<- c(0.01,0.02,0.5)
h.tfor.grid<- exp(seq(log(5*10^(-3)),log(5),length.out=10))
h.ksd.est<- matrix(0,nrow=length(n.grid),ncol = reps)
h.ksd.or<-h.ksd.est
h.tfor.est<-h.ksd.est
risk.ksd.odd<- h.ksd.est
risk.ksd.suc<- h.ksd.est
risk.ksdor.odd<- h.ksd.est
risk.ksdor.suc<- h.ksd.est
risk.km.odd<- h.ksd.est
risk.km.suc<- h.ksd.est
risk.gauskm.odd<-h.ksd.est
risk.gauskm.suc<-h.ksd.est
risk.naive.odd<- h.ksd.est
risk.naive.suc<- h.ksd.est
risk.tfor.odd<- h.ksd.est
risk.tfor.suc<- h.ksd.est

cl <- makeCluster(8) # create a cluster with 2 cores
registerDoParallel(cl) # register the cluster

for (N in 1:10){
  n<- n.grid[N]
  m<- rep(10,n)
 
  result<-foreach(r = 1:reps,.packages=c("CVXR","Rfast","Rlinsolve","MASS","REBayes"))%dopar%{
    
    set.seed(10*r)
    q<-runif(n,0,1)
    set.seed(r^2)
    odd<- (q<=0.8)*runif(n,0.5,0.5)+(q>0.8)*rgamma(n,shape=1,scale=2)
    sp<- odd/(1+odd)
    set.seed(r)
    X<- rbinom(n,m,sp)
    
    # Begin estimation of h est and h orc
    ##########################################
    ure.delta<- matrix(0,nrow=length(h.ksd.grid),ncol=1)
    loss.delta<- matrix(0,nrow=length(h.ksd.grid),ncol=1)
    for(i in 1:length(h.ksd.grid)){
      h.ksd<-h.ksd.grid[i]
      out<-tryCatch(ksd_binom_k0_new(X,h.ksd,n,m),
                    error=function(e) NULL)
      step1.p<- out$delt.p
      Xp<-out$xp
      mp<-out$mp
      step1<- out$delt
      out.ure<-ure_binom_new(X,Xp,step1,step1.p,m,mp,h.ksd,0)
      ure.delta[i]<-out.ure$ure*(!is.null(out))+10^12*(is.null(out))
      loss.delta[i]<-sum((step1-odd)^2)*(!is.null(step1))+10^12*(is.null(step1))
    }
    idx.ksd.or<- which(loss.delta == min(loss.delta))
    idx.ksd<- which(ure.delta == min(ure.delta))
    h.ksd.est<-h.ksd.grid[idx.ksd]
    h.ksd.or<- h.ksd.grid[idx.ksd.or]
    
    #--------------- NEB h OR -----------------------
    out.ksd<- ksd_binom_k0_new(X,h.ksd.or,n,m)
    w.hat.or<-out.ksd$w.hat
    delta.ksdor.odd<-out.ksd$delt
    delta.ksdor.suc<- delta.ksdor.odd/(1+delta.ksdor.odd)
    risk.ksdor.odd<-sum((delta.ksdor.odd-odd)^2)
    risk.ksdor.suc<-sum((delta.ksdor.suc-sp)^2)
    
    #--------------- NEB hest -----------------------
    out.ksd<- ksd_binom_k0_new(X,h.ksd.est,n,m)
    w.hat<-out.ksd$w.hat
    delta.ksd.odd<-out.ksd$delt
    delta.ksd.suc<-delta.ksd.odd/(1+delta.ksd.odd)
    risk.ksd.odd<-sum((delta.ksd.odd-odd)^2)
    risk.ksd.suc<-sum((delta.ksd.suc-sp)^2)
    
    #------------- KM ----------------
    out.km<- Bmix(X,m,v=300,collapse=FALSE)
    km.dy<- predict(out.km,X,Loss=2,m)
    km.odd<- km.dy/(1-km.dy)
    risk.km.odd<- sum((km.odd-odd)^2)
    risk.km.suc<- sum((km.dy-sp)^2)
    
    # #--------------- Normal approximation -----------------------
    Y<- asin(sqrt((X+0.25)/(m+0.5)))
    sigma<- sqrt(1/(4*m))
    out.gauss.km<- GLmix(Y,sigma=sigma)
    delta.gauskm.suc<- (sin(out.gauss.km$dy))^2
    delta.gauskm.odd<- delta.gauskm.suc/(1-delta.gauskm.suc)
    risk.gauskm.suc<-sum((delta.gauskm.suc-sp)^2)
    risk.gauskm.odd<-sum((delta.gauskm.odd-odd)^2)
    
    # #--------------- TF Oracle -----------------------
    loss.tfor<- matrix(0,nrow=length(h.tfor.grid),ncol=1)
    for(i in 1:length(h.tfor.grid)){
      h.tfor<-h.tfor.grid[i]
      step1<- log(tf_binom(X,h.tfor,n,m))
      loss.tfor[i]<-sum((step1-log(odd))^2)
    }
    idx.tfor<- which(loss.tfor == min(loss.tfor))
    h.tfor.est<-mean(h.tfor.grid[idx.tfor])
    delta.tfor.odd<- tf_binom(X,h.tfor.est,n,m)
    delta.tfor.suc<- delta.tfor.odd/(1+delta.tfor.odd)
    
    risk.tfor.suc<-sum((delta.tfor.suc-sp)^2)
    risk.tfor.odd<-sum((delta.tfor.odd-odd)^2)
    
    # #------------- Naive --------------------
    naive.suc<- X/m
    naive.odd<- naive.suc/(1.5-naive.suc)
    risk.naive.odd<- sum((naive.odd-odd)^2)
    risk.naive.suc<- sum((naive.suc-sp)^2)
   
    return(list("risk.km.odd"=risk.km.odd,"risk.gauskm.odd"=risk.gauskm.odd,"risk.tf.odd"=risk.tfor.odd,
                "risk.naive.odd"=risk.naive.odd,"risk.ksd.odd"=risk.ksd.odd,
                "risk.ksdor.odd"=risk.ksdor.odd,"risk.km.suc"=risk.km.suc,"risk.gauskm.suc"=risk.gauskm.suc,"risk.tf.suc"=risk.tfor.suc,
                "risk.naive.suc"=risk.naive.suc,"risk.ksd.suc"=risk.ksd.suc,
                "risk.ksdor.suc"=risk.ksdor.suc))
    
  }
  risk.ksd.odd[N,]<-sapply(1:reps, function(i) result[[i]]$risk.ksd.odd)
  risk.ksdor.odd[N,]<-sapply(1:reps, function(i) result[[i]]$risk.ksdor.odd)
  risk.km.odd[N,]<-sapply(1:reps, function(i) result[[i]]$risk.km.odd)
  risk.gauskm.odd[N,]<-sapply(1:reps, function(i) result[[i]]$risk.gauskm.odd)
  risk.tfor.odd[N,]<-sapply(1:reps, function(i) result[[i]]$risk.tf.odd)
  risk.naive.odd[N,]<-sapply(1:reps, function(i) result[[i]]$risk.naive.odd)
  risk.ksd.suc[N,]<-sapply(1:reps, function(i) result[[i]]$risk.ksd.suc)
  risk.ksdor.suc[N,]<-sapply(1:reps, function(i) result[[i]]$risk.ksdor.suc)
  risk.km.suc[N,]<-sapply(1:reps, function(i) result[[i]]$risk.km.suc)
  risk.gauskm.suc[N,]<-sapply(1:reps, function(i) result[[i]]$risk.gauskm.suc)
  risk.tfor.suc[N,]<-sapply(1:reps, function(i) result[[i]]$risk.tf.suc)
  risk.naive.suc[N,]<-sapply(1:reps, function(i) result[[i]]$risk.naive.suc)
  print(N)
  
}
stopCluster(cl)
registerDoSEQ()

risk.est.odd<- cbind(rowMeans(risk.km.odd),rowMeans(risk.gauskm.odd),
                     rowMeans(risk.naive.odd),rowMeans(risk.tfor.odd),
                     rowMeans(risk.ksd.odd),rowMeans(risk.ksdor.odd))

risk.est.suc<- cbind(rowMeans(risk.km.suc),rowMeans(risk.gauskm.suc),
                     rowMeans(risk.naive.suc),rowMeans(risk.tfor.suc),
                     rowMeans(risk.ksd.suc),rowMeans(risk.ksdor.suc))

plotdata.odd<- as.data.frame(c(c(risk.est.odd[,1]/n.grid),c(risk.est.odd[,2]/n.grid),
                           c(risk.est.odd[,5]/n.grid),c(risk.est.odd[,6]/n.grid)))
names(plotdata.odd)<- 'risk'
plotdata.odd$n<- rep(n.grid,4)
plotdata.odd$type<- as.factor(c(rep('KM',length(n.grid)),
                                rep('TF Gauss',length(n.grid)),
                                rep('NEB',length(n.grid)),
                            rep('NEB OR',length(n.grid))))

g1<-ggplot() +
  geom_line(data=plotdata.odd,aes(x=n,y=risk,color=type))+
  geom_point(data=plotdata.odd,aes(x=n,y=risk,shape=type,color=type))+
  xlab(expression(n))+ylab("risk")+theme_bw()+
  theme(legend.position="top",legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

