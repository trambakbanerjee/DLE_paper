

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
n.grid<- seq(500,5000,by=500)
reps<- 100

#------------------ Initializations --------------------------

h.ksd.grid<- exp(seq(log(1*10^(-4)),log(0.1),length.out=10))
h.tf.grid<- exp(seq(log(5*10^(-3)),log(5),length.out=10))
h.ksd.est<- matrix(0,nrow=length(n.grid),ncol = reps)
h.ksd.or<- h.ksd.est
h.tf.est<- h.ksd.est
risk.ksd<- h.ksd.est
risk.ksd.or<- h.ksd.est
risk.km<- h.ksd.est
risk.naive<-h.ksd.est
risk.gauss.km<-h.ksd.est
risk.tf<-h.ksd.est

cl <- makeCluster(8) # create a cluster with 2 cores
registerDoParallel(cl) # register the cluster

for (N in 1:length(n.grid)){
  n<- n.grid[N]
  w.hat<- matrix(1,nrow=n,ncol=reps)
  w.hat.or<-w.hat
  delta.ksd<- matrix(0,nrow=n,ncol=reps)
  delta.ksd.or<-delta.ksd
  delta.gauss.km<-delta.ksd
  delta.tf<-delta.ksd

  result<-foreach(r = 1:reps,.packages=c("CVXR","Rfast","Rlinsolve","MASS","REBayes"))%dopar%{
    
    set.seed(10*r)
    q<-runif(n,0,1)
    set.seed(10*r)
    lam<- runif(n,0.5,15)
    set.seed(r)
    X<- rpois(n,lam)
    
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
    risk.ksd.or<-sum(((delta.ksd.or-lam)^2)/lam)
    
    #--------------- NEB h est -----------------------
    out.ksd<- ksd_poiss_k1_new(X,h.ksd.est,n)
    w.hat<-out.ksd$w.hat
    delta.ksd<-out.ksd$delt
    risk.ksd<-sum(((delta.ksd-lam)^2)/lam)
   
    #------------- KM ----------------
    out.km<- Pmix(X,v=300,exposure=rep(1,n))
    km.dy<- predict(out.km,X,Loss=2,newexposure=rep(1,n))
    risk.km<- sum(((km.dy-lam)^2)/lam)
    
    #--------------- Normal approximation -----------------------
    Y<- 2*sqrt(X+0.25)
    sigma<- rep(1,n)
    out.gauss.km<- GLmix(Y,sigma=sigma)
    delta.gauss.km<- 0.25*(out.gauss.km$dy^2)
    risk.gauss.km<-sum(((delta.gauss.km-lam)^2)/lam)
    
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
    risk.tf<-sum(((delta.tf-lam)^2)/lam)
    
    #--------------- Naive Estimator -----------------------
    risk.naive<-sum(((X-lam)^2)/lam)
    
    
    return(list("risk.km"=risk.km,"risk.gauss.km"=risk.gauss.km,"risk.tf"=risk.tf,
                "risk.naive"=risk.naive,"risk.ksd"=risk.ksd,
                "risk.ksd.or"=risk.ksd.or))
    
  }
  risk.ksd[N,]<-sapply(1:reps, function(i) result[[i]]$risk.ksd)
  risk.ksd.or[N,]<-sapply(1:reps, function(i) result[[i]]$risk.ksd.or)
  risk.km[N,]<-sapply(1:reps, function(i) result[[i]]$risk.km)
  risk.gauss.km[N,]<-sapply(1:reps, function(i) result[[i]]$risk.gauss.km)
  risk.tf[N,]<-sapply(1:reps, function(i) result[[i]]$risk.tf)
  risk.naive[N,]<-sapply(1:reps, function(i) result[[i]]$risk.naive)
  print(N)
  
}
stopCluster(cl)
registerDoSEQ()
risk.est<- cbind(rowMeans(risk.km),rowMeans(risk.gauss.km),
                 rowMeans(risk.tf),rowMeans(risk.naive),
                 rowMeans(risk.ksd),rowMeans(risk.ksd.or))
plotdata<- as.data.frame(c(c(risk.est[,1]/n.grid),c(risk.est[,2]/n.grid),
                           c(risk.est[,3]/n.grid),c(risk.est[,5]/n.grid),
                           c(risk.est[,6]/n.grid)))
names(plotdata)<- 'risk'
plotdata$n<- rep(n.grid,5)
plotdata$type<- as.factor(c(rep('KM',length(n.grid)),rep('TF Gauss',length(n.grid)),
                            rep('TF OR',length(n.grid)),
                            rep('NEB',length(n.grid)), 
                            rep('NEB OR',length(n.grid))))

library(ggplot2)
g1<-ggplot() +
  geom_line(data=plotdata,aes(x=n,y=risk,color=type))+
  geom_point(data=plotdata,aes(x=n,y=risk,shape=type,color=type))+
  xlab(expression(n))+ylab("risk")+theme_bw()+
  theme(legend.position="top",legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))



