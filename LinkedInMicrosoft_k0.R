
library(readr)
library(REBayes)
library(ggplot2)
library(gridExtra)
library(truncdist)
suppressMessages(suppressWarnings(library(CVXR)))

source('ksd_lib_dev.R')

LinkedIn_Microsoft <- read_csv('LinkedIn_Microsoft.csv')

dat<-setNames(data.frame(matrix(ncol = 17, nrow = 20702)), c("id", "TT1", "TT2","TT3",
      "TT4","TT5","TT6","TT7","TT8","PT1","PT2","PT3","PT4","PT5","PT6","PT7","PT8"))

dat$id<-LinkedIn_Microsoft$IDLink
for(i in 2:17){
  
  s<-i+8*(i-2)
  e<-i+8*(i-1)
  dat[,i]<-rowSums(LinkedIn_Microsoft[,s:e])
  
}
idx<-(dat$TT1>0)
dat<-dat[idx,]

modeldat<-setNames(data.frame(matrix(ncol = 17, nrow = dim(dat)[1])), c("id", "TT1", "TT2","TT3",
          "TT4","TT5","TT6","TT7","TT8","PT1","PT2","PT3","PT4","PT5","PT6","PT7","PT8"))
modeldat$id<-dat$id
modeldat$TT1<-rep(1,dim(modeldat)[1])

for(i in 1:dim(modeldat)[1]){
  
  modeldat[i,3:17]<-1*(diff(as.numeric(dat[i,2:17]))>0)
}
idx<- (rowSums(modeldat[,10:17])>0)
modeldat<-modeldat[idx,]

#-------------- The Binomial count data ----------------
X<- rowSums(modeldat[,2:9])
Y<- rowSums(modeldat[,10:17])
n<- dim(modeldat)[1]
m<- rep(8,n)
sp<-Y/(m)
odd<- sp/(1-sp)
idx<-which(odd==Inf)
odd[idx]<-max(odd[-idx])

#------------------ Initializations --------------------------
reps=1
h.ksd.grid<- exp(seq(log(5*10^(-3)),log(5),length.out=10))
h.tfor.grid<- exp(seq(log(1*10^(-3)),log(5),length.out=10))

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
risk.ksdor.odd<-sum((delta.ksdor.odd-odd)^2)/n

#--------------- NEB hest -----------------------
out.ksd<- ksd_binom_k0_new(X,h.ksd.est,n,m)
w.hat<-out.ksd$w.hat
delta.ksd.odd<-out.ksd$delt
risk.ksd.odd<-sum((delta.ksd.odd-odd)^2)/n

#------------- KM ----------------
out.km<- Bmix(X,m,v=300,collapse=FALSE)
km.dy<- predict(out.km,X,Loss=2,m)
km.odd<- km.dy/(1-km.dy)
risk.km.odd<- sum((km.odd-odd)^2)/n

# #--------------- Normal approximation -----------------------
XX<- asin(sqrt((X+0.25)/(m+0.5)))
sigma<- sqrt(1/(4*m))
out.gauss.km<- GLmix(XX,sigma=sigma)
delta.gauskm.suc<- (sin(out.gauss.km$dy))^2
delta.gauskm.odd<- delta.gauskm.suc/(1-delta.gauskm.suc)
risk.gauskm.odd<-sum((delta.gauskm.odd-odd)^2)/n

# #------------- Naive --------------------
naive.suc<- X/m+0.5
naive.odd<- naive.suc/(1-naive.suc)
risk.naive.odd<- sum((naive.odd-odd)^2)/n

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
risk.tfor.odd<-sum((delta.tfor.odd-odd)^2)/n



     