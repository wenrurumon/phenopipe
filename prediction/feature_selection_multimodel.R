
library(data.table)
library(dplyr)
library(ggplot2)
library(xgboost)

####################################
#Methylation
####################################

rm(list=ls())
setwd('/Users/wenrurumon/Documents/postdoc/gaoyuan_predict')
setwd('methylation')
raw <- lapply(dir(pattern='.csv'),fread)
x1 <- as.matrix(raw[[1]][,-1])
dimnames(x1) <- list(raw[[1]]$V1,colnames(x1))
y1 <- raw[[3]]$AMS.2018
names(y1) <- colnames(x1)
x2 <- as.matrix(raw[[2]][,-1])
dimnames(x2) <- list(raw[[2]]$V1,colnames(x2))
y2 <- raw[[4]]$AMS.1993
names(y2) <- colnames(x2)
X <- cbind(x1,x2)
Y <- c(y1,y2)
X <- t(X[,!is.na(Y),drop=F])
Y <- Y[!is.na(Y)]

#SFS

y <- Y
x <- X
breaks <- cut(1:ncol(x),breaks=20)
x.sel1 <- lapply(unique(breaks),function(i){
  gc()
  xi <- x[,breaks==i,drop=F]
  xi.sd <- apply(xi,2,sd)
  xi.m0 <- colMeans(xi[y==0,])
  xi.m1 <- colMeans(xi[y==1,])
  xi.d0 <- t(t(xi)-xi.m0)
  xi.d1 <- t(t(xi)-xi.m1)
  xi.lda <- (abs(xi.d1)<abs(xi.d0))+0
  xi.d <- abs(xi.m1-xi.m0)/xi.sd
  xi.pp <- colMeans(xi.lda==y)/(max(mean(y),1-mean(y)))
  xi.d>0.1&xi.pp>1.1
})
x.sel1 <- do.call(c,x.sel1)
x <- x[,x.sel1,drop=F]

#SFS2

rlt.sfs <- lapply(1:1000,function(i){
  set.seed(i)
  j <- sample(1:length(y),length(y),replace=T)
  xj <- x[j,]
  yj <- y[j]
  x.sd <- apply(xj,2,sd)
  x.m0 <- colMeans(xj[yj==0,])
  x.m1 <- colMeans(xj[yj==1,])
  x.d0 <- t(t(x)-x.m0)
  x.d1 <- t(t(x)-x.m1)
  x.lda <- (abs(x.d1)<abs(x.d0))+0
  x.d <- abs(x.m1-x.m0)/x.sd
  x.pp <- colMeans(x.lda[-unique(j),]==y[-unique(j)])/max(mean(y[-unique(j)]),1-mean(y[-unique(j)]))
  cbind(d=x.d,pp=x.pp,sel=(x.d>0.1&x.pp>1))
})
rlt.sfs <- -sort(-rowMeans(sapply(rlt.sfs,function(x){x[,3]})))
x <- x[,match(names(rlt.sfs),colnames(x)),drop=F]

#MFS

mvfsi <- function(y,xk){
  sapply(1:100,function(i){
    set.seed(i)
    j <- sample(1:length(y),length(y),replace=T)
    mfile <- data.frame(y=y,xk)
    model <- MASS::lda(y~.,data=mfile[j,,drop=F])
    yj <- max(mean(y[unique(j)]),1-mean(y[unique(j)]))
    yj2 <- max(mean(y[-unique(j)]),1-mean(y[-unique(j)]))
    model.train <- predict(model,mfile[unique(j),,drop=F])$class==y[unique(j)]
    model.test <- predict(model,mfile[-unique(j),,drop=F])$class==y[-unique(j)]
    rlt <- c(mean(model.train),mean(model.test))
    c(rlt,rlt/c(yj,yj2))
  }) %>% t
}

xk <- x[,1,drop=F]
l0 <- mvfsi(y,xk)
history <- list(l0)
for(i in 2:ncol(x)){
  xi <- cbind(xk,x[,i,drop=F])
  li <- mvfsi(y,xi)
  if(mean(li[,2]>l0[,2])>0.5){
    xk <- xi
    l0 <- li
  }
  history <- c(history,list(li))
  print(c(i,ncol(xk),mean(l0[,2])))
}
rlt.mfs <- data.frame(id=1:ncol(x),t(sapply(history,colMeans)))
temp <- melt(rlt.mfs %>% select(id=1,Train=2,Valid=3),id=1)
ggplot() + geom_line(
  data=temp, aes(x=id,y=value*100,colour=variable)
) + labs(x='#Variables Scanned',y='Accurency%',colour='')
save(x,y,xk,history,file='methylation.rda')

####################################
#Expression
####################################

rm(list=ls())
setwd('/Users/wenrurumon/Documents/postdoc/gaoyuan_predict')
setwd('expression')
raw <- lapply(dir(pattern='.csv'),fread)
x1 <- as.matrix(raw[[1]][,-1])
dimnames(x1) <- list(raw[[1]]$X1,colnames(raw[[1]])[-1])
y1 <- raw[[2]]$AMS.2018
names(y1) <- colnames(x1)
x2 <- as.matrix(raw[[3]][,-1])
dimnames(x2) <- list(raw[[3]]$X,colnames(raw[[3]])[-1])
y2 <- raw[[4]]$AMS.1993
names(y2) <- colnames(x2)
X <- cbind(x1,x2)
Y <- c(y1,y2)
X <- t(X[,!is.na(Y),drop=F])
Y <- Y[!is.na(Y)]

#SFS

y <- Y
x <- X
breaks <- cut(1:ncol(x),breaks=20)
x.sel1 <- lapply(unique(breaks),function(i){
  gc()
  xi <- x[,breaks==i,drop=F]
  xi.sd <- apply(xi,2,sd)
  xi.m0 <- colMeans(xi[y==0,])
  xi.m1 <- colMeans(xi[y==1,])
  xi.d0 <- t(t(xi)-xi.m0)
  xi.d1 <- t(t(xi)-xi.m1)
  xi.lda <- (abs(xi.d1)<abs(xi.d0))+0
  xi.d <- abs(xi.m1-xi.m0)/xi.sd
  xi.pp <- colMeans(xi.lda==y)/(max(mean(y),1-mean(y)))
  xi.d>0.1&xi.pp>1.1
})
x.sel1 <- do.call(c,x.sel1)
x <- x[,x.sel1,drop=F]

#SFS2

rlt.sfs <- lapply(1:1000,function(i){
  set.seed(i)
  j <- sample(1:length(y),length(y),replace=T)
  xj <- x[j,]
  yj <- y[j]
  x.sd <- apply(xj,2,sd)
  x.m0 <- colMeans(xj[yj==0,])
  x.m1 <- colMeans(xj[yj==1,])
  x.d0 <- t(t(x)-x.m0)
  x.d1 <- t(t(x)-x.m1)
  x.lda <- (abs(x.d1)<abs(x.d0))+0
  x.d <- abs(x.m1-x.m0)/x.sd
  x.pp <- colMeans(x.lda[-unique(j),]==y[-unique(j)])/max(mean(y[-unique(j)]),1-mean(y[-unique(j)]))
  cbind(d=x.d,pp=x.pp,sel=(x.d>0.1&x.pp>1))
})
rlt.sfs <- -sort(-rowMeans(sapply(rlt.sfs,function(x){x[,3]})))
x <- x[,match(names(rlt.sfs),colnames(x)),drop=F]

#MFS

mvfsi <- function(y,xk){
  sapply(1:100,function(i){
    set.seed(i)
    j <- sample(1:length(y),length(y),replace=T)
    mfile <- data.frame(y=y,xk)
    model <- MASS::lda(y~.,data=mfile[j,,drop=F])
    yj <- max(mean(y[unique(j)]),1-mean(y[unique(j)]))
    yj2 <- max(mean(y[-unique(j)]),1-mean(y[-unique(j)]))
    model.train <- predict(model,mfile[unique(j),,drop=F])$class==y[unique(j)]
    model.test <- predict(model,mfile[-unique(j),,drop=F])$class==y[-unique(j)]
    rlt <- c(mean(model.train),mean(model.test))
    c(rlt,rlt/c(yj,yj2))
  }) %>% t
}

xk <- x[,1,drop=F]
l0 <- mvfsi(y,xk)
history <- list(l0)
for(i in 2:ncol(x)){
  xi <- cbind(xk,x[,i,drop=F])
  li <- mvfsi(y,xi)
  if(mean(li[,2]>l0[,2])>0.5){
    xk <- xi
    l0 <- li
  }
  history <- c(history,list(li))
  print(c(i,ncol(xk),mean(l0[,2])))
}
rlt.mfs <- data.frame(id=1:ncol(x),t(sapply(history,colMeans)))
temp <- melt(rlt.mfs %>% select(id=1,Train=2,Valid=3),id=1)
ggplot() + geom_line(
  data=temp, aes(x=id,y=value*100,colour=variable)
) + labs(x='#Variables Scanned',y='Accurency%',colour='')
save(x,y,xk,history,file='expression.rda')

####################################
#OMIC
####################################

rm(list=ls())
setwd('/Users/wenrurumon/Documents/postdoc/gaoyuan_predict/expression')
load('expression.rda')
temp.expr <- list(x=x,y=y,xk=xk,history=history)
setwd('/Users/wenrurumon/Documents/postdoc/gaoyuan_predict/methylation')
temp.methy <- list(x=x,y=y,xk=xk,history=history)

setwd('/Users/wenrurumon/Documents/postdoc/gaoyuan_predict')
setwd('omic')
raw <- lapply(dir(pattern='.csv'),fread)
raw <- raw[[1]]
X <- as.matrix(raw[,-1:-2])
id <- gsub('2019HA','A19HA',raw$ID)
y <- temp.expr$y[match(id[id%in%names(temp.expr$y)],names(temp.expr$y))]
x <- X[match(id[id%in%names(temp.expr$y)],id),,drop=F]
x <- apply(x,2,as.numeric)
x <- apply(x,2,function(x){ifelse(is.na(x),median(x,na.rm=T),x)})
x <- x[,colnames(x)%in%c(colnames(temp.methy$x),colnames(temp.expr$x)),drop=F]

#SFS2

rlt.sfs <- lapply(1:1000,function(i){
  set.seed(i)
  j <- sample(1:length(y),length(y),replace=T)
  xj <- x[j,]
  yj <- y[j]
  x.sd <- apply(xj,2,sd)
  x.m0 <- colMeans(xj[yj==0,])
  x.m1 <- colMeans(xj[yj==1,])
  x.d0 <- t(t(x)-x.m0)
  x.d1 <- t(t(x)-x.m1)
  x.lda <- (abs(x.d1)<abs(x.d0))+0
  x.d <- abs(x.m1-x.m0)/x.sd
  x.pp <- colMeans(x.lda[-unique(j),]==y[-unique(j)])/max(mean(y[-unique(j)]),1-mean(y[-unique(j)]))
  cbind(d=x.d,pp=x.pp,sel=(x.d>0.1&x.pp>1))
})
rlt.sfs <- -sort(-rowMeans(sapply(rlt.sfs,function(x){x[,3]})))
x <- x[,match(names(rlt.sfs),colnames(x)),drop=F]

#MFS

mvfsi <- function(y,xk){
  sapply(1:100,function(i){
    set.seed(i)
    j <- sample(1:length(y),length(y),replace=T)
    mfile <- data.frame(y=y,xk)
    model <- MASS::lda(y~.,data=mfile[j,,drop=F])
    yj <- max(mean(y[unique(j)]),1-mean(y[unique(j)]))
    yj2 <- max(mean(y[-unique(j)]),1-mean(y[-unique(j)]))
    model.train <- predict(model,mfile[unique(j),,drop=F])$class==y[unique(j)]
    model.test <- predict(model,mfile[-unique(j),,drop=F])$class==y[-unique(j)]
    rlt <- c(mean(model.train),mean(model.test))
    c(rlt,rlt/c(yj,yj2))
  }) %>% t
}
xk <- x[,1,drop=F]
l0 <- mvfsi(y,xk)
history <- list(l0)
for(i in 2:ncol(x)){
  xi <- cbind(xk,x[,i,drop=F])
  li <- mvfsi(y,xi)
  if(mean(li[,2]>l0[,2])>0.5){
    xk <- xi
    l0 <- li
  }
  history <- c(history,list(li))
  print(c(i,ncol(xk),mean(l0[,2])))
}
rlt.mfs <- data.frame(id=1:ncol(x),t(sapply(history,colMeans)))
temp <- melt(rlt.mfs %>% select(id=1,Train=2,Valid=3),id=1)
ggplot() + geom_line(
  data=temp, aes(x=id,y=value*100,colour=variable)
) + labs(x='#Variables Scanned',y='Accurency%',colour='')
save(x,y,xk,history,file='omic.rda')
