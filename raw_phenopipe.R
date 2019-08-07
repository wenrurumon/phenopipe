
#############################
# Package
#############################

rm(list=ls())
setwd('/Users/wenrurumon/Documents/gaoyuan/phenox')
library(data.table)
library(dplyr)
library(fda)
library(MASS)
library(flare)
library(corpcor)
library(igraph)
library(reshape2)

#Basic Functions
init <- function(x,...) UseMethod('init')
lhead <- function(x){lapply(x,head)}

#############################
# Class presult
#############################

setClass('presult',slot=list(
  name = 'character',
  result = 'list'
))
setMethod('show','presult',function(object){
  if(object@name=='cca'){
    s <- list(do.call(rbind,object@result))
    names(s) <- object@name    
  }
  print(s)
})
as.data.frame.presult <- function(x){
  do.call(rbind,x@result)
}
as.matrix.presult <- function(x){
  acast(do.call(rbind,x@result),xi~xj,value.var='pvalue')
}

#############################
# Class Pheno
#############################

#Class

setClass('pheno',slot=list(
  data = 'array',
  nodes = 'character',
  dim = 'integer',
  clust = 'character'
))
setMethod('show','pheno',function(object){print(class(object))})

#init
init.character <- function(...){
  tmp <- c(...)
  d <- lapply(tmp,fread)
  nodes <- table(unlist(lapply(d,colnames)))
  d <- sapply(
    lapply(d,function(x){
      as.matrix(x)[,match(names(nodes),colnames(x)),drop=F]
    }), 
    identity, simplify="array")
  colnames(d) <- names(nodes)
  dimnames(d)[[3]] <- tmp
  p <- init(d,clust=names(nodes))
  p
}
init.array <- function(x,clust=NULL){
  if(is.null(dimnames(x))){
    nodes <- rep('!Uncoded',3)
  } else {
    nodes <- sapply(dimnames(x),function(i){
      if(is.null(i)){
        return('!Uncoded')
      }else{
        return(paste(i,collapse=', '))  
      }
    })
  }
  if(is.null(clust)){
    clust <- colnames(x)
  }
  new('pheno',
      data=x,
      dim=dim(x),
      nodes=nodes,
      clust=clust)
}
init.pheno <- function(p,data=NULL,nodes=NULL,clust=NULL,
                       X=NULL,Y=NULL,Z=NULL){
  if(!is.null(data)){p@data <- data}
  if(is.null(X)){X <- 1:p@dim[1]}
  if(is.null(Y)){Y <- 1:p@dim[2]}
  if(is.null(Z)){Z <- 1:p@dim[3]}
  p@data <- p@data[X,Y,Z,drop=F]
  if(!is.null(clust)){
    clust <- getclust(clust)[Y]
  } else {
    clust <- p@clust[Y]
  }
  init(p@data,clust=clust)
}

#############################
# Model Feature
#############################

feature <- function(x,...) UseMethod('feature')
feature.pheno <- function(x,clust=NULL){
  if(!is.null(clust)){
    x@clust <- clust
  }
  rlt <- lapply(unique(x@clust),function(i){
    i <- (x@data[,x@clust==i,,drop=F])
    i <- array(i,dim=c(nrow(i),prod(dim(i)[-1])))
    return(i)
  })
  names(rlt) <- unique(x@clust)
  rlt[sapply(rlt,length)>0]
}

kick <- function(x,...) UseMethod('kick')
fill <- function(x,...) UseMethod('fill')
kick.matrix <- function(x,fun=list(
  function(x){colMeans(is.na(x))<1},
  function(x){apply(x,2,var,na.rm=T)>0})){
  if(!is.null(fun)){
    for(f in fun){
      x <- x[,f(x),drop=F]
    }
  }
  x
}
kick.list <- function(x,fun=list(
  function(x){colMeans(is.na(x))<1},
  function(x){apply(x,2,var,na.rm=T)>0})){
  lapply(x,kick,fun=fun)
}
fill.matrix <- function(x){
  apply(x,2,function(i){
    i[is.na(i)] <- median(i,na.rm=T)
    i
  })
}
fill.list <- function(x){
  lapply(x,fill)
}

#PCA

pca <- function(x,...){UseMethod('pca')}
pca1 <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0|rank==ncol(A)){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
    d <- d[d > 1e-8]
  }
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}
pca.matrix <- function(A,prop=0.99,ifscale=TRUE){
  rlt <- pca1(A,rank=0,ifscale=ifscale)
  if(prop==1){
    return(rlt)
  }
  rlt <- pca1(A,rank=which(rlt$prop>=prop)[1],ifscale=TRUE)
  return(rlt)
}
pca.list <- function(X,prop=0.99,ifscale=T){
  lapply(lapply(X,pca,ifscale=ifscale,prop=prop),function(x){x$X})
}

#Prerun
# x <- init(paste0('pheno',1:4,'.csv')) %>% init(Y=-49)
# x %>% feature %>% kick %>% fill %>% lhead %>% head
# x %>% feature %>% kick %>% fill %>% pca(prop=0.8) %>% lhead %>% head
    
#############################
# Model CCA
cca <- function(...) UseMethod('cca')
#############################

p_ginv_sq <- function(X,p){
  X.eigen = eigen(X);
  X.rank = sum(X.eigen$values>1e-8);
  X.value = X.eigen$values[1:X.rank]^(-1*p);
  if (length(X.value)==1){
    D = as.matrix(X.value);
  }else{
    D = diag(X.value);
  }
  rlt = X.eigen$vectors[,1:X.rank] %*% D %*% t(X.eigen$vectors[,1:X.rank]);
  return(rlt);
}
mrank <- function(X){
  X.svd = svd(X);
  X.rank = sum(X.svd$d>1e-6);
  return(X.rank);
}
mrank_sq <- function(X){
  X.eigen = eigen(X);
  X.rank = sum(Re(X.eigen$values)>1e-6);
  return(X.rank);
}
CCA_chisq_test <- function(rho,n,p,q){
  tstat = -1*n*sum(log(1-rho^2));
  p_value = pchisq(tstat,(p*q),lower.tail=FALSE);
  return(p_value);          
}
ccap <- function(A,B){
  n = nrow(A);
  p = mrank(A);
  q = mrank(B);
  if (p <= q){
    X = A;
    Y = B;
  }else{
    X = B;
    Y = A;
  }
  R = p_ginv_sq(cov(X),0.5) %*% cov(X,Y) %*% p_ginv_sq(cov(Y),1) %*% cov(Y,X) %*% p_ginv_sq(cov(X),0.5);
  k = mrank_sq(R);
  d = Re(eigen(R)$values);
  rho = d[1:k]^(0.5);
  rho[rho >= 0.9999]=0.9;
  chisq_p = CCA_chisq_test(rho,n,p,q);
  return(c(pvalue=chisq_p,df=p*q));
}
cca.list <- function(...){
  tmp <- list(...)
  if(length(tmp)==2){
    x1 <- tmp[[1]]
    x2 <- tmp[[2]]
  } else {
    x1 <- x2 <- tmp[[1]]
  }
  n1 <- names(x1)
  n2 <- names(x2)
  rlt <- lapply(x1,function(i){
    t(sapply(x2,function(j){
      ccap(i,j)
    }))
  })
  rlt <- lapply(1:length(rlt),function(i){
    data.table(xi=n1[i],xj=rownames(rlt[[i]]),rlt[[i]])
  })
  names(rlt) <- n1
  rlt
  # new('presult',name='cca',result=rlt)
}
cca.pheno <- function(...){
  tmp <- lapply(list(...),function(x){
    x %>% feature %>% kick %>% fill
  })
  if(length(tmp)==2){
    x1 <- tmp[[1]] 
    x2 <- tmp[[2]] 
  } else {
    x1 <- x2 <- tmp[[1]]
  }
  new('presult',name='cca',result=(cca(x1,x2)))
}

#Prerun

# rlt_cca <- cca(init(x))

#############################
# Model Clust
clust <- function(x,k,model) UseMethod('clust')
getclust <- function(x) UseMethod('getclust')
getclust.character <- function(x){x}
getclust.numeric <- function(x){paste(x)}
getclust.clust <- function(x){x@clust}
setClass('clust',slot=list(model='character',data='matrix',result='list',clust='character',score='numeric'))
#############################

fc <- function(x){
  w<-as.vector(t(x))[t(x)>0]
  x <- graph_from_adjacency_matrix(x>0,mode='undirected')
  fc <- membership(fastgreedy.community(x,weight=w))
  fc[] <- match(fc,unique(fc))
  fc
}

plotnet <- function(x,mode='undirected'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=1)
}

plotclust <- function(x,...) UseMethod('plotclust')
plotclust.matrix <- function(x,membership=NULL){
  diag(x) <- 0
  G <- graph_from_adjacency_matrix(x>0,mode='undirected')
  if(is.null(membership)){membership=rep(1,ncol(x))}
  plot(create.communities(G, membership), 
       as.directed(G),
       layout=layout.kamada.kawai(as.undirected(G)),
       edge.arrow.size=0,
       vertex.size=.3,
       vertex.label.cex=1,
       edge.width=.1)
}
plotclust.clust <- function(x){
  rlt <- x
  membership <- getclust(x)
  membership <- match(membership,unique(membership))
  x <- (x@data < (0.05/(length(x@data))))
  diag(x) <- 0
  G <- graph_from_adjacency_matrix(x>0,mode='undirected')
  plot(create.communities(G, membership), 
       as.directed(G),
       layout=layout.kamada.kawai(as.undirected(G)),
       edge.arrow.size=0,
       vertex.size=.3,
       vertex.label.cex=1,
       edge.width=.1)
  return(rlt)
}

fc <- function(x){
  w<-as.vector(t(x))[t(x)>0]
  x <- graph_from_adjacency_matrix(x>0,mode='undirected')
  fc <- membership(fastgreedy.community(x,weight=w))
  fc[] <- match(fc,unique(fc))
  fc
}
fci <- function(x){
  fc0 <- fc(x)
  if(length(unique(fc0))<=2){
    return(fc0)
  }
  s <- sapply(unique(fc0),function(i){
    mean(x[fc0==i,fc0!=i])
  })
  s <- ifelse(is.na(match(fc0,(unique(fc0)[s==min(s)]))),0,match(fc0,(unique(fc0)[s==min(s)])))+1
  names(s) <- names(fc0)
  return(s)
}
fc2 <- function(x,fc0){
  n <- colnames(x)
  loss <- sapply(unique(fc0),function(i){
    xi <- x[fc0==i,fc0==i,drop=F]
    fc1 <- fci(xi)
    if(length(fc1)==1){
      return(Inf)
    } else {
      1-sum(xi[outer(fc1,fc1,'==')])/sum(xi)
    }
  })
  xs <- lapply(unique(fc0),function(i){
    xi <- x[fc0==i,fc0==i,drop=F]
    fci(xi)
  })
  for(i in which(loss!=min(loss))){
    ni <- names(xs[[i]])
    xs[[i]] <- rep(1,length(xs[[i]]))
    names(xs[[i]]) <- ni
  }
  for(i in 2:length(xs)){
    xs[[i]] <- xs[[i]] + max(xs[[i-1]])
  }
  unlist(xs)[match(n,names(unlist(xs)))]
}
fc3 <- function(x,k){
  fc0 <- fci(x)
  ki <- length(unique(fc0))
  while(ki<k){
    fc0 <- fc2(x,fc0)
    ki <- length(unique(fc0))
  }
  return(fc0)
}

clust <- function(x,k,model){
  x <- as.matrix(x)
  if(model=='hclust'){
    r <- hclust(as.dist(x))
    clust <- cutree(r,k=k)
  } else if(model=='kmeans'){
    r <- kmeans(x,k)
    clust <- r$cluster
  } else if(model=='netclust'){
    r <- fc3((x<(0.05/length(x))),k=k)
    clust <- r
  }
  clust <- sapply(1:k,function(i){
    paste(names(clust)[clust==i],collapse=',')
  })[clust]
  score <- mean(x[outer(clust,clust,'==')])/mean(x)
  r <- list(result=list(r),clust=clust,score=score)
  return(new('clust',model=model,data=x,result=r$result,clust=paste(r$clust),score=r$score))
}

#Prerun
# rlt_hclust <- sapply(sapply(2:20,clust,x=rlt_cca,model='hclust'),function(x){x@score}) %>% plot.ts
# rlt_kmeans<- sapply(sapply(2:20,clust,x=rlt_cca,model='kmeans'),function(x){x@score}) %>% lines(col=2)
# rlt_netclust<- sapply(sapply(2:20,clust,x=rlt_cca,model='netclust'),function(x){x@score}) %>% lines(col=3)

########################################################################
########################################################################

##################
# PhenoPipe
##################

#input data and feature selection
pipe <- init('pheno1.csv','pheno2.csv','pheno3.csv','pheno4.csv') %>% init(Y=-49)
clust_cca <- cca(pipe) %>% clust(k=8,model='netclust') %>% plotclust

#cluster the phenoytpes with submodule detection model
pipe <- pipe %>% init(clust=clust_cca)
pipe %>% feature %>% kick %>% fill %>% pca(prop=0.8) 



