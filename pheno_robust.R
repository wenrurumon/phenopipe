
############################
# Def Fucntion
############################

rm(list=ls())
setwd('/Users/wenrurumon/Documents/gaoyuan')
library(data.table)
library(dplyr)
library(fda)
library(MASS)
library(flare)
library(corpcor)
library(ggplot2)
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
cca <- function(A,B){
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
  return(c("chisq_p"=chisq_p,"df"=p*q));
}
qqplot <- function(p_value){
  n = length(p_value);
  exp = -log10((c(1:n)-0.5)/n);
  rgen = -log10(sort(p_value));
  plot(exp,rgen,xlab="-log10(Expect)",ylab="-log10(Real)");
  abline(0,1,col="red")
}
ccap <- function(l1,l2){
  rlt <- sapply(l2,function(x2){
    sapply(l1,function(x1){
      cca(x1,x2)[[1]]
    })
  })
  dimnames(rlt) <- list(names(l1),sapply(l2,function(x){colnames(x)[1]}))
  rlt
}
qpca <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
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
qpca2 <- function(x,p=0.99){
  A <- qpca(x)
  A <- qpca(x,rank=which(A$prop>p)[1])$X
  A
}
library(igraph)
plotnet <- function(x,mode='undirected'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=1)
}
fc <- function(x){
  w<-as.vector(t(x))[t(x)>0]
  x <- graph_from_adjacency_matrix(x>0,mode='undirected')
  fc <- membership(fastgreedy.community(x,weight=w))
  fc[] <- match(fc,unique(fc))
  fc
}
fc2 <- function(x){
  x.mat <- (x<0.01/length(x))+0
  diag(x.mat) <- 0
  x.score <- -log(x)
  # x.score <- log(2^x)
  x.score[x.mat==0] <- 0
  x.score[x.score==Inf] <- max(x.score[x.score!=Inf]*2)
  x.score <- x.score/max(x.score)
  x.g <- graph_from_adjacency_matrix(x.mat,mode='directed')
  E(x.g)$weight <- as.vector(x.score)[x.mat>0]
  x.g <- as.undirected(x.g)
  # plotclust(x.mat,rlt <- fastgreedy.community(x.g)$membership,main=main)
  list(network = x.mat,
       cluster = fastgreedy.community(x.g)$membership)
}
plotclust <- function(x,membership=NULL,main=NULL){
  G <- graph_from_adjacency_matrix(x>0,mode='undirected')
  if(is.null(membership)){membership=rep(1,ncol(x))}
  plot(create.communities(G, membership), 
       # as.undirected(G), 
       as.directed(G),
       layout=layout.kamada.kawai(as.undirected(G)),
       edge.arrow.size=.1,
       vertex.size=.3,
       vertex.label.cex=1,
       edge.width=.1,
       main=main)
}
plotclust2 <- function(p.clust){
  library(networkD3)
  g <- p.clust$network>0
  g <- apply(g,2,function(x){
    names(which(x))
  })
  g2 <- p.clust$cluster
  names(g2) <- names(g)
  g2[] <- sapply(unique(g2),function(i){
    paste(names(which.max(sapply(g[g2==i],length))),collapse=',')
  })[match(g2,unique(g2))]
  tmp <- matrix(0,0,3)
  colnames(tmp) <- c('source','target','value')
  for(i in 1:length(g)){
    tmp <- rbind(tmp,cbind(names(g)[i],names(g)[i],TRUE))
    if(length(g[[i]])>0){tmp <- rbind(tmp,cbind(names(g)[i],g[[i]],TRUE))}
  }
  plink3 <- as.data.frame(tmp)
  pnode2 <- data.frame(name=names(g2),group=g2,size=1)
  # return(apply(pnode2,2,paste))
  plink3$source <- match(plink3$source,pnode2$name)-1
  plink3$target <- match(plink3$target,pnode2$name)-1
  forceNetwork(Links = plink3, Nodes = pnode2, Source = "source",
               Target = "target", Value = "value", NodeID = "name",
               Nodesize = "size",linkColour = "#999",
               radiusCalculation = "Math.sqrt(d.nodesize)+6",
               Group = "group", opacity =20,charge=-10, legend = T
               ,zoom=T,opacityNoHove=100) 
}

#linear fix

linear.fix <- function(X,lambda=1e-15){
  X0 <- 0
  X2 <- apply(X,2,function(x){
    x[is.na(x)] <- mean(x,na.rm=T)
    x
  })
  X2 <- sapply(1:ncol(X2),function(i){
    predict(lm(X2[,i]~X2[,-i]))
  })
  X2[!is.na(X)] <- X[!is.na(X)]
  l <- mean(abs(X2-X0))
  if(l<=lambda){
    return(X2)
  }
  X0 <- X2
}
minmax <- function(x){
  (x-min(x))/(max(x)-min(x))
}
scale2 <- function(x){
  (x-mean(x))/sd(x)
}

############################
# Data
############################

f <- paste0('new.pheno',1:4,'.txt')
f <- lapply(f,read.table,header=T)
for(i in 1:4){
  rownames(f[[i]]) <- f[[i]]$HID
  f[[i]] <- t(f[[i]][,-1])
}
phenos <- f
phenolist <- names(which(table(unlist(lapply(phenos,colnames)))==4))
phenos <- lapply(phenos,function(x){
  x[,match(phenolist,colnames(x))]
})
names(phenos) <- paste0('pheno',1:4)
phenos <- lapply(1:ncol(phenos[[1]]),function(i){
  (cbind(phenos$pheno1[,i],phenos$pheno2[,i],phenos$pheno3[,i],phenos$pheno4[,i]))
})
names(phenos) <- phenolist
phnoes <- phenos[-7]
phenos$LLS <- phenos$headache+phenos$dizziness+phenos$fatigue+phenos$difficulty_sleep+phenos$GI_symptoms
raw <- phenos
phenos <- lapply(phenos,linear.fix)
# p <- lapply(phenos,minmax)[-7]
p <- lapply(phenos,scale2)[-7]
for(i in 1:length(p)){
  colnames(p[[i]]) <- paste0(names(p)[i],1:4)
}

############################
# Cluster RBST
############################

test <- sapply(1:1000,function(i){
  print(i)
  set.seed(i); sel <- sample(1:822,600)
  pi <- lapply(p,function(x){x[sel,]})
  p.cca <- ccap(pi,pi)
  dimnames(p.cca) <- list(names(pi),names(pi))
  p.clust <- fc2(p.cca)
  p.clust$cluster
})
test <- t(test)
colnames(test) <- names(p)
write.csv(test,'cluster_sampling_600.csv')

############################
# Cluster
############################

# p.cca <- ccap(lapply(p,scale),lapply(p,scale))
p.cca <- ccap(p,p)
dimnames(p.cca) <- list(names(p),names(p))
p.clust <- fc2(p.cca)
plotclust2(p.clust)
mp <- lapply(unique(p.clust$cluster),function(i){
  x <- p[which(p.clust$cluster==i)]
  lapply(1:4,function(j){
    sapply(x,function(x){x[,j]})
  })
})
score <- function(x,y){
  mean((x-y)^2)
}
mp.score <- t(sapply(mp,function(x){
  sapply(x,function(xi){
    sapply(x,function(xj){
      score(xi,xj)
    })
  })[,1]
}))
rownames(mp.score) <- sapply(unique(p.clust$cluster),function(i){
  paste(names(which.max(rowSums(p.clust$network[p.clust$cluster==i,,drop=F]))),collapse=',')
})

p.bench <- lapply(1:4,function(j){
  sapply(p,function(x){x[,j]})
})
rlt <- t(t(mp.score) / sapply(p.bench,function(xi){
  sapply(p.bench,function(xj){
    score(xi,xj)
  })
})[,1])
rlt[,1] <- 0
rlt <- melt(rlt)
colnames(rlt) <- c('cluster','stage','value')

hmap <- ggplot(rlt, aes(stage, cluster)) + geom_tile(aes(fill = value),colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")
hmap

paste(rownames(mp.score),
      sapply(unique(p.clust$cluster),function(i){
        paste0(colnames(p.clust$network)[p.clust$cluster==i],collapse=', ')
      }),
      sep=':')
