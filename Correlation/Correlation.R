pearson = function(x,y) {
  n=length(x)
  mX=mean(x)
  mY=mean(y)
  R=sum((x-mX)*(y-mY)) / 
    sqrt(sum((x-mX)^2)*sum((y-mY)^2))
  T=R/sqrt(1-R^2)*sqrt(n-2)
  p=2*pt(abs(T), n-2, lower.tail=FALSE)
  list('corelation'=R, 'p-value'=p)
}
spearman=function(x,y) {
  ranks.s=function(x) {
    n=length(x)
    new.x=sort(x)
    new.x=as.data.frame(new.x)
    colnames(new.x)='obs'
    x=as.data.frame(x)
    x$rank.s=0
    colnames(x)='obs'
    new.n=length(unique(x))
    if (n==new.n) {
      new.x$ranks.s=1:n
      for (i in 1:n) {
        x$rank.s[i]=new.x$ranks.s[new.x$obs==x$obs[i]]
      }   
    } else {
      new.x=cbind(new.x, 1:n)
      colnames(new.x)=c('obs','rank1')
      new.x=as.data.frame(new.x)
      tmp=aggregate(new.x$rank1, by=list(new.x$obs),
                    FUN='mean')
      colnames(tmp)=c('obs','rank.s')
      for (i in 1:n) {
        x$rank.s[i]=tmp$rank.s[tmp$obs==x$obs[i]]
      } 
    }
    return(ranks.s=x$rank.s)
  }
  R=ranks.s(x)
  S=ranks.s(y)
  n=length(x)
  
  Rs=1-6*sum((R-S)^2)/(n*(n^2-1))
  Ts=Rs/sqrt(1-Rs^2)*sqrt(n-2)
  p=2*pt(abs(Ts),n-2,lower.tail=FALSE)
  list('corelation'=Rs, 'p-value'=p)
}
ranks <- function(x,y){
  X2 <- rank(x)
  Y2 <- rank(y)
  ranked <- data.frame(X2, Y2)
  ranked <- ranked[order(ranked$X2), ]
  return(ranked)
}
kendall <- function(x,y) {
  ranked2 <- ranks(x,y)
  concordant <- 0
  discordant <- 0
  for (i in 1:length(ranked2$Y2)) {
    concordant <- concordant + sum(ranked2$Y2[i]<ranked2$Y2[(i+1):length(ranked2$Y2)], na.rm = TRUE)
    discordant <- discordant + sum(ranked2$Y2[i]>ranked2$Y2[(i+1):length(ranked2$Y2)], na.rm = TRUE)
  }
  Tau <- (concordant-discordant)/(concordant+discordant)
  return(Tau)
}

Corelation <- function(x,
                       y,
                       method = c("pearson","spearman","kendall")) {
  options(warn=-1)
  if (!is.numeric(x) | !is.numeric(y)) {
    stop("Wrong X or Y, use only vector!",call. = FALSE)
  }
  if (length(x)!=length(y)){
    stop("Length of variables not equal!", call. = FALSE)
  }
  if (method == "pearson") {
    cat('Pearson corelation:\n')
    pearson(x,y)
  }
  else if (method == "spearman") {
    cat('Spearman corelation:\n')
    spearman(x,y)
  }
  else if (method == "kendall") {
    cat('Kendall corelation:\n')
    kendall(x,y)
  }
  else {
    stop("Choose one of given method!", call. = FALSE)
  }
}

