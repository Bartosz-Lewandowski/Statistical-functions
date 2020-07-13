one_way <- function(Y,X, post.hoc = FALSE) {
  count <- c()
  for (i in 1:length(unique(X))){
    count <- c(count,as.matrix(table(X))[i])
  }
  n <- length(Y)
  k <- length(unique(X))
  srY = mean(Y)
  srG = aggregate(Y, by = list(X), FUN = "mean")
  colnames(srG) = c("X", "srG")
  SSt = sum((Y - srY)^2)
  SSb = sum((srG$srG-srY)^2*count)
  SSw = SSt - SSb
  q.value <- qtukey(p=0.95,k,n-k)
  F_statistic = (SSb/(k-1))/(SSw/(n-k))
  p.value = pf(F_statistic,k-1,n-k, lower.tail = FALSE)
  star <- c()
  if (p.value > 0.1){
    star <- ""
  } else if (p.value > 0.05 & p.value < 0.1){
    star <- "."
  } else if (p.value > 0.01 & p.value < 0.05){
    star <- "*"
  } else if (p.value > 0.001 & p.value < 0.01){
    star <- "**"
  } else if (p.value < 0.001){
    star <- "***"
  }
  df <- data.frame(Df=c(k,n-k),
                   'Sum Sq'=c(SSb,SSw),
                   'F value'=c(round(F_statistic,2),""),
                   "p.value"=c(round(p.value,2),""),
                   dec=c(star,""))
  row.names(df) <- c('X','Residuals')
  cat('One-way ANOVA\n')
  print(format(df, digit=4))
  cat('\n---\n')
  cat('Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1')
  if (post.hoc==TRUE){
  tukey.test(q.value, SSw, srG, k, n)
  }
}

two_way <- function(Y,X,X2){
  n <- length(Y)
  k <- length(unique(X))
  k2 <- length(unique(X2))
  count <- c()
  for (i in 1:k){
    count <- c(count,as.matrix(table(X))[i])
  }
  count2 <- c()
  for(i in 1:k2){
    count2 <- c(count2,as.matrix(table(X2))[i])
  }
  srY = mean(Y)
  mAlpha = aggregate(Y, by = list(X), FUN = "mean")
  colnames(mAlpha) = c("X", "srAlpha")
  mBetta = aggregate(Y, by = list(X2), FUN = "mean")
  colnames(mBetta) = c("X2", "srBetta")
  SSt = sum((Y - srY)^2)
  SSa = sum((mAlpha$srAlpha - srY)^2*count)
  SSb = sum((mBetta$srBetta - srY)^2*count2)
  SSe = SSt - SSa - SSb
  Fa = (SSa / (k - 1)) / (SSe / (-1 + n - k))
  Fb = (SSb / (k2 - 1)) / (SSe / (-1+n-k))
  p1 = pf(Fa,-1+n-k, k-1, lower.tail = FALSE)
  p2 = pf(Fb,-1+n-k, k2-1, lower.tail = FALSE)
  p.value <- c(p1,p2)
  star <- c(1:2)
  for (i in 1:2){
    if (p.value[i] > 0.1){
      star[i] <- ""
    } else if (p.value[i] > 0.05 & p.value[i] < 0.1){
      star[i] <- "."
    } else if (p.value[i] > 0.01 & p.value[i] < 0.05){
      star[i] <- "*"
    } else if (p.value[i] > 0.001 & p.value[i] < 0.01){
      star[i] <- "**"
    } else if (p.value[i] < 0.001){
      star[i] <- "***"
    }
  }
  df <- data.frame(Df=c(k-1,k2-1,-1+n-k),
                   'Sum Sq'=c(SSa,SSb,SSe),
                   'F value'=c(round(Fa,2),round(Fb,2),""),
                   "p.value"=c(round(p1,2),round(p2,2),""),
                   dec=c(star,""))
  row.names(df) <- c('X','X2','Residuals')
  cat('Two-way ANOVA\n')
  print(format(df,digits = 2))
  cat('\n---')
  cat('\nSignif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1')
}

tukey.test <- function(q.value, SSw, srG, k, n) {
  MSw <- SSw/(n-k)
  tukey.value <- q.value * sqrt(MSw / (n-k))
  means <- as.vector(srG$srG)
  name <- as.vector(srG$X)
  result <- data.frame(diag(NA, nrow=length(name), ncol=length(name)))
  rownames(result) <- name
  colnames(result) <- name

  for (i in 1:length(name)){
    for (j in 1:length(name)){
      if (i != j){
        name1 <- name[i]
        name2 <- name[j]
        diff <- abs(means[i] - means[j])
        if (tukey.value > diff){
          result[name1, name2] <- ""
        } else {
          result[name1, name2] <-  "***"
        }
      }
    }
  }
  cat('\n---')
  cat('\n Tukey HSD result\n')
  print(result)
  cat("\n")
  cat("Tukey Honestly Signficant Difference ", tukey.value, "\n")
}

anova.test <- function(Y,X,X2 = NULL, post.hoc=FALSE) {
  if (is.numeric(Y) & is.factor(X)) {
    if (is.null(X2)) {
      one_way(Y,X, post.hoc)
    }
    else if (is.factor(X2)){
      two_way(Y,X,X2)
      if(post.hoc==TRUE){
        cat('\nPost hoc not avaliable, try to one-way anova.')
      }
    }
    else {
      stop("X2 has to be factor, change variable and try again!", call. = FALSE)
    }
  }
  else {
    stop("Wrong variables entered!",call. = FALSE)
  }
}

data(iris)
m1 = lm(Sepal.Length~Species, data = iris)
summary(aov(m1))
iris$season = rbinom(150,1,prob=0.5)
m2 = lm(Sepal.Length~Species + factor(season), data = iris)
summary(aov(m2))


anova.test(iris$Sepal.Length,iris$Species, post.hoc = TRUE)
anova.test(iris$Sepal.Length,iris$Species,iris$season)
