T_test <- function(X,
                   Y = NULL,
                   mu0 = 0,
                   alpha = 0.05,
                   test = c('one-sample','independent','dependent'),
                   alternative = c("greater","less","two-sided")){
  options(warn=-1)
  if (length(X)==1 | class(X)!="numeric") {
    stop("X must be numeric vector.", call. = FALSE)
  }
  if (length(mu0)!=1 | class(mu0)!="numeric") {
    stop("'mu0' must be single number ", call. = FALSE)
  }
  if (length(alpha)!=1 | class(alpha)!="numeric" | alpha>1 | alpha<0) {
    stop("'alpha' must be number between 0 and 1 ", call. = FALSE)
  }
  if (test=='one-sample' && is.null(Y) | test==1)  {
    m=mean(X)
    s=sd(X)
    n=length(X)
    T=sqrt(n)*(m-mu0)/s
    cat('One sample T test:\n')
    if (alternative == "greater") {
      p=pt(T,n-1, lower.tail=F) 
    } else if (alternative == "less") {
      p=pt(T,n-1) 
    } else if (alternative == "two-sided") {
      p=2*pt(abs(T),n-1, lower.tail=F) 
    } else {
        stop("Incorrect alternative hypothesis entered, enter 'greater','less', or 'two-sided' in 'alternative' parametr!",call. = FALSE) 
      }
    if (p>alpha) {
      cat('True mean is equal to:', mu0 ,'\n')
    } else {
      cat('True mean is not equal to:', mu0,'\n')
    }
  }
  else if (test=='independent' | test==2 | test=='one-sample' && !is.null(Y)) {
    if (length(Y)==1 | class(Y)!="numeric") {
      stop("Y must be numeric vector.", call. = FALSE)
    }
    if(is.numeric(Y)) {
      cat('T test for two independent samples:\n')
      m1=mean(X)
      m2=mean(Y)
      v1=var(X)
      v2=var(Y)
      n1=length(X)
      n2=length(Y)
      T=(m1-m2)/sqrt((n1-1)*v1+(n2-1)*v2)*sqrt(n1*n2*(n1+n2-2)/(n1+n2))
      if (alternative=="greater") {
        p=pt(T,n1+n2-2, lower.tail=F) 
      } else if (alternative=="less") {
        p=pt(T,n1+n2-2) 
      } else if(alternative == "two-sided") {
        p=2*pt(abs(T),n1+n2-2, lower.tail=F) 
      }
        else {
          stop("Incorrect alternative hypothesis entered, enter 'greater','less', or 'two-sided' in 'alternative' parametr!",call. = FALSE) 
        }
      if (p>alpha) {
        cat('True mean of X is equal to true mean of Y.\n')
      } else {
        cat('True mean of X is not equal to true mean of Y.\n')
      }
    }
    else{
      stop("You must input Y vector!", call. = FALSE)
    }
  }
  else if (test=='dependent' | test==3) {
    if(length(X)!=length(Y)) { 
      stop("Not equal length of X and Y",call. = FALSE)
    }
    else {
    cat('T test for two dependent samples:\n')
    D=X-Y
    m=mean(D)
    s=sd(D)
    n=length(D)
    T=sqrt(n)*(m)/s
    p=pt(T,n-1, lower.tail=T) 
    }
    if (p>alpha) {
      cat('True mean of X is equal to true mean of Y.\n')
    } else {
      cat('True mean of X is not equal to true mean of Y.\n')
    }
  }
  else{
    stop("Incorrect test entered, enter 'one-sample','independent', or 'dependent' in 'test' parametr!", call. = FALSE)
  }
  list(Statistic=T, p_value=p)
  
}
