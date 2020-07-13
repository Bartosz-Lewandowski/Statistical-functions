statistics = function(y = NULL,
                      y2=NULL,
                      x=NULL,
                      x2=NULL,
                      ...){
  options(warn=-1)
  
  if (length(list(...)) > 0 || !is.null(y2) || !is.null(x) || !is.null(x2)){
    features = TRUE
  } else {
    features = FALSE
  }
  
  if (is.null(x) && !is.null(x2)){
    x <- x2
    x2 <- NULL
  }

# Student T test ----------------------------------------------------------
  T_test <- function(y,
                     y2 = NULL,
                     mu0 = 0,
                     alpha = 0.05,
                     test = c('one-sample','independent','dependent'),
                     alternative = c("greater","less","two-sided")){
    if (length(y)==1 | class(y)!="numeric") {
      stop("y must be numeric vector.", call. = FALSE)
    }
    if (length(mu0)!=1 | class(mu0)!="numeric") {
      stop("'mu0' must be single number ", call. = FALSE)
    }
    if (length(alpha)!=1 | class(alpha)!="numeric" | alpha>1 | alpha<0) {
      stop("'alpha' must be number between 0 and 1 ", call. = FALSE)
    }
    if (test=='one-sample' && is.null(y2) | test==1)  {
      m=mean(y)
      s=sd(y)
      n=length(y)
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
    else if (test=='independent' | test==2 | test=='one-sample' && !is.null(y2)) {
      if (length(y2)==1 | class(y2)!="numeric") {
        stop("y2 must be numeric vector.", call. = FALSE)
      }
      if(is.numeric(y2)) {
        cat('T test for two independent samples:\n')
        m1=mean(y)
        m2=mean(y2)
        v1=var(y)
        v2=var(y2)
        n1=length(y)
        n2=length(y2)
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
          cat('True mean of y is equal to true mean of y2.\n')
        } else {
          cat('True mean of y is not equal to true mean of y2.\n')
        }
      }
      else{
        stop("You must input y2 vector!", call. = FALSE)
      }
    }
    else if (test=='dependent' | test==3) {
      if(length(y)!=length(y2)) { 
        stop("Not equal length of y and y2",call. = FALSE)
      }
      else {
        cat('T test for two dependent samples:\n')
        D=y-y2
        m=mean(D)
        s=sd(D)
        n=length(D)
        T=sqrt(n)*(m)/s
        p=pt(T,n-1, lower.tail=T) 
      }
      if (p>alpha) {
        cat('True mean of y is equal to true mean of y2.\n')
      } else {
        cat('True mean of y is not equal to true mean of y2.\n')
      }
    }
    else{
      stop("Incorrect test entered, enter 'one-sample','independent', or 'dependent' in 'test' parametr!", call. = FALSE)
    }
    cat('\n')
    print(list(Statistic=T, p_value=p))
    
  }
  
  

# Linear Regression -------------------------------------------------------
  linear_regression <- function(y,...) {
    data <- data.frame(...)
    k <- ncol(data)
    n <- length(y)
    X <- matrix(1,n,k+1)
    for (i in 1:k){
      X[, i+1] <- data[,i]
    }
    Y <- as.matrix(y)
    alpha = solve(t(X)%*%X)%*%t(X)%*%Y
    Y2 = X%*%alpha
    residua = Y -Y2
    S2e = sum(residua^2)/(n-k-1)
    S2ee = sum(residua^2)/(n)
    D = S2e*solve(t(X)%*%X)
    Salpha = sqrt(diag(D))
    testT = alpha/Salpha
    p.value = 2*pt(abs(testT), n-k-1, lower.tail=FALSE)
    star <- c(0:k+1)
    for (i in 0:k+1){
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
    R2=1 - sum(residua^2)/sum((y-mean(y))^2)
    Ftest=(R2/(1-R2))*(n-k-1)/(k)
    p.valueF = pf(Ftest, k,n-k-1, lower.tail=FALSE)
    output <- data.frame ('Name'=z,
                          Estimate.Std.=alpha,
                          Std.Error=Salpha,
                          't value'=testT,
                          'Pr'=p.value,
                          Dec = star)
    cat("\nResiduals:\n\n", summary(residua),"\n\n")
    cat("Coefficients:\n")
    print(format(output,justify = "left",digit=2),right = FALSE)
    cat("---\nSignif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n\n")
    cat("Residual standard error:", round(sqrt(S2e), 4),  "on", n-k-1,  "degrees of freedom\nR-squared:", round(R2, 4), "\nF-statistic:", round(Ftest,1), "on", k, "and", n-k-1, "DF, F-statistics p-value:", p.valueF)
    return(p.value)
    
  }
  
  

# Correlation -------------------------------------------------------------
  pearson = function(x,y) {
    n=length(x)
    mX=mean(x)
    mY=mean(y)
    R=sum((x-mX)*(y-mY)) / 
      sqrt(sum((x-mX)^2)*sum((y-mY)^2))
    T=R/sqrt(1-R^2)*sqrt(n-2)
    p=2*pt(abs(T), n-2, lower.tail=FALSE)
    print(list('corelation'=R, 'p-value'=p))
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
    print(list('corelation'=Rs, 'p-value'=p))
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
    print(Tau)
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
  

# ANOVA -------------------------------------------------------------------
  one_way <- function(y,x, post.hoc = FALSE) {
    count <- c()
    for (i in 1:length(unique(x))){
      count <- c(count,as.matrix(table(x))[i])
    }
    n <- length(y)
    k <- length(unique(x))
    srY = mean(y)
    srG = aggregate(y, by = list(x), FUN = "mean")
    colnames(srG) = c("X", "srG")
    SSt = sum((y - srY)^2)
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
  
  two_way <- function(y,x,x2){
    n <- length(y)
    k <- length(unique(x))
    k2 <- length(unique(x2))
    count <- c()
    for (i in 1:k){
      count <- c(count,as.matrix(table(x))[i])
    }
    count2 <- c()
    for(i in 1:k2){
      count2 <- c(count2,as.matrix(table(x2))[i])
    }
    srY = mean(y)
    mAlpha = aggregate(y, by = list(x), FUN = "mean")
    colnames(mAlpha) = c("X", "srAlpha")
    mBetta = aggregate(y, by = list(x2), FUN = "mean")
    colnames(mBetta) = c("X2", "srBetta")
    SSt = sum((y - srY)^2)
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
  
  anova.test <- function(y,x,x2 = NULL, post.hoc=FALSE) {
    if (is.numeric(y) & is.factor(x)) {
      if (is.null(x2)) {
        one_way(y,x, post.hoc)
      }
      else if (is.factor(x2)){
        two_way(y,x,x2)
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

# Program menu ------------------------------------------------------------

  cat('Welcome in statistical program!.\nWhat You want to do?\n\n')
  
value <- ''
while (value != 'exit' && value != 'e') {
  if(!is.null(y)){
  cat('Write "1" if you want to make student t test.\n')
  }
  if (!is.null(y2)){
  cat('Write "2", if you want to make linear regression.\n')
  cat('Write "3", if you want to make correlation test.\n')
  }
  if (!is.null(x) || !is.null(x2)){
  cat('Write "4", if you want to make ANOVA test.\n')
  }
  cat('Write "help" to get help.\n')
  cat('Write "exit", to close program.\n')
  value <- readline()
  if (value == 'exit' | value == 'e') {
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt))
    stop()
  } else if (value == '1'){
    cat('Write "1", to make one-sampled t test.\n')
    if (is.numeric(y2)){
    cat('Write "2", to make t test for independent variables.\n')
    cat('Write "3", to make t test for dependent variables.\n')
    }
    value2 <- readline()
    if (value2 == '1'){
      test = 'one-sample'
      mu0 <- as.numeric(readline(prompt="mu0: "))
      alpha <- as.numeric(readline(prompt="Alfa: "))
      alternative <- readline(prompt="Alternative hypothesis (greater,less,two-sided): ")
      T_test(y, mu0 = mu0, alpha = alpha, test = test, alternative = alternative)
    
    } else if (value2 == '2') {
      if (!is.null(y2)){
      test = 'independent'
      alpha <- as.numeric(readline(prompt="Alfa: "))
      alternative <- readline(prompt="Alternative hypothesis (greater,less,two-sided): ")
      T_test(y, y2, alpha = alpha, test = test, alternative = alternative)
      } else {
        stop("Unknown command!",call. = FALSE)
      }
      
    } else if (value2 == '3') {
      if (!is.null(y2)){
      test = 'dependent'
      alpha <- as.numeric(readline(prompt="Alfa: "))
      T_test(y, y2, alpha = alpha, test = test)
      } else {
        stop("Unknown command!",call. = FALSE)
      }
      
    } else {
      stop("Unknown command!",call. = FALSE)
    }

  } else if (value == '2'){
    if (features == TRUE){
    question = ''
    
    if (is.null(x2) & is.null(x)){
      tmp <- data.frame(y2,...)
    } else if (is.null(x2)){
      tmp <- data.frame(y2,x,...)
    } else if (is.null(x)){
      tmp <- data.frame(y2,x2,...)
    } else {
      tmp <- data.frame(y2,x,x2,...)
    }
    
    z <- c("Intercept",colnames(tmp))
    lin <- linear_regression(y,tmp)
    while (question!='No' | question!='N' | question!='NO' | question!='no'){
      if (dim(tmp)[2]==1){
        break()
      }
      cat("\n\nWould You like to eliminate variable with the highest p-value?\nType 'Yes' or 'No'\n")
      question = readline()
      if(question=='No' | question=='N' | question=='NO' | question=='no'){
        break()
      } else if (question=="Yes" | question=="Y" | question == "YES" | question == "yes"){
        p_value <- lin
        highest <- which.max(p_value)
        tmp <- tmp[,-(highest-1)]
        z <- z[-(highest)]
        lin <- linear_regression(y,tmp)
        if(is.null(dim(tmp))){
          cat("\n\nYou have eliminated every possible variable.\n")
          break()
        }
        
      }
    }
    } else {
      stop("Unknown command!",call. = FALSE)
    }
  } else if (value == '3'){
    if (!is.null(y2)){
    cat("Write '1', to use Pearson's method\n")
    cat("Write '2', to use Spearman's method\n")
    cat("Write '3', to use Kendall's method\n")
    value3 <- readline()
    if(value3 == '1'){
      method = 'pearson'
    } else if(value3 == '2'){
      method = 'spearman'
    } else if(value3 == '3'){
      method = 'kendall'
    } else {
      cat("Unknown option")
    }
    Corelation(y, y2, method)
    } else {
      stop("Unknown command!",call. = FALSE)
    }
  } else if (value == '4'){
    if (!is.null(x) || !is.null(x2)){
    value4 = ''
    cat("Type '1' if you want to make one-way ANOVA.\n")
    if (!is.null(x) && !is.null(x2)){
    cat("Type '2' if you want to make two-way ANOVA. \n")
    }
    value4 <- readline()
    if (value4 == '1'){
      value5 <- ''
      cat ("Do you want to make post-hoc test (Tuley-HSD)?\n")
      cat("Type 'Yes' or 'No'")
      value5 <- readline()
      if (value5 == 'Yes' | value5 == 'Y' | value5 == 'yes' | value5 == 'YES'){
        anova.test(y,x,post.hoc = TRUE)
      }
      else {
        anova.test(y,x,post.hoc = FALSE)
      }
    } else if (value4 == 2){
      if (!is.null(x2)){
      anova.test(y,x,x2,post.hoc = FALSE)
      } else {
        stop("Unknown command!",call. = FALSE)
      }
    } else {
      stop("Unknown command!",call. = FALSE)
    }
    } else {
      stop("Unknown command!",call. = FALSE)
    }
  } else if (value == 'help'| value == 'HELP' | value == 'Help'){
    cat('Welcome!\n')
    cat('To use this function you have to entered variables as argument.\n')
    cat(' y - numerical varible (this variable is required),\n y2 - second numerical variable,\n x - factorial variable,\n x2 - second factorial variable,\n')
    cat('Next you choose what test you want to do.\n')
    cat('To make T test you must enter at least y. If you want to make t test for two variables, you must enter y, and y2.\n')
    cat('To make linear regression you must enrer y, and second variable. Dependent variable is y.\n')
    cat('To make correlation you must enter y, and y2.\n')
    cat('To make anova test you must enter y, x or x2, depending on which (one-way or two-way) test you want to make.\n')
    cat('No matter how many variables you entered it will take only y,y2,x and x2 to tests, only linear regression will take every variable inputed.\n')
  }
   else {
    cat("Please type 1,2,3,4 or exit.\n")
  }
  cat('\n\n')
}
}


x <- rnorm(100, 45, 10)
y <- rnorm(100, 45, 10)
z <- rnorm(100, 45, 10)
f <- iris[, 1]
r <- iris[, 2]
l <- iris[, 3]
p <- iris[, 4]
q <- iris[, 5]

#(y, y2, x, x2, ...)
statistics(f, l, q, x2=NULL, r, p)

