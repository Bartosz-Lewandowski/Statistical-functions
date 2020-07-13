linear_regression <- function(y,...,r.stud=FALSE) {
  data <- data.frame(...)
  k <- ncol(data)
  z <- colnames(data)
  z <- c("Intercept",z)
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
  ri <- residua/sqrt(S2ee)
  ti <- ri*sqrt(((n-k-2)/(n-k-1-(ri^2))))
  r <- as.vector(abs(ti))
  values <- c()
  indexes <- c()
  for (i in 1:length(r)) {
    if (r[i] > 2.0) {
      values <- c(values,r[i])
      indexes <- c(indexes,i)
      }
  }
  outliers <- data.frame(Index = indexes,
                         Score = values)
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

  if (r.stud==TRUE) {
    cat("\n\n---")
    cat("\nStudentized residual test\nHighest outliers:\n")
    print(format(outliers[order(-outliers$Score),], digit=3), row.names = FALSE)
    }
  
  }

attach(iris)
linear_regression(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width, r.stud = TRUE)
model=lm(Sepal.Length~.-Species, data=iris)
summary(model)
reszty.stud <- rstudent(model)
outlier.iris <- reszty.stud[which.max(abs(reszty.stud))]
outlier.iris
