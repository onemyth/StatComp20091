## ----random, echo=TRUE--------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
plot(lm.D9)

## ----table, echo=TRUE---------------------------------------------------------
knitr::kable(head(iris))

## -----------------------------------------------------------------------------
n<-1000
uy<-runif(n)
ux<-2/((1-uy)^(1/2))#F(x)=1-(2/x)^2
hist(ux,prob=TRUE,main=expression(f(x)==8/x^3),freq=F,breaks=100)
y<-seq(0.01,50,0.1)
lines(y,8/y^3,col="red")

## -----------------------------------------------------------------------------
n<-1000
u1<-runif(n,min=-1,max=1)#u1
dim(u1)<-c(1,n)
u2<-runif(n,min=-1,max=1)#u2
dim(u2)<-c(1,n)
u3<-runif(n,min=-1,max=1)#u3
dim(u3)<-c(1,n)
fy<-array(0,dim=c(1,n))
i=1
while (i<1001) {
   a=abs(u1[1,i])
   b=abs(u2[1,i])
   c=abs(u3[1,i])
   if ((c>=b) & (c>=a))
      fy[1,i]<-u2[1,i]
   else
      fy[1,i]<-u3[1,i]
   i<-i+1
}
hist(fy,prob=T,main=expression(f(x)==3/4(1-x^2)))
y<-seq(-1,1,0.01)
lines(y,3*(1-y^2)/4,col="red")

## -----------------------------------------------------------------------------
n<-1000
uy<-runif(n)
ux<-2/((1-uy)^(1/4))-2#F(x)=1-(2/y)^2
hist(ux,prob=TRUE,main=expression(f(x)==64/(2+x)^5),freq=F)
y<-seq(0,50,0.1)
lines(y,64/(2+y)^5,col="red")

## -----------------------------------------------------------------------------
n<-1e5;
x<-runif(n,min=0,max=pi/3);
Ex<-mean(sin(x))*pi/3;
print(c(Ex,1-cos(pi/3)));

## -----------------------------------------------------------------------------
n<-1e5
ux<-runif(n/2,min=0,max=1);
uy<-runif(n,min=0,max=1);
ex<-(exp(ux)+exp(1-ux))/2;# antithetic variate approach-s1
ey<-exp(uy);#  simple Monte Carlo method-s2
meex<-mean(ex);# the mean of s1
meey<-mean(ey);# the mean of s2
vaex<-var(ex);# the variance of s1
vaey<-var(ey);# the variance of s2
pr<-(vaey-vaex)/vaey; # empirical estimate of the percent reduction in variance
print(c(meex,meey,vaex,vaey,pr));

## -----------------------------------------------------------------------------
x<-seq(1,20,by=0.2);
y1<-x^2*exp(-x^2/2)/sqrt(2*pi);
y2<-exp(1-x);
y3<-7/x^8;
plot(x,y1,type="l",col = "black",xlim=c(0,22),ylim = c(0, 1.2));
lines(x,y2,type="l",col = "blue");
lines(x,y3,type="l",col = "red")
legend(12,1,c("The function","f1","f2"),col=c("black","blue","red"),text.col=c("black","blue","red"),lty=c(1,1,1))

## -----------------------------------------------------------------------------
n=10000
t<-runif(n,0,1)
x1<-1-log(1-t)
x2<-1/(1-t)^(1/7)
y11<-(x1)^2*exp(-(x1)^2/2)/sqrt(2*pi)
y12<-(x2)^2*exp(-(x2)^2/2)/sqrt(2*pi)
y2<-exp(1-(x1))
y3<-7/(x2)^8
m1<-mean(y11/y2)
m2<-mean(y12/y3)
sd1<-sd(y11/y2)
sd2<-sd(y12/y2)
s=0.400626#True integral value
print(c(s,m1,m2,sd1,sd2))

## -----------------------------------------------------------------------------
n=10000#The number of random numbers
k=5#Number of intervals
r=n/k#Number of iterations per cell
N=50#Number of iterations
T2<- numeric(k)
T3<-numeric(k)
estimates <- matrix(0, N, 1)

for (i in 1:N) {
  for (j in 1:k){
    t<-runif(r,0,1)
    x<--log(exp(-(j-1)/5)-(exp((1-j)/5)-exp(-j/5))*t)
    T2[j]<-mean((exp((1-j)/5)-exp((-j)/5))/(1+x^2))
  }
  estimates[i, 1] <- sum(T2)
}
ss<-mean(estimates)
sd1=sd(estimates)
s=0.5247971#True integral value
print(c(ss,s,sd1))

## -----------------------------------------------------------------------------
m=1e5
n=0
k=9
u=1
v=1
for (i in c(1:m)){
  x=rlnorm(k,u,v)
  y=log(x)
  mm=mean(y)
  ss=sd(y)
  l1=mm-0.7687*ss
  l2=mm+0.7687*ss
  if (l1<u && u<l2)
  {
    n=n+1
  }
}
ECP=n/m
print(c('The empirical confidence level ECP is:',ECP))
print(c('The actual confidence level CL is:',0.95))

## -----------------------------------------------------------------------------
m=1e5
n=20
k=0
u=2
v=2
for (i in c(1:m)){
  x=rchisq(n,u)
  mm=mean(x)
  ss=sd(x)
  l1=mm-0.468*ss
  l2=mm+0.468*ss
  if (l1<u && u<l2)
  {
    k=k+1
  }
}
ECP=k/m
print(c('The empirical confidence level ECP is:',ECP))
print(c('The actual confidence level CL is:',0.95))

## -----------------------------------------------------------------------------
sk<-function(x){
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  return(m3/m2^1.5)
}

alpha<-0.1
b<-c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)
n<-40#Sample size
m<-2500#Number of repetitions
cv<-qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
sktests<-numeric(m)
pwr<-numeric(length(b))
for (j in 1:length(b)){
for (i in 1:m){
  x<-rbeta(n,b[j],b[j])
  sktests[i]<-as.integer(abs(sk(x)>=cv))
  }
pwr[j]<-mean(sktests)
}
plot(b,pwr,type="b",xlab=bquote(b),ylim=c(0,0.05))
abline(h=.1,lty=3)
se<-sqrt(pwr*(1-pwr)/m)
lines(b,pwr+se,lty=3)
lines(b,pwr-se,lty=3)

## -----------------------------------------------------------------------------
sk<-function(x){
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  return(m3/m2^1.5)
}

alpha<-0.1
b<-c(seq(1,100,by=1))
n<-40#Sample size
m<-2500#Number of repetitions
cv<-qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
sktests<-numeric(m)
pwr<-numeric(length(b))
for (j in 1:length(b)){
   for (i in 1:m){
    x<-rt(n,df=b[j])
    sktests[i]<-as.integer(abs(sk(x)>=cv))
   }
   pwr[j]<-mean(sktests)
}
plot(b,pwr,type="b",xlab=bquote(b),ylim=c(0,0.6))
abline(h=.1,lty=3)
se<-sqrt(pwr*(1-pwr)/m)
lines(b,pwr+se,lty=3)
lines(b,pwr-se,lty=3)

## -----------------------------------------------------------------------------
count5test<-function(x,y){
  X<-x-mean(x)
  Y<-y-mean(y)
  outx<-sum(X>max(Y))+sum(X<min(Y))
  outy<-sum(Y>max(X))+sum(Y<min(X))
  return(as.integer(max(c(outx,outy))>5))
}

a<-0.055
sigma1<-1
sigma2<-1.5
n<-c(10,100,5000)
k<-length(n)
m<-3000
power<-array(2*k, c(k, 2))
for (i in 1:k){
  power[i,1]<-mean(replicate(m,expr={
    x1<-rnorm(n[i],0,sigma1)
    y1<-rnorm(n[i],0,sigma2)
    count5test(x1,y1)
  }))
  
  power[i,2]<-mean(replicate(m,expr={
    x2<-rnorm(n[i],0,sigma1)
    y2<-rnorm(n[i],0,sigma2)
    s<-var(x2)/var(y2)
    l1<-qf(1-a/2,n[i]-1,n[i]-1)
    l2<-1/l1
    as.numeric((s<l2)|(s>l1))
  }))
  
}

print(power)

## -----------------------------------------------------------------------------
rmvn.eigen <- function(n, mu, Sigma) {
  # generate random vectors from MVN(mu, Sigma)
  # dimension is inferred from mu and Sigma
  d <- length(mu)
  ev <- eigen(Sigma, symmetric = TRUE)
  lambda <- ev$values
  V <- ev$vectors
  C <- V %*% diag(sqrt(lambda)) %*% t(V)
  Z <- matrix(rnorm(n*d), nrow = n, ncol = d)
  X <- Z %*% C + matrix(mu, n, d, byrow = TRUE)
  X
}

alpha<-0.05
cv<-qchisq(1-alpha,df=4)
m<-1000#Number of repetitions for different sample sizes
A<-array(c(5,3,3,7),dim=c(2,2))
n <- c(10, 20, 30, 50, 100, 500)
mu<-array(c(0,0),dim=c(1,2))
sv<-numeric(length(n))
for (i in 1:length(n)){
    cc<-numeric(m)
    for (k in 1:m){
      zX<-numeric(2*n[i])
      dim(zX)<-c(n[i],2)
      X<-numeric(n[i]*2)
      dim(X)<-c(n[i],2)
      X<-rmvn.eigen(n[i],mu,A)
      zX[,1]<-X[,1]-mean(X[,1])
      zX[,2]<-X[,2]-mean(X[,2])
      la<-sum(zX[,1]*zX[,1])/(n[i]-1)
      lb<-sum(zX[,2]*zX[,2])/(n[i]-1)
      lc<-sum(zX[,1]*zX[,2])/(n[i]-1)
      L<-solve(array(c(la,lc,lc,lb),dim=c(2,2)))
      dim(zX)<-c(n[i],2)
      s<-zX%*%L%*%t(zX)
      count<-sum(s^3)/n[i]/6
      cc[k]<-as.integer(count>cv)
    }
  sv[i]<-mean(cc)
} 
    
for (d in 1:length(n)) {
  
  print(paste0("when n=",n[d],",the empirical estimate for the first type of error rate is:" ,sv[d]))
  
}

## -----------------------------------------------------------------------------
rmvn.eigen <- function(n, mu, Sigma) {
  # generate random vectors from MVN(mu, Sigma)
  # dimension is inferred from mu and Sigma
  d <- length(mu)
  ev <- eigen(Sigma, symmetric = TRUE)
  lambda <- ev$values
  V <- ev$vectors
  C <- V %*% diag(sqrt(lambda)) %*% t(V)
  Z <- matrix(rnorm(n*d), nrow = n, ncol = d)
  X <- Z %*% C + matrix(mu, n, d, byrow = TRUE)
  X
}

alpha<-0.1
cv<-qchisq(1-alpha,df=4)
m<-1000#Number of duplications of different Epsilon value samples
B<-numeric(8)
dim(B)<-c(2,2,2)
B[,,1]<-array(c(5,3,3,7),dim=c(2,2))
B[,,2]<-array(c(60,1,1,10),dim=c(2,2))
n <-30
epsilon<-c(seq(0,1,0.05))
mu<-array(c(0,0),dim=c(1,2))
sv<-numeric(length(epsilon))
for (i in 1:length(epsilon)){
    e<-epsilon[i]
    cc<-numeric(m)
    for (k in 1:m){
      zX<-numeric(2*n)
      dim(zX)<-c(n,2)
      X<-numeric(n*2)
      dim(X)<-c(n,2)
      a<-sample(c(1,2),replace=TRUE,size=n,prob=c(e,1-e))
      for (j in 1:n){
       A<-B[,,a[j]]
       X[j,]<-rmvn.eigen(1,mu,A)
      }
       zX[,1]<-X[,1]-mean(X[,1])
       zX[,2]<-X[,2]-mean(X[,2])
       la<-sum(zX[,1]*zX[,1])/(n-1)
       lb<-sum(zX[,2]*zX[,2])/(n-1)
       lc<-sum(zX[,1]*zX[,2])/(n-1)
       L<-solve(array(c(la,lc,lc,lb),dim=c(2,2)))
       dim(zX)<-c(n,2)
       s<-zX%*%L%*%t(zX)
       count<-sum(s^3)/n/6
       cc[k]<-as.integer(count>cv)
    }
  sv[i]<-mean(cc)
} 

plot(epsilon,sv,type="b",xlab=bquote(epsilon),ylim=c(0,1))
abline(h=0.1,lty=3)

## -----------------------------------------------------------------------------
library(bootstrap)
n<-dim(law)[1]
dd<-cor(law[,1],law[,2])
d<-numeric(n)
for (i in 1:n){
  A<-law[(1:n)[-i],]
  d[i]<-cor(A[,1],A[,2])
}
bias<-(n-1)*(mean(d)-dd)
sD<-sqrt((n-1)*mean((d-mean(d))^2))
print(paste0("The bias of the correlation coefficient is:",bias,",The standard deviation of the correlation coefficient is:",sD))

## -----------------------------------------------------------------------------
library(boot)
#bootstrap
rn<-2000
P<-aircondit$hours
dim(P)<-c(dim(aircondit)[1],1)
n<-nrow(P)
r<-numeric(n) 
for (i in 1:n){
  a<-sample(1:n,size=n,replace=TRUE)
  A<-P[a]
  r[i]<-mean(A)
}

#compute the confidence interval
stat <- function(dat, index) {
  dat[index]}

boot.obj <-boot(r, statistic = stat, R=rn)
boot.ci(boot.obj, conf=0.95,type=c("norm","basic","perc", "bca"))

## -----------------------------------------------------------------------------
data(scor, package = "bootstrap")
s<-cov(scor)
l<-eigen(s)$values
l1<-max(l)
theta.hat<-l1/sum(l)
n <- nrow(scor)
print (theta.hat)
#jackknife calculate the bias and standard deviation
theta.jack <- numeric(n)
for (i in 1:n){
   JA<-cov(scor[-i,])
   g<-eigen(JA)$values
   g1<-max(g)
   theta.jack[i]<-g1/sum(g)}
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
se<-sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))
print(paste0("The bias is:",bias,",the standard deviation is:",se)) 

## -----------------------------------------------------------------------------
library(DAAG)
l<-length(ironslag$chemical)
D<-numeric(l*2)
dim(D)<-c(l,2)
X<-ironslag$chemical
Y<-ironslag$magnetic
e1<-e2<-e3<-e4<-numeric(l*(l-1))
n<-1
for (i in 1:l){
  k<-(1:l)[-i]
  for (j in k){
    index<-c(i,j)
    x<-X[-index]
    y<-Y[-index]
    
    L1<-lm(y~x)
    yhat1 <- L1$coef[1] + L1$coef[2] * X[index]
    e1[n] <- sum(Y[index] - yhat1)^2/2
    
    L2 <- lm(y ~ x + I(x^2))
    yhat2 <- L2$coef[1] + L2$coef[2] * X[index] +
      L2$coef[3] * X[index]^2
    e2[n] <- sum(Y[index] - yhat2)^2/2
    
    L3 <- lm(log(y) ~ x)
    logyhat3 <- L3$coef[1] + L3$coef[2] * X[index]
    yhat3 <- exp(logyhat3)
    e3[n] <- sum(Y[index] - yhat3)^2/2
    
    L4 <- lm(log(y) ~ log(x))
    logyhat4 <- L4$coef[1] + L4$coef[2] * log(X[index])
    yhat4 <- exp(logyhat4)
    e4[n] <- sum(Y[index] - yhat4)^2/2
    n<-n+1
  }
}
 
print(c(mean(e1),mean(e2),mean(e3),mean(e4))) 

## -----------------------------------------------------------------------------
count5test<-function(x,y){
  X<-x-mean(x)
  Y<-y-mean(y)
  outx<-sum(X>max(Y))+sum(X<min(Y))
  outy<-sum(Y>max(X))+sum(Y<min(X))
  return(as.integer(max(c(outx,outy))>5))
}

n1<-5
n2<-15
mu1<-mu2<-0
sigma1<-sigma2<-1
m<-1000#Number of repetitions
x<-rnorm(n1,mu1,sigma1)
y<-rnorm(n2,mu2,sigma2)
z<-c(x,y)
K<-1:length(z)
D<-numeric(m)
for (i in (1:m) ){
  k<-sample(K,size=(n1+n2)/2,replace=FALSE)
  x1<-z[k]
  y1<-z[-k]
  D[i]<-count5test(x1,y1)
}
p<-mean(D)
print(p)

## -----------------------------------------------------------------------------
library(MASS)
library(boot)
library(RANN)
library(energy)
library(Ball)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
A<-matrix(data=c(1,0,0,1),nrow=2)
a<-c(1,1)
B<-matrix(data=c(2,-1,-1,2),nrow=2)
b<-c(1,1)
m <- 1e3; 
k<-3;
p<-2;
mu <- 0.5;
set.seed(12345)
n1 <- n2 <- 50; 
R<-999; 
n <- n1+n2; 
N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- mvrnorm(n1,a,A);
  y <- mvrnorm(n2,b,B);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=R,seed=i*12345)$p.value
}

alpha <- 0.1;
pow <- colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
library(MASS)
library(boot)
library(RANN)
library(energy)
library(Ball)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
A<-matrix(data=c(1,0,0,1),nrow=2)
a<-c(1,1)
B<-matrix(data=c(1.3,0.5,0.5,0.7),nrow=2)
b<-c(0.1,0.5)
m <- 1e3;
k<-3;
p<-2;
set.seed(12345)
n1 <- n2 <- 50; 
R<-999; 
n <- n1+n2; 
N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- mvrnorm(n1,a,A);
  y <- mvrnorm(n2,b,B);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=R,seed=i*12345)$p.value
}

alpha <- 0.1;
pow <- colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
library(boot)
library(RANN)
library(energy)
library(Ball)
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

m <- 1e3; 
k<-3; p<-2; 
mu <- 0.5; set.seed(12345)
n1 <- n2 <- 50; 
R<-999; 
n <- n1+n2; 
N = c(n1,n2)
q<-0.7#Proportional parameter of mixed distribution
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rt(n1*p,df=1),ncol=p);
  y<-numeric(n1*p)
  u<-runif(n1*p)
  l<-which(u<=q)
  y[l]<-rnorm(length(l));
  y[-l]<-rnorm(n1*p-length(l),mean=mu);
  y<-matrix(y,ncol=p)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=R,seed=i*12345)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
library(MASS)
library(boot)
library(RANN)
library(energy)
library(Ball)
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
A<-matrix(data=c(1,0,0,1),nrow=2)
a<-c(1,1)
B<-matrix(data=c(1,0,0,3),nrow=2)
b<-c(1,1)
m <- 1e3; 
k<-3;
p<-2;
mu <- 0.5;
set.seed(12345)
n1 <-20;
n2<-200;
R<-999;
n <- n1+n2; 
N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- mvrnorm(n1,a,A);
  y <- mvrnorm(n2,b,B);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=R,seed=i*12345)$p.value
}

alpha <- 0.1;
pow <- colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
LS<-function(x){
  f<-exp(-abs(x))/2
  return(f)
}

rw.Metropolis<-function(sigma,x0,N){
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for (i in 2:N){
   y<-rnorm(1,x[i-1],sigma)
   if (u[i] <= LS(y) / LS(x[i-1]))
  {x[i] <- y}
else
 {x[i] <- x[i-1]
  k<-k+1
 }
  }
  return(list(x=x,k=k))
}


N<-5000
sigma<-c(0.1,0.5,1,5,20,50)
x0<-10
rw1<-rw.Metropolis(sigma[1],x0,N)
rw2<-rw.Metropolis(sigma[2],x0,N)
rw3<-rw.Metropolis(sigma[3],x0,N)
rw4<-rw.Metropolis(sigma[4],x0,N)
rw5<-rw.Metropolis(sigma[5],x0,N)
rw6<-rw.Metropolis(sigma[6],x0,N)

lk<-(N-c(rw1$k,rw2$k,rw3$k,rw4$k,rw5$k,rw6$k))/N#Acceptance rate
for (l in 1:length(lk)){
  print(paste('The',l,'th chain s corresponding normal standard deviation is:',sigma[l],',Accept rate is :',lk[l]))
}

#par(mfrow = c(3, 2))

plot(rw1$x,type='l',xlab=expression(paste(sigma,'=0.1')))
plot(rw2$x,type='l',xlab=expression(paste(sigma,'=0.5'))) 
plot(rw3$x,type='l',xlab=expression(paste(sigma,'=1')))
plot(rw4$x,type='l',xlab=expression(paste(sigma,'=5')))
plot(rw5$x,type='l',xlab=expression(paste(sigma,'=20')))
plot(rw6$x,type='l',xlab=expression(paste(sigma,'=50')))

mtext(expression(paste("The performance of different ",sigma)), side =3, line = -2.5, outer = TRUE)

## -----------------------------------------------------------------------------
set.seed(3426)
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

LS<-function(x){
  f<-exp(-abs(x))/2
  return(f)
}

rw.Metropolis<-function(sigma,x0,N){
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for (i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    if (u[i] <= LS(y) / LS(x[i-1]))
    {x[i] <- y}
    else
    {x[i] <- x[i-1]
    }
  }
  return(x)
}


RR<-1.2
sigma<-0.5#Consider only 0.5 for a different initial value
Rr<-2
x0 <- c(-10, -5, 5, 10)#The initial value
N<-15000#Let's do it 15000 times, and then let's do the index
rw<-matrix(rep(0,length(x0)*N),nrow=length(x0))

for (rr in 1:length(x0)){
    sig<-x0[rr]
    rw[rr,]<-rw.Metropolis(sigma,sig,N)
}

psi <- t(apply(rw, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))


#plot psi for the six chains
b<-1000
par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)
lines(c(0,N),c(1.2,1.2),type='b')

n<-5500#Look at the graph, let's start at 5500
IN<-numeric(N)
while (Rr>=RR){
  n<-n+1
  pp<-psi[,1:n]
  IN[n]<-Gelman.Rubin(pp)
  Rr<-IN[n]
}

print(paste('When we get to the ',n,',the index we want is less than 1.2'))

print('Here is the mean of each chain:')

par(mfrow=c(1,1))
plot(psi[1,],type='l',col=1,lwd=1.3,xlab="",ylab="",ylim=c((min(psi)-2),(max(psi)+2)))
par(new=TRUE)
plot(psi[2,],type='l',col=2,lwd=1.3,xlab="",ylab="",ylim=c((min(psi)-2),(max(psi)+2)))
par(new=TRUE)
plot(psi[3,],type='l',col=3,lwd=1.3,xlab="",ylab="",ylim=c((min(psi)-2),(max(psi)+2)))
par(new=TRUE)
plot(psi[4,],type='l',col=4,lwd=1.3,xlab="",ylab="",ylim=c((min(psi)-2),(max(psi)+2)))

legend('topright',legend = c('x0=-10','x0=-5', 'x0=5','x0=10'), col = c(1,2,3,4),lty=c(1,2,3,4))
mtext("Mean of each chain", side =3, line = -2.5, outer = TRUE)

## -----------------------------------------------------------------------------
tF<-function(x,n){
   exp(lgamma((n+1)/2)-lgamma(n/2))*(1+x^2/n)^(-(n+1)/2)/(sqrt(n*pi))
}
KKK<-c(4:25,100,500,1000)
#par(mfrow=c(5,length(KK)/5))
for (dd in 1:5){
  KK<-KKK[(5*(dd-1)+1):(5*dd)]
  #par(mfrow=c(1,5))
  for (k in  1:length(KK)){
   kk<-KK[k]
   l<-0.01
   x<-seq(0,sqrt(kk)-l,by=l)
   y<-numeric(length(x))
   for (p in 1:length(x)){
      pp<-x[p]
      l1<-sqrt(pp^2*(kk-1)/(kk-pp^2))
      l2<-sqrt(pp^2*kk/(kk+1-pp^2))
      la<-integrate(tF,lower=l1,upper=Inf,n=(kk-1))
      lb<-integrate(tF,lower=l2,upper=Inf,n=kk)
      y[p]<-la[1]$value-lb[1]$value
   }
   plot(x,y,'l',xlim=c(0,sqrt(kk)),ylim=c(min(y)*1.2,max(y)*1.2),xlab=paste("The corresponding value of k is",kk))
   abline(h=0,lty=3)
  }
}

## -----------------------------------------------------------------------------
tF<-function(x,n){
   exp(lgamma((n+1)/2)-lgamma(n/2))*(1+x^2/n)^(-(n+1)/2)/(sqrt(n*pi))
}
KKK<-c(4:25,100,500,1000) 
tt<-numeric(length(KKK))#Storage root
n<-numeric(length(KKK))
for (k in  1:length(KKK)){
   kk<-KKK[k]
   l<-0.01
   x1<-l#The left endpoint of the interval where the root is
   x2<-sqrt(kk)-l#The right endpoint of the interval where the root is
   barr<-10^(-6)#The threshold value
   y<-numeric(2)
   
   while (abs(x1-x2)>barr){
      tt[k]<-(x1+x2)/2#The midpoint of the interval
      
      al<-sqrt((x1)^2*(kk-1)/(kk-(x1)^2))
      bl<-sqrt((x1)^2*kk/(kk+1-(x1)^2))
      
      cl<-sqrt((tt[k])^2*(kk-1)/(kk-(tt[k])^2))
      dl<-sqrt((tt[k])^2*kk/(kk+1-(tt[k])^2))
      
      ala<-integrate(tF,lower=al,upper=Inf,n=(kk-1))
      blb<-integrate(tF,lower=bl,upper=Inf,n=kk)
      
      cla<-integrate(tF,lower=cl,upper=Inf,n=(kk-1))
      dlb<-integrate(tF,lower=dl,upper=Inf,n=kk)
      y[1]<-ala[1]$value-blb[1]$value
      y[2]<-cla[1]$value-dlb[1]$value
      if(y[1]*y[2]<=0){
         x2<-tt[k]
      }
      else{
         x1<-tt[k] 
      }
      n[k]<-n[k]+1
   }     
}

tty<-numeric(length(tt))


for (dd in 1:5){
   l<-5*(dd-1)+1
   KK<-KKK[l:(l+4)]
   #par(mfrow=c(1,5))
   for (k in  1:length(KK)){
      kk<-KK[k]
      ll<-0.01
      xx<-seq(0,sqrt(kk)-ll,by=ll)
      yy<-numeric(length(xx))
      for (p in 1:length(xx)){
         pp<-xx[p]
         l1<-sqrt(pp^2*(kk-1)/(kk-pp^2))
         l2<-sqrt(pp^2*kk/(kk+1-pp^2))
         la<-integrate(tF,lower=l1,upper=Inf,n=(kk-1))
         lb<-integrate(tF,lower=l2,upper=Inf,n=kk)
         yy[p]<-la[1]$value-lb[1]$value
      }
      
      tt1<-sqrt(tt[l+k-1]^2*(kk-1)/(kk-tt[l+k-1]^2))
      tt2<-sqrt(tt[l+k-1]^2*kk/(kk+1-tt[l+k-1]^2))
      tta<-integrate(tF,lower=tt1,upper=Inf,n=(kk-1))
      ttb<-integrate(tF,lower=tt2,upper=Inf,n=kk)
      tty[l+k-1]<-tta[1]$value-ttb[1]$value
      
      plot(xx,yy,'l',xlab="",xlim=c(0,sqrt(kk)),ylim=c(min(yy)*1.2,max(yy)*1.2))
      abline(h=0,lty=3)
      par(new=TRUE)
      plot(tt[l+k-1],tty[l+k-1],'b',col='red',xlab=paste("The corresponding k value is:",kk),ylab="",xlim=c(0,sqrt(kk)),ylim=c(min(yy)*1.2,max(yy)*1.2))
   }
    if(dd==1){
    mtext("The location of each root (red circle)", side =3, line = -2.5, outer = TRUE)
     }
}

for (ind in 1:length(tt)){
   print(paste('when k=',KKK[ind],',the corresponding root is:',tt[ind]))
}

## -----------------------------------------------------------------------------
n<-15
L<-matrix(rep(0,n*3),nrow=n)#Three parameters
R<-numeric(n)#Maximum likelihood for observed data
L[1,1]<-L[1,2]<-0.1
L[1,3]<-0.8
bar<-1e-8
k=1
Ba<-sqrt((L[k,1])^2+(L[k,1])^2+(L[k,1])^2)
while(Ba>bar){
  p<-L[k,1]
  q<-L[k,2]
  r<-L[k,3]
  R[k]<-444*log(p^2+2*p*r)+132*log(q^2+2*q*r)+722*log(r)+63*log(2*p*q)
  M<-444+63+444*p^2/(p^2+2*p*r)
  N<-132+63+132*q^2/(q^2+2*q*r)
  S<-2*361+444+132-132*q^2/(q^2+2*q*r)-444*p^2/(p^2+2*p*r)
  L[k+1,1]<-M/(M+N+S)
  L[k+1,2]<-N/(M+N+S)
  L[k+1,3]<-S/(M+N+S)
  Ba<-sqrt((L[k+1,1]-L[k,1])^2+(L[k+1,2]-L[k,2])^2+(L[k+1,3]-L[k,3])^2)
  k<-k+1
}

print("The three parameter values (p,q,r) estimated each time are respectively:")
print(L[1:(k-1),])
print("The logarithmic likelihood function value of the observed data is as follows:")

plot(R[1:(k-1)],type='b',col='red',xlab="Number of calculations",ylab="The logarithmic likelihood function value of the observed data")

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

# lapply
L1<- lapply(formulas, lm, data = mtcars)
L2<- lapply(formulas, function(x) lm(formula = x, data = mtcars))

# circulation
L3 <- vector("list", length(formulas))
n<-seq_along(formulas)
for (i in n){
  L3[[i]] <- lm(formulas[[i]], data = mtcars)
}

print("The results of lapply:")
print(L1)
print("The results of loop:")
print(L3)

## -----------------------------------------------------------------------------
trials <- replicate(
  100, 
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

# Using sapply () and anonymous functions:
sapply(trials, function(x) x[["p.value"]])
# Do not use anonymous functions:
sapply(trials, "[[", "p.value")

## -----------------------------------------------------------------------------
testlist <- list(iris, mtcars, cars)

lapply(testlist, function(x) vapply(x, mean, numeric(1)))

lmapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE){return(simplify2array(out))}
  out
}

lmapply(testlist, mean, numeric(1))

## -----------------------------------------------------------------------------
library(Rcpp)
dir_cpp <- '../src/'
#Can create source file in Rstudio
sourceCpp(paste0(dir_cpp,"Metropolis.cpp"))
#library(StatComp20091)
library(microbenchmark)

LSR<-function(x){
  f<-exp(-abs(x))/2
  return(f)
}

rw.Metropolis<-function(sigma,x0,N){
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for (i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    if (u[i] <= LSR(y) / LSR(x[i-1]))
    {x[i] <- y}
    else
    {x[i] <- x[i-1]
    k<-k+1
    }
  }
  return(list(x=x,k=k))
}


N<-5000
sigma<-c(0.5,1,5,10)
x0<-10
lk<-matrix(rep(0,length(sigma)*2),ncol=2)

for ( i in c(1,2,3,4)){
  ts<-microbenchmark(rw1<-rw.Metropolis(sigma[i],x0,N),rw2<-Metropolis(sigma[i],x0,N))
 
  

  #par(mfrow = c(1, 2))
  plot(rw1$x,type='l',col='red',xlab=paste('R:',expression(sigma),'=',sigma[i])) 
  plot(rw2[,1],type='l',col='blue',xlab=paste('C++:',expression(sigma),'=',sigma[i]))
  
  mtext("The performance of different language ", side =3, line = -2.5, outer = TRUE)
   print(paste('when sigma=',sigma[i],',the running time of the corresponding two cases is:'))
  RS<-summary(ts)[,c(1,3,5,6)]
  print(RS)
  lk[i,]<-(c(N-(rw1$k),N-rw2[1,3]))/N#Accept rate

  print(paste('The normal standard deviation of the first chain in R is:',sigma[i],',accept rate is :',lk[i,1],';'))
  print(paste('The normal standard deviation of the first chain in C++:',sigma[i],',accept rate is :',lk[i,2]))
  
}

