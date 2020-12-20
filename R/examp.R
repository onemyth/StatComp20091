#' @title Bobble
#' @param L The matrix with n rows and 1 column for sorting
#' @param n  The rows of the matrix
#' @return The index of results
#' @examples
#' \dontrun{
#' L<-matrix(rnorm(20,1,4),nrow=20)
#' n<-20
#' LL<-numeric(n)
#' Ind<-Bobble(L,n)
#' Res<-numeric(n)
#' for ( k in(1:n)){
#'   ids<-Ind[k]
#'   Res[k]<-L[ids]
#'   LL[k]<-L[k,1]
#' }
#' print("the original order is:\n")
#' print(LL)
#' print("the sorted order is:\n")
#' print(Res)
#' }
#' @import Rcpp
#' @useDynLib StatComp20091
#' @export
Bobble<-function(L,n){
  requireNamespace('Rcpp', quietly = TRUE)
  cppFunction('NumericMatrix Bobble(NumericMatrix x,int m){
  NumericMatrix mat(m, 1);
  NumericMatrix Index(m,1);
  for (int k=0;k<m;k++){
    mat(k,0)=x(k,0);
    Index(k,0)=k+1;
  }
  int i,j,t;
  double l;
  for(i=0;i<(m-1);i++){
     for(j=0;j<(m-1-i);j++){
        if(mat(j,0)>mat((j+1),0)){
        l=mat(j,0);mat(j,0)=mat((j+1),0);mat((j+1),0)=l;
        t=Index(j,0);Index(j,0)=Index((j+1),0);Index((j+1),0)=t;
        }}}
  return(Index);
}')
LL<-numeric(n)
Ind<-Bobble(L,n)
return(Ind)
}


#' @title  compute factorials under recursion
#' @param n for n in n!
#' @return n the value of n!
#' @examples
#' \dontrun{
#' y<-fac(5)
#' }
#' @import Rcpp
#' @useDynLib StatComp20091
#' @export
fac<-function(n){
requireNamespace('Rcpp', quietly = TRUE)
cppFunction('int fac(int m){
  int P;
  if((m==0)||(m==1)){
    P=1;
  }
  else {P=m*fac(m-1);}
  return(P);
}')
y<-fac(n)
return(y)
}


#' @title Multiplicative kernel
#' @description  calculate the distribution based on the Multiplicative kernel
#' @param a1 with b1 to determine the bar of x
#' @param a2 with b2 to determine the bar of y
#' @param b1 with a1 to determine the bar of x
#' @param b2 with a2 to determine the bar of y
#' @param ha the bandwidth of x
#' @param hb the bandwidth of y
#'
#' @return the result of Multiplicative kernel
#' @examples
#' \dontrun{
#' library(ks)
#' library(rgl)
#' xx <- as.matrix(faithful)
#' kdex <- kde(xx)
#' data("faithful")
#' ef<-faithful$eruptions
#' wf<-faithful$waiting
#' h1<-0.3
#' h2<-6
#' x<-unlist(kdex$eval.points[1])
#' y<-unlist(kdex$eval.points[2])
#' n<-length(ef)
#' n1<-length(x)
#' n2<-length(y)
#' z1<-matrix(rep(0,n1*n2),nrow=n1)
#' for ( i in 1:n1){
#'   for (j in 1:n2){
#'     k1<-0
#'     for (l in 1:n){
#'       k1<-k1+fcount(ef[l],wf[l],x[i],y[j],h1,h2)
#'     }
#'     z1[i,j]<-k1/n
#'   }
#' }
#' persp(x, y, z1,phi = 30, theta = 40, col = 'black', border = 0
#'       ,xlab='eruptions',ylab='waiting'
#'       ,main='Results under Multiplicative kernel')
#'  }
#' @export
fcount<-function(a1,a2,b1,b2,ha,hb){#Multiplicative kernel
  s<-0
  bar1<-(b1-a1)/ha
  bar2<-(b2-a2)/hb
  bar<-((-1<=bar1)&(bar1<1)&(-1<=bar2)&(bar2<1))
  if (bar){s<-1}
  return(s/(ha*hb))
}


#' @title the multivariate Epanechnikov (spherical) 
#' @description  calculate the distribution based on the multivariate Epanechnikov  kernel
#' @param a1 with b1 to determine the bar of x
#' @param a2 with b2 to determine the bar of y
#' @param b1 with a1 to determine the bar of x
#' @param b2 with a2 to determine the bar of y
#' @param ha the bandwidth of x
#' @param hb the bandwidth of y
#' @return the result of  the multivariate Epanechnikov (spherical) kernel
#' @examples
#' \dontrun{
#' library(ks)
#' library(rgl)
#' xx <- as.matrix(faithful)
#' kdex <- kde(xx)
#' data("faithful")
#' ef<-faithful$eruptions
#' wf<-faithful$waiting
#' h1<-0.3
#' h2<-6
#' x<-unlist(kdex$eval.points[1])
#' y<-unlist(kdex$eval.points[2])
#' n<-length(ef)
#' n1<-length(x)
#' n2<-length(y)
#' z2<-matrix(rep(0,n1*n2),nrow=n1)
#' for ( i in 1:n1){
#'   for (j in 1:n2){
#'     k2<-0
#'     for (l in 1:n){
#'       k2<-k2+fEcount(ef[l],wf[l],x[i],y[j],h1,h2)
#'     }
#'     z2[i,j]<-k2/n
#'   }
#' }
#' persp(x, y, z2,phi = 30, theta = 40, col = 'red', border = 0,
#'      xlab='eruptions',ylab='waiting'
#'       ,main='Results under  Multivariate Epanechnikov (spherical)')
#'  }
#' @export
fEcount<-function(a1,a2,b1,b2,ha,hb){#The multivariate Epanechnikov (spherical)
  s<-0
  bar1<-(b1-a1)/ha
  bar2<-(b2-a2)/hb
  bar3<-c(bar1,bar2)
  bar<-t(bar3)%*%(bar3)
  if(bar<=1){
    s<-(1-bar)
  }
  return(s/(ha*hb))
}



#' @title The multivariate Epanechnikov (multiplicative)
#' @description  calculate the distribution based on the multivariate Epanechnikov kernel 
#' @param a1 with b1 to determine the bar of x
#' @param a2 with b2 to determine the bar of y
#' @param b1 with a1 to determine the bar of x
#' @param b2 with a2 to determine the bar of y
#' @param ha the bandwidth of x
#' @param hb the bandwidth of y
#' @return the result of he multivariate Epanechnikov (multiplicative) kernel
#' @examples
#' \dontrun{
#' library(ks)
#' library(rgl)
#' xx <- as.matrix(faithful)
#' kdex <- kde(xx)
#' data("faithful")
#' ef<-faithful$eruptions
#' wf<-faithful$waiting
#' h1<-0.3
#' h2<-6
#' x<-unlist(kdex$eval.points[1])
#' y<-unlist(kdex$eval.points[2])
#' n<-length(ef)
#' n1<-length(x)
#' n2<-length(y)
#' z3<-matrix(rep(0,n1*n2),nrow=n1)
#' for ( i in 1:n1){
#'   for (j in 1:n2){
#'     k3<-0
#'     for (l in 1:n){
#'       k3<-k3+fE2count(ef[l],wf[l],x[i],y[j],h1,h2)
#'     }
#'    z3[i,j]<-k3/n
#'   }
#' }
#' persp(x, y, z3,phi = 30, theta = 40, col = 'blue', border = 0
#'       ,xlab='eruptions',ylab='waiting'
#'       ,main='Results under  Multivariate Epanechnikov (multiplicative)')
#'  }
#' @export
fE2count<-function(a1,a2,b1,b2,ha,hb){#The multivariate Epanechnikov (multiplicative)
  s<-0
  bar1<-(b1-a1)/ha
  bar2<-(b2-a2)/hb
  bar<-((abs(bar1)<=1)&(abs(bar2)<=1))
  if (bar){
    s<-(3/4)^2*(1-(bar1)^2)*(1-(bar2)^2)
  }
  return(s/(ha*hb))
}