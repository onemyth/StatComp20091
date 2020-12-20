#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs  sampler using Rcpp
//' @param sigma the sigma value
//' @param x0 the start point
//' @param N the number of circles
//' @return a matrix of results
//' @export
// [[Rcpp::export]]
NumericMatrix Metropolis(double sigma ,float x0,int N) {
  
  NumericMatrix mat(N, 3);
  mat(0,0)=x0;
  int k=0;
 
  for (int i=2;i<N;i++){
  	double u=::Rf_runif(0,1); 
  	int oo=i-1;
  	double y=::Rf_rnorm(mat(oo,0),sigma);
  	mat(i-1,1)=y;
  	if(u<=(exp(fabs(mat(oo,0))-fabs(y)))){
  		mat(i,0)=y;
  	}
  	else{
  		mat(i,0)=mat(oo,0);
  		k=k+1;
  	}
  }
  mat(0,2)=k;
  return(mat);
}
