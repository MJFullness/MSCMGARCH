
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "Copulas.h"
// [[Rcpp::plugins(cpp11)]]



// // [[Rcpp::export]]
// SEXP BiCop_cpp(int family, double par ,double par2){
//   
//   // Obtaining namespace of Matrix package
//   Rcpp::Environment pkg = Rcpp::Environment::namespace_env("VineCopula");
//   
//   
//   // Picking up Matrix() function from Matrix package
//   Rcpp::Function f = pkg["BiCop"];
//   SEXP res = f(family, par, par2);
//   // Executing Matrix( m, sparse = TRIE )
//   return res;
// }
// 
// // [[Rcpp::export]]
// double BiCopHfunc2_cpp(double returns1, double returns2 , SEXP family){
//   
//   // Obtaining namespace of Matrix package
//   Rcpp::Environment pkg = Rcpp::Environment::namespace_env("VineCopula");
//   
//   
//   // Picking up Matrix() function from Matrix package
//   Rcpp::Function f = pkg["BiCopHfunc2"];
//   Rcpp::NumericVector res = f(returns1, returns2, family);
//   // Executing Matrix( m, sparse = TRIE )
//   return Rcpp::as<double>(res);
// }
// 
// // [[Rcpp::export]]
// double BiCopHfunc1_cpp(double returns1, double returns2 , SEXP family){
//   
//   // Obtaining namespace of Matrix package
//   Rcpp::Environment pkg = Rcpp::Environment::namespace_env("VineCopula");
//   
//   
//   // Picking up Matrix() function from Matrix package
//   Rcpp::Function f = pkg["BiCopHfunc1"];
//   Rcpp::NumericVector res = f(returns1, returns2, family);
//   // Executing Matrix( m, sparse = TRIE )
//   return Rcpp::as<double>(res);
// }



// [[Rcpp::export]]
double loglike_Normal_Gumbel_GumbelSurvival(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  double par_gumbel =theta[6];
  double par_gumbelSurvival =theta[7];
  //check for identification
  if(valid_bekk(C,A,G)==false  || p13 < 0 || p13 > 1  || p31 < 0 || p31 > 1  || p23 < 0 || p23 > 1  || p11+p12>1 || p21+p22>1 || p32+p33>1 ||p11<0 || p12<0 || p11>1 || p12>1 || p21<0 || p22<0 || p21>1 || p22>1  || p33<0 || p32<0 || p33>1 || p32>1||par_gumbelSurvival<=1 || par_gumbelSurvival >=100 || par_gumbel<= 1|| par_gumbel>=100){
    return -1e25;
  }
  
 
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",0,par_gumbel));
  //  
  arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double gumbel_det = arma::det(cor_gumbel_chol);
  
  //set_seed(3592);
  //arma::mat cor_gumbelSurvival=arma::cor(rbicop_cpp(100000,"gumbel",180,par_gumbelSurvival));
  arma::mat cor_gumbelSurvival=cor_Gumbel(par_gumbelSurvival);
  arma::mat cor_gumbelSurvival_chol=arma::chol(cor_gumbelSurvival).t();
  
  double gumbelSurvival_det = arma::det(cor_gumbelSurvival_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double H_det=arma::det(H_eigen_inv);
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  
  
  
  double llv=log(arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*H_det);;
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(i).t());
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*H_det);
    double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*gumbelPDF(p_norm_1(0),p_norm_1(1),par_gumbel)*H_det* gumbel_det);
    double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*gumbelPDF(1-p_norm_2(0),1-p_norm_2(1),par_gumbelSurvival)*H_det*gumbelSurvival_det);
    //double llv_temp2 = (p1t*p12+p2t*p22+p3t*p32)*d_norm_1*dbicop_cpp(p_norm_1,"gumbel",0,par_gumbel)*H_det* gumbel_det;
    //double llv_temp3 = (p1t*p13+p2t*p23+p3t*p33)*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",180,par_gumbelSurvival)*H_det*gumbelSurvival_det;
    
    double llv_temp = llv_temp1+llv_temp2+llv_temp3;
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = 1-p1t-p2t;
    llv+=log(llv_temp);
    
    
  }
  
  return llv;
}

// [[Rcpp::export]]
double loglike_Normal_Gumbel_Clayton(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  double par_gumbel =theta[6];
  double par_clayton =theta[7];
  //check for identification
  if(valid_bekk(C,A,G)==false  || p13 < 0 || p13 > 1  || p31 < 0 || p31 > 1  || p23 < 0 || p23 > 1  || p11+p12>1 || p21+p22>1 || p32+p33>1 ||p11<0 || p12<0 || p11>1 || p12>1 || p21<0 || p22<0 || p21>1 || p22>1  || p33<0 || p32<0 || p33>1 || p32>1||par_clayton<=0 || par_clayton >=100 || par_gumbel<= 1|| par_gumbel>=100){
    return -1e25;
  }
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",0,par_gumbel));
  
  arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double gumbel_det = arma::det(cor_gumbel_chol);
  
  // set_seed(3592);
  // arma::mat cor_clayton=arma::cor(rbicop_cpp(100000,"clayton",0,par_clayton));
  arma::mat cor_clayton=cor_Clayton(par_clayton);
  arma::mat cor_clayton_chol=arma::chol(cor_clayton).t();
  
  double clayton_det = arma::det(cor_clayton_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double H_det=arma::det(H_eigen_inv);
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  
  
  
  double llv=log(arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*H_det);;
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t());
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*H_det);
    //double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*dbicop_cpp(p_norm_1,"gumbel",0,par_gumbel)*H_det* gumbel_det);
    //double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*dbicop_cpp(p_norm_2,"clayton",0,par_clayton)*H_det*clayton_det);
    double llv_temp2 = (p1t*p12+p2t*p22+p3t*p32)*d_norm_1*gumbelPDF(p_norm_1(0),p_norm_1(1),par_gumbel)*H_det* gumbel_det;
    double llv_temp3 = (p1t*p13+p2t*p23+p3t*p33)*d_norm_2*claytonPDF(p_norm_2(0),p_norm_2(1),par_clayton)*H_det*clayton_det;
    
    double llv_temp = llv_temp1+llv_temp2+llv_temp3;
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = 1-p1t-p2t;
    llv+=log(llv_temp);
    
    
  }
  
  return llv;
}

// [[Rcpp::export]]
double loglike_Normal_GumbelS_ClaytonS(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  double par_gumbel =theta[6];
  double par_clayton =theta[7];
  //check for identification
  if(valid_bekk(C,A,G)==false  || p13 < 0 || p13 > 1  || p31 < 0 || p31 > 1  || p23 < 0 || p23 > 1  || p11+p12>1 || p21+p22>1 || p32+p33>1 ||p11<0 || p12<0 || p11>1 || p12>1 || p21<0 || p22<0 || p21>1 || p22>1  || p33<0 || p32<0 || p33>1 || p32>1||par_clayton<=0 || par_clayton >=100 || par_gumbel<= 1|| par_gumbel>=100){
    return -1e25;
  }
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  //set_seed(3592);
  //arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",180,par_gumbel));
     
  arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double gumbel_det = arma::det(cor_gumbel_chol);
  
  // set_seed(3592);
  // arma::mat cor_clayton=arma::cor(rbicop_cpp(100000,"clayton",180,par_clayton));
  arma::mat cor_clayton=cor_Clayton(par_clayton);
  arma::mat cor_clayton_chol=arma::chol(cor_clayton).t();
  
  double clayton_det = arma::det(cor_clayton_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double H_det=arma::det(H_eigen_inv);
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  
  
  
  double llv=log(arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*H_det);;
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t());
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*H_det);
    //double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*dbicop_cpp(p_norm_1,"gumbel",180,par_gumbel)*H_det* gumbel_det);
    //double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*dbicop_cpp(p_norm_2,"clayton",180,par_clayton)*H_det*clayton_det);
    double llv_temp2 = (p1t*p12+p2t*p22+p3t*p32)*d_norm_1*gumbelPDF(1-p_norm_1(0),1-p_norm_1(1),par_gumbel)*H_det* gumbel_det;
    double llv_temp3 = (p1t*p13+p2t*p23+p3t*p33)*d_norm_2*claytonPDF(1-p_norm_2(0),1-p_norm_2(1),par_clayton)*H_det*clayton_det;
    
    double llv_temp = llv_temp1+llv_temp2+llv_temp3;
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = 1-p1t-p2t;
    llv+=log(llv_temp);
    
    
  }
  
  return llv;
}


// [[Rcpp::export]]
double loglike_Normal_Clayton_ClaytonSurvival(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  double par_clayton =theta[6];
  double par_claytonSurvival =theta[7];
  //check for identification
  if(valid_bekk(C,A,G)==false  || p13 >= 1||   p13 <= 0|| p23 >= 1|| p21+p22>=1 || p32+p33>=1 || p11+p12>=1 || p23 <= 0||  p31 >= 1||  p31 <= 0 || p11<=0 || p12<=0 || p11>=1 || p12>=1 || p21<=0 || p22<=0 || p21>=1 || p22>=1  || p33<=0 || p32<=0 || p33>=1 || p32>=1||par_claytonSurvival<=0 || par_claytonSurvival >=100 || par_clayton<= 0|| par_clayton>=100){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_clay = arma::cor(rbicop_cpp(100000,"clayton",0,par_clayton));
  arma::mat cor_clay=cor_Clayton(par_clayton);
  arma::mat cor_clay_chol=arma::chol(cor_clay).t();
  
  double clay_det = arma::det(cor_clay_chol);
  
  arma::mat cor_claytonSurvival=cor_Clayton(par_claytonSurvival);
  // set_seed(3592);
  // arma::mat cor_claytonSurvival = arma::cor(rbicop_cpp(100000,"clayton",180,par_claytonSurvival));
  arma::mat cor_claytonSurvival_chol=arma::chol(cor_claytonSurvival).t();
  
  double claytonSurvival_det = arma::det(cor_claytonSurvival_chol);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  double H_det=arma::det(H_eigen_inv);
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;

  
  double llv = log(arma::as_scalar(arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t())))*H_det);
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_claytonSurvival_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_claytonSurvival_chol*H_eigen_inv*r.row(i).t());
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*H_det );
    double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*claytonPDF(p_norm_1(0),p_norm_1(1),par_clayton)*H_det* clay_det);
    double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*claytonPDF(1-p_norm_2(0),1-p_norm_2(1),par_claytonSurvival)*H_det*claytonSurvival_det);
    // double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*dbicop_cpp(p_norm_1,"clayton",0,par_clayton)*H_det* clay_det);
    // double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*dbicop_cpp(p_norm_2,"clayton",180,par_claytonSurvival)*H_det*claytonSurvival_det);
     
    double llv_temp = llv_temp1+llv_temp2+llv_temp3;
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = 1-p1t-p2t;
    llv+=log(llv_temp);
    
    
  }
  
  return llv;
}


// [[Rcpp::export]]
double loglike_Clayton_Gumbel(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_clayton =theta[2];
  double par_gumbel =theta[3];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1 || par_clayton<=0 || par_clayton >=17 || par_gumbel<= 1|| par_gumbel>=17){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  set_seed(3592);
  arma::mat cor_clayton=arma::cor(BiCopSim_cpp(100000,3,par_clayton));
  arma::mat cor_clayton_chol=arma::chol(cor_clayton).t();
  set_seed(3592);
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,4,par_gumbel));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
 
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));

  
//initial state
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*arma::det(cor_clayton_chol)*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),3,par_clayton)));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
      
    arma::vec p_norm_1=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
      
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*arma::det(H_eigen_inv)*arma::det(cor_clayton_chol)*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),3,par_clayton)+(p1t*(1-p)+(1-p1t)*(q))*arma::det(H_eigen_inv)*arma::det(cor_gumbel_chol)*d_norm_2*BiCopPDF_cpp(p_norm_2(0),p_norm_2(1),4,par_gumbel));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*arma::det(H_eigen_inv)*arma::det(cor_clayton_chol)*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),3,par_clayton))/llv_temp;
      
    llv+=log(llv_temp);
     
    
      }
  return llv;
}


// [[Rcpp::export]]
double loglike_Clayton_Gumbel90(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_clayton =theta[2];
  double par_gumbel =theta[3];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1 || par_clayton<=0 || par_clayton >=17 || par_gumbel<= -17|| par_gumbel>=-1){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  set_seed(3592);
  arma::mat cor_clayton=arma::cor(BiCopSim_cpp(100000,3,par_clayton));
  arma::mat cor_clayton_chol=arma::chol(cor_clayton).t();
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,24,par_gumbel));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double clay_det=arma::det(cor_clayton_chol);
  double gumb_det=arma::det(cor_gumbel_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*clay_det*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),3,par_clayton)));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double H_det =arma::det(H_eigen_inv);
    
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    double llv_1 =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*clay_det*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),3,par_clayton));
    double llv_2 =arma::as_scalar((p1t*(1-p)+(1-p1t)*(q))*H_det*gumb_det*d_norm_2*BiCopPDF_cpp(p_norm_2(0),p_norm_2(1),24,par_gumbel));
    double llv_temp = llv_1+llv_2;
    p1t=llv_1/llv_temp;
    llv+=log(llv_temp);
    
    
  }
  return llv;
}

// [[Rcpp::export]]
double loglike_Normal_Gumbel(const arma::vec& bekk, const arma::vec& theta, const  arma::mat& r) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }

  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_gumbel =theta[2];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_gumbel<= 1|| par_gumbel>=100){
    return -1e25;
  }

  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;

  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();


  //compute covariance matrix

  
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",0,par_gumbel));
  arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  double gumb_det=arma::det(cor_gumbel_chol);

  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));



  double p1t = 1;

  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
  //double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));

  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());

  double llv = log(arma::as_scalar(d_norm_1));


  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;

    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double h_det= arma::det(H_eigen_inv);


    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*h_det;;

    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));

    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    double llv_temp =(p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*gumbelPDF(p_norm_2(0),p_norm_2(1),par_gumbel);
    
    //double llv_temp =(p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",0,par_gumbel);
    // double llv_temp =(p1t*p+(1-p1t)*(1-q))*d_norm_1+std::exp(std::log(p1t*(1-p)+(1-p1t)*(q))+std::log(h_det)+std::log(gumb_det)+std::log(d_norm_2)+gumbelPDF_log(p_norm_2(0),p_norm_2(1),par_gumbel));
    // double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",0,par_gumbel));

    p1t=(p1t*p+(1-p1t)*(1-q))*d_norm_1/llv_temp;

    llv+=log(llv_temp);


  }
  
  return llv;
}

// // [[Rcpp::export]]
// double loglike_Normal_Gumbel1(const arma::vec& bekk, const arma::vec& theta, const  arma::mat& r) {
//   // Log-Likelihood function
//   
//   // convert to matrices
//   int n = r.n_cols;
//   // Length of each series
//   int NoOBs = r.n_rows;
//   int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
//   arma::mat C = arma::zeros(n, n);
//   
//   int index = 0;
//   
//   for(int i = 0; i < n; i++){
//     for (int j = i; j < n; j++) {
//       C(j, i) = bekk(index);
//       index += 1;
//     }
//   }
//   
//   arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
//   arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
//   double p = theta[0];
//   double q = theta[1];
//   double par_gumbel =theta[2];
//   //check for identification
//   if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_gumbel<= 1|| par_gumbel>=17){
//     return -1e25;
//   }
//   
//   // compute inital H
//   arma::mat H = (r.t() * r) / r.n_rows;
//   
//   arma::mat CC  = C * C.t();
//   arma::mat At  = A.t();
//   arma::mat Gt  = G.t();
//   
//   
//   //compute covariance matrix
//   
//   
//   
//   arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
//   arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
//   double gumb_det=arma::det(cor_gumbel_chol);
//   
//   arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
//   
//   
//   
//   double p1t = 1;
//   
//   double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
//   //double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
//   
//   arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
//   arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
//   
//   double llv = log(arma::as_scalar(d_norm_1));
//   
//   
//   for (int i = 1; i < NoOBs; i++) {
//     H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
//     
//     arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
//     double h_det= arma::det(H_eigen_inv);
//     
//     
//     double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*h_det;;
//     
//     double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
//     
//     arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
//     
//     double llv_temp1 =(p1t*p+(1-p1t)*(1-q))*d_norm_1;
//     double llv_temp2 =(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",0,par_gumbel);
//     
//     double llv_temp =llv_temp1+llv_temp2;
//     
//     // double llv_temp =(p1t*p+(1-p1t)*(1-q))*d_norm_1+std::exp(std::log(p1t*(1-p)+(1-p1t)*(q))+std::log(h_det)+std::log(gumb_det)+std::log(d_norm_2)+gumbelPDF_log(p_norm_2(0),p_norm_2(1),par_gumbel));
//     // double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",0,par_gumbel));
//     
//     p1t=(p1t*p+(1-p1t)*(1-q))*d_norm_1/llv_temp;
//    
//     llv+=llv_temp;
//     
//     
//   }
//   
//   return llv;
// }



// [[Rcpp::export]]
double loglike_Normal_GumbelS(const arma::vec& bekk, const arma::vec& theta, const  arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_gumbel =theta[2];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_gumbel<= 1|| par_gumbel>=100){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  
 
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",180,par_gumbel));
  arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  double gumb_det=arma::det(cor_gumbel_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
  //double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double h_det= arma::det(H_eigen_inv);
    
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*arma::det(H_eigen_inv);;
    
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    //double llv_temp =(p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",180,par_gumbel);
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*gumbelPDF(1-p_norm_2(0),1-p_norm_2(1),par_gumbel));
    
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1)/llv_temp;
    
    llv+=log(llv_temp);
    
    
  }
  
  return llv;
}

// [[Rcpp::export]]
double loglike_Normal_Gumbel1(const arma::vec& bekk, const arma::vec& theta, const  arma::mat& r) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }

  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_gumbel =theta[2];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_gumbel<= 1|| par_gumbel>=16){
    return -1e25;
  }

  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;

  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();


  //compute covariance matrix


  //arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
  set_seed(3592);
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,4,par_gumbel));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  double gumb_det=arma::det(cor_gumbel_chol);

  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));



  double p1t = 1;

  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
  //double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));

  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());

  double llv = log(arma::as_scalar(d_norm_1)*arma::det(H_eigen_inv));


  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;

    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double h_det= arma::det(H_eigen_inv);


    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*arma::det(H_eigen_inv);;

    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));

    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());

    //double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*gumbelPDF(p_norm_2(0),p_norm_2(1),par_gumbel));
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*BiCopPDF_cpp(p_norm_2(0),p_norm_2(1),4,par_gumbel));

    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1)/llv_temp;

    llv+=log(llv_temp);


  }
  return llv;
}


// // [[Rcpp::export]]
// double loglike_Normal_Gumbel(const arma::vec& bekk, const arma::vec& theta, const  arma::mat& r) {
//   // Log-Likelihood function
//   
//   // convert to matrices
//   int n = r.n_cols;
//   // Length of each series
//   int NoOBs = r.n_rows;
//   int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
//   arma::mat C = arma::zeros(n, n);
//   
//   int index = 0;
//   
//   for(int i = 0; i < n; i++){
//     for (int j = i; j < n; j++) {
//       C(j, i) = bekk(index);
//       index += 1;
//     }
//   }
//   
//   arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
//   arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
//   double p = theta[0];
//   double q = theta[1];
//   double par_gumbel =theta[2];
//   //check for identification
//   if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_gumbel<= 1|| par_gumbel>=50){
//     return -1e25;
//   }
//   
//   // compute inital H
//   arma::mat H = (r.t() * r) / r.n_rows;
//   
//   arma::mat CC  = C * C.t();
//   arma::mat At  = A.t();
//   arma::mat Gt  = G.t();
//   
//   
//   //compute covariance matrix
//   
//   set_seed(3592);
//   arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",0,par_gumbel));
//   //arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",0,par_gumbel));
//   arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
//   double gumb_det=arma::det(cor_gumbel_chol);
//   
//   arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
//   
//   
//   
//   double p1t = 1;
//   
//   double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
//   //double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
//   
//   arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
//   arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
//   
//   double llv = log(arma::as_scalar(d_norm_1));
//   
//   
//   for (int i = 1; i < NoOBs; i++) {
//     H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
//     
//     arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
//     double h_det= arma::det(H_eigen_inv);
//     
//     
//     double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*arma::det(H_eigen_inv);;
//     
//     double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
//     
//     arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
//     
//     double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",0,par_gumbel));
//     //double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*h_det*gumb_det*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",0,par_gumbel));
//     
//     p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1)/llv_temp;
//     
//     llv+=log(llv_temp);
//     
//     
//   }
//   return llv;
// }


// [[Rcpp::export]]
double loglike_Normal_GumbelS1(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }

  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_gumbel =theta[2];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_gumbel<= 1|| par_gumbel>=100){
    return -1e25;
  }

  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;

  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();


  //compute covariance matrix

  set_seed(3592);
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,4,par_gumbel));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();


  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));



  double p1t = 1;

  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;

  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));

  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());

  double llv = log(arma::as_scalar(d_norm_1));


  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;

    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));



    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*arma::det(H_eigen_inv);;
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));

    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());

    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*arma::det(H_eigen_inv)*arma::det(cor_gumbel_chol)*d_norm_2*BiCopPDF_cpp(1-p_norm_2(0),1-p_norm_2(1),4,par_gumbel));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1)/llv_temp;

    llv+=log(llv_temp);


  }
  return llv;
}

// [[Rcpp::export]]
double loglike_Normal_Frank(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_frank =theta[2];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_frank<= -35 || par_frank>=35 || par_frank==0 ){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  
  set_seed(3592);
  arma::mat cor_frank=arma::cor(rbicop_cpp(100000,"frank",0,par_frank));
  arma::mat cor_frank_chol=arma::chol(cor_frank).t();
  
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
  
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_frank_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*arma::det(H_eigen_inv);;
    double d_norm_2=arma::prod(dnorm_cpp(cor_frank_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_frank_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*arma::det(H_eigen_inv)*arma::det(cor_frank_chol)*d_norm_2*dbicop_cpp(p_norm_2,"frank",0,par_frank));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1)/llv_temp;
    
    llv+=log(llv_temp);
    
    
  }
  return llv;
}


// [[Rcpp::export]]
double loglike_Normal_Clayton1(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }

  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_clayton =theta[2];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_clayton<= 0|| par_clayton>=100){
    return -1e25;
  }

  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;

  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();


  //compute covariance matrix

  set_seed(3592);
  // arma::mat cor_clay=arma::cor(BiCopSim_cpp(100000,13,par_clayton));
  arma::mat cor_clay=cor_Clayton(par_clayton);
  arma::mat cor_clay_chol=arma::chol(cor_clay).t();


  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));



  double p1t = 1;

  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));

  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(0).t());

  double llv = log(arma::as_scalar(d_norm_1));


  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;

    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));



    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*arma::det(H_eigen_inv);;
    double d_norm_2=arma::prod(dnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t()));

    arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t());

    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*arma::det(H_eigen_inv)*arma::det(cor_clay_chol)*d_norm_2*BiCopPDF_cpp(1-p_norm_2(0),1-p_norm_2(1),13,par_clayton));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1)/llv_temp;

    llv+=log(llv_temp);


  }
  return llv;
}
// [[Rcpp::export]]
double loglike_Normal_ClaytonS(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_clayton =theta[2];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_clayton<= 0|| par_clayton>=100){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  
 //set_seed(3592);
 // arma::mat cor_clay=arma::cor(rbicop_cpp(100000,"clayton",180,par_clayton));
  arma::mat cor_clay=cor_Clayton(par_clayton);
  arma::mat cor_clay_chol=arma::chol(cor_clay).t();
  
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*arma::det(H_eigen_inv);;
    double d_norm_2=arma::prod(dnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*arma::det(H_eigen_inv)*arma::det(cor_clay_chol)*d_norm_2*claytonPDF(1-p_norm_2(0),1-p_norm_2(1),par_clayton));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1)/llv_temp;
    
    llv+=log(llv_temp);
    
    
  }
  
  return llv;
}
// [[Rcpp::export]]
double loglike_Normal_Clayton(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_clayton =theta[2];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_clayton<= 0|| par_clayton>=100){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  
  // set_seed(3592);
  // arma::mat cor_clay=arma::cor(rbicop_cpp(100000,"clayton",0,par_clayton));
  arma::mat cor_clay=cor_Clayton(par_clayton);
  arma::mat cor_clay_chol=arma::chol(cor_clay).t();
  
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*arma::det(H_eigen_inv);;
    double d_norm_2=arma::prod(dnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*arma::det(H_eigen_inv)*arma::det(cor_clay_chol)*d_norm_2*claytonPDF(p_norm_2(0),p_norm_2(1),par_clayton));
    //double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*arma::det(H_eigen_inv)*arma::det(cor_clay_chol)*d_norm_2*dbicop_cpp(p_norm_2,"clayton",0,par_clayton));
    
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1)/llv_temp;
    
    llv+=log(llv_temp);
    
    
  }
 
  return llv;
}

// 
// 
// // [[Rcpp::export]]
// double loglike_copula_Copula_3_rot(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r,arma::vec& type1, arma::vec& type2) {
//   // Log-Likelihood function
//   
//   // convert to matrices
//   int n = r.n_cols;
//   // Length of each series
//   int NoOBs = r.n_rows;
//   int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
//   arma::mat C = arma::zeros(n, n);
//   
//   int index = 0;
//   
//   for(int i = 0; i < n; i++){
//     for (int j = i; j < n; j++) {
//       C(j, i) = bekk(index);
//       index += 1;
//     }
//   }
//   
//   arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
//   arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
//   double p = theta[0];
//   double q = theta[1];
//   double par_c12 =theta[2];
//   double par_c23 =theta[3];
//   double par_c13 =theta[4];
//   
//   double par_c12_2 =theta[5];
//   double par_c23_2 =theta[6];
//   double par_c13_2 =theta[7];
//   
//   SEXP c12 = BiCop_cpp(type1[0],par_c12);
//   SEXP c23 = BiCop_cpp(type1[1],par_c23);
//   SEXP c13 = BiCop_cpp(type1[2],par_c13);
//   
//   SEXP c12_2 = BiCop_cpp(type2[0],par_c12_2);
//   SEXP c23_2 = BiCop_cpp(type2[1],par_c13_2);
//   SEXP c13_2 = BiCop_cpp(type2[2],par_c23_2);
//   
//   arma::mat rotation1=rot_mat(theta.subvec(8,10),n);
//   arma::mat rotation2=rot_mat(theta.subvec(11,13),n);
//   
//   //check for clayton or gumbel
//   arma::vec lower_bounds_1 = {0,0,0};
//   arma::vec lower_bounds_2 = {0,0,0};
//   for (int i=0; i<n;i++){
//     if(type1[i]==4){
//       lower_bounds_1[i]=1;
//     }
//     if(type2[i]==4){
//       lower_bounds_1[i]=1;
//     }
//   }
//   
//   //check for identification
//   if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1 || par_c12<=  lower_bounds_1[0]|| par_c12>=17 || par_c13<=  lower_bounds_1[1]|| par_c1>=17 || par_c23<=  lower_bounds_1[2]|| par_c12_>=17  || par_c12_2<=  lower_bounds_2[0]|| par_c12_2>=17 || par_c13_2<=  lower_bounds_2[1]|| par_c13_2>=17 || par_c23_2<=  lower_bounds_2[2]|| par_c23_2>=17){
//     return -1e25;
//   }
//   
//   
//   
//   
//   // compute inital H
//   arma::mat H = (r.t() * r) / r.n_rows;
//   
//   arma::mat CC  = C * C.t();
//   arma::mat At  = A.t();
//   arma::mat Gt  = G.t();
//   
//   
//   set_seed(3592);
//   arma::mat cor_c12=arma::cor(BiCopSim_cpp(100000,c12));
//   arma::mat cor_c13=arma::cor(BiCopSim_cpp(100000,c13));
//   arma::mat cor_c23=arma::cor(BiCopSim_cpp(100000,c23));
//   arma::mat cor_c12_2=arma::cor(BiCopSim_cpp(100000,c12_2));
//   arma::mat cor_c13_2=arma::cor(BiCopSim_cpp(100000,c13_2));
//   arma::mat cor_c23_2=arma::cor(BiCopSim_cpp(100000,c23_2));
//   
//   arma::mat cor_cop1=arma::eye(n,n);
//   cor_cop1(0,1)=cor_c12(0,1);
//   cor_cop(1,0)=cor_c12(0,1);
//   cor_cop1(0,2)=cor_c13(0,1);
//   cor_cop1(2,0)=cor_c13(0,1);
//   cor_cop1(1,2)=cor_c23(0,1);
//   cor_cop1(2,1)=cor_c23(0,1);
//   
//   arma::mat cor_cop2=arma::eye(n,n);
//   cor_cop2(0,1)=cor_c12_2(0,1);
//   cor_cop2(1,0)=cor_c12_2(0,1);
//   cor_cop2(0,2)=cor_c13_2(0,1);
//   cor_cop2(2,0)=cor_c13_2(0,1);
//   cor_cop2(1,2)=cor_c23_2(0,1);
//   cor_cop2(2,1)=cor_c23_2(0,1);
//   
//   arma::mat cor_cop_chol=arma::chol(cor_cop1).t();
//   
//   double cop_chol_det= arma::det(cor_cop_chol);
//   
//   arma::mat cor_cop_chol2=arma::chol(cor_cop2).t();
//   
//   
//   double cop_chol_det2= arma::det(cor_cop_chol2);
//   
//   
//   arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
//   
//   
//   
//   double p1t = 1;
//   
//   double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
//   // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
//   
//   arma::vec p_norm_1=pnorm_cpp(cor_cop_chol*rotation1*H_eigen_inv*r.row(0).t());
//   arma::vec p_norm_2=pnorm_cpp(cor_cop_chol2*rotation2*H_eigen_inv*r.row(0).t());
//   
//   double llv = log(arma::as_scalar((p1t*(1-p)+(1-p1t)*q)*BiCopPDF_cpp(BiCopHfunc2_cpp(p_norm_1[1],p_norm_1[2],C12),BiCopHfunc1_cpp(p_norm_1[2],p_norm_1[3],C23),C13)*BiCopPDF_cpp(p_norm_1[1],p_norm_1[2],C12)*BiCopPDF_cpp(p_norm_1[2],p_norm_1[3],C23)*d_norm_1*H_det*cop_chol_det));
//   
//   
//   
//   for (int i = 1; i < NoOBs; i++) {
//     H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
//     
//     arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
//     
//     
//     
//     double d_norm_1=arma::prod(dnorm_cpp(cor_cop_chol*rotation1*H_eigen_inv*r.row(i).t()));
//     double d_norm_2=arma::prod(dnorm_cpp(cor_clay_chol2*rotation2*H_eigen_inv*r.row(i).t()));
//     
//     arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*rotation1*H_eigen_inv*r.row(i).t());
//     arma::vec p_norm_2=pnorm_cpp(cor_clay_chol2*rotation2*H_eigen_inv*r.row(i).t());
//     
//     double llv_1 = arma::as_scalar((p1t*(p)+(1-p1t)*(1-q))*BiCopPDF_cpp(BiCopHfunc2_cpp(p_norm_1[1],p_norm_1[2],C12),BiCopHfunc1_cpp(p_norm_1[2],p_norm_1[3],C23),C13)*BiCopPDF_cpp(p_norm_1[1],p_norm_1[2],C12)*BiCopPDF_cpp(p_norm_1[2],p_norm_1[3],C23)*d_norm_1*H_det*cop_chol_det);
//     double llv_2 = arma::as_scalar((p1t*(1-p)+(1-p1t)*q)*BiCopPDF_cpp(BiCopHfunc2_cpp(p_norm_2[1],p_norm_2[2],C12_2),BiCopHfunc1_cpp(p_norm_2[2],p_norm_2[3],C23_2),C13_2)*BiCopPDF_cpp(p_norm_2[1],p_norm_2[2],C12_2)*BiCopPDF_cpp(p_norm_2[2],p_norm_2[3],C23_2)*d_norm_2*H_det*cop_chol_det2);
//     
//     double llv_temp = llv_1+llv_2
//       p1t=llv_1/llv_temp;
//     
//     llv+=log(llv_temp);
//     
//     
//   }
//   return llv;
// }
// 
// 
// 
// 
// // [[Rcpp::export]]
// double loglike_Normal_Copula_3_rot(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r,arma::vec& type) {
//   // Log-Likelihood function
//   
//   // convert to matrices
//   int n = r.n_cols;
//   // Length of each series
//   int NoOBs = r.n_rows;
//   int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
//   arma::mat C = arma::zeros(n, n);
//   
//   int index = 0;
//   
//   for(int i = 0; i < n; i++){
//     for (int j = i; j < n; j++) {
//       C(j, i) = bekk(index);
//       index += 1;
//     }
//   }
//   
//   arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
//   arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
//   double p = theta[0];
//   double q = theta[1];
//   double par_c12 =theta[2];
//   double par_c23 =theta[3];
//   double par_c13 =theta[4];
//   
//   arma::mat rotation=rot_mat(theta.subvec(5,7),n);
//   
//   
//   
//   SEXP c12 = BiCop_cpp(type[0],par_c12);
//   SEXP c23 = BiCop_cpp(type[1],par_c23);
//   SEXP c13 = BiCop_cpp(type[2],par_c13);
//   
//   //check for identification
//   if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_clayton<= 0|| par_clayton>=17){
//     return -1e25;
//   }
//   
//   
//   
//   // compute inital H
//   arma::mat H = (r.t() * r) / r.n_rows;
//   
//   arma::mat CC  = C * C.t();
//   arma::mat At  = A.t();
//   arma::mat Gt  = G.t();
//   
//   BiCopPDF(BiCopHfunc2(pnorm(series2[1]),pnorm(series2[2]),C_2_12),BiCopHfunc1(pnorm(series2[2]),pnorm(series2[3]),C_2_23),C_2_13)*BiCopPDF(pnorm(series2[1]),pnorm(series2[2]),C_2_12)*BiCopPDF(pnorm(series2[2]),pnorm(series2[3]),C_2_23)*dnorm(series2[1])*dnorm(series2[2])^2*dnorm(series2[3])*abs(det(Decomposed_cov_matrix2)*determinant_H_inverse)
//     //compute covariance matrix
//     
//     set_seed(3592);
//   arma::mat cor_c12=arma::cor(BiCopSim_cpp(100000,c12));
//   arma::mat cor_c13=arma::cor(BiCopSim_cpp(100000,c13));
//   arma::mat cor_c23=arma::cor(BiCopSim_cpp(100000,c23));
//   
//   arma::mat cor_cop=arma::eye(n,n);
//   cor_cop(0,1)=cor_c12(0,1);
//   cor_cop(1,0)=cor_c12(0,1);
//   cor_cop(0,2)=cor_c13(0,1);
//   cor_cop(2,0)=cor_c13(0,1);
//   cor_cop(1,2)=cor_c23(0,1);
//   cor_cop(2,1)=cor_c23(0,1);
//   
//   arma::mat cor_cop_chol=arma::chol(cor_cop).t();
//   
//   double cop_chol_det= arma::det(cor_cop_chol);
//   
//   arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
//   
//   
//   
//   double p1t = 1;
//   
//   double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
//   // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
//   
//   arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
//   arma::vec p_norm_2=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(0).t());
//   
//   double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
//   
//   
//   for (int i = 1; i < NoOBs; i++) {
//     H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
//     
//     arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
//     
//     double H_det=arma::det(H_eigen_inv);
//     
//     double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
//     double d_norm_2=arma::prod(dnorm_cpp(cor_clay_chol*rotation*H_eigen_inv*r.row(i).t()));
//     
//     arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*rotation*H_eigen_inv*r.row(i).t());
//     
//     double llv_1 = arma::as_scalar((p1t*p+(1-p1t)*(1-q))*arma::det(H_eigen_inv)*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*arma::det(H_eigen_inv)*arma::det(cor_clay_chol)*d_norm_2*BiCopPDF_cpp(p_norm_2(0),p_norm_2(1),3,par_clayton));
//     double llv_2 = arma::as_scalar((p1t*(1-p)+(1-p1t)*q)*BiCopPDF_cpp(BiCopHfunc2_cpp(p_norm_2[1],p_norm_2[2],C12),BiCopHfunc1_cpp(p_norm_2[2],p_norm_2[3],C23),C13)*BiCopPDF_cpp(p_norm_2[1],p_norm_2[2],C12)*BiCopPDF_cpp(p_norm_2[2],p_norm_2[3],C23)*d_norm_2*H_det*cop_chol_det);
//     
//     double llv_temp = llv_1+llv_2
//       p1t=llv_1/llv_temp;
//     
//     llv+=log(llv_temp);
//     
//     
//   }
//   return llv;
// }


// [[Rcpp::export]]
double loglike_Normal_ClaytonS1(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_clayton =theta[2];
  //check for identification
  if(valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1|| par_clayton<= 0|| par_clayton>=100){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  
  set_seed(3592);
  arma::mat cor_clay=arma::cor(BiCopSim_cpp(100000,13,par_clayton));
  arma::mat cor_clay_chol=arma::chol(cor_clay).t();
  
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*arma::det(H_eigen_inv);;
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*arma::det(H_eigen_inv);;
    double d_norm_2=arma::prod(dnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_clay_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*arma::det(H_eigen_inv)*arma::det(cor_clay_chol)*d_norm_2*BiCopPDF_cpp(p_norm_2(0),p_norm_2(1),13,par_clayton));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*d_norm_1)/llv_temp;
    
    llv+=log(llv_temp);
    
    
  }
  return llv;
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_clayton_gumbel(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=0;
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  double best_val = loglike_Clayton_Gumbel(bekk,theta,r);
  //set the seed
  
  
  // Generating random values theta
  while(l<(100/nc)){
    for(int i=0; i<(std::pow(n,2)-n+2);i++){
      theta_candidate[i]= arma::randn()*0.04+theta_mu[i];
  }
  double llv = loglike_Clayton_Gumbel(bekk,theta_candidate,r);
  if(llv > -1e25){
     l++;
   }
  if(llv>best_val){
     best_val=llv;
     theta_optim=theta_candidate;
     theta_mu=theta_optim;
      
      }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_gumbel(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=1;
  
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta;
  double llv = loglike_Normal_Gumbel(bekk,theta_candidate,r);
  double best_val = llv;
  //set the seed
  
  
  // Generating random values theta
  while(l<500){
    
      theta_candidate[0]= arma::randn()*0.05+theta_mu[0];
      theta_candidate[1]= arma::randn()*0.05+theta_mu[1];
      theta_candidate[2]= arma::randn()*8+theta_mu[2];
    
    llv = loglike_Normal_Gumbel(bekk,theta_candidate,r);
   
    
      l++;
    
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      
      
      theta_mu=theta_optim;
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}


// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_gumbel1(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=1;
  int j=0;
  
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta;
  
  double best_val = -1e24;
  //set the seed
  
  
  // Generating random values theta
  while(l<100 && j<5){
    
    theta_candidate[0]= arma::randn()*0.02+theta_mu[0];
    theta_candidate[1]= arma::randn()*0.02+theta_mu[1];
    theta_candidate[2]= arma::randn()*0.05+theta_mu[2];
    
    double llv = loglike_Normal_Gumbel1(bekk,theta_candidate,r);
    
    
    l++;
    
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      j++;
      theta_mu=theta_candidate;
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_gumbelS(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=1;
  int j=1;
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta;
  double llv = loglike_Normal_GumbelS(bekk,theta_candidate,r);
  double best_val = llv;
  //set the seed
  
  
  // Generating random values theta
  while(l<500){
    theta_candidate[0]= arma::randn()*0.05+theta_mu[0];
    theta_candidate[1]= arma::randn()*0.05+theta_mu[1];
    theta_candidate[2]= arma::randn()*8+theta_mu[2];
    llv = loglike_Normal_GumbelS(bekk,theta_candidate,r);
    
      l++;
    
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      
      
        theta_mu=theta_optim;
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_clayton(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=1;
  
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta;
  double llv = loglike_Normal_Clayton(bekk,theta_candidate,r);
  double best_val = llv;
  //set the seed
  
  
  // Generating random values theta
  while(l<500){
    theta_candidate[0]= arma::randn()*0.05+theta_mu[0];
    theta_candidate[1]= arma::randn()*0.05+theta_mu[1];
    theta_candidate[2]= arma::randn()*8+theta_mu[2];
    llv = loglike_Normal_Clayton(bekk,theta_candidate,r);
    
      l++;
    
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      
      
        theta_mu=theta_optim;
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_claytonS(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=1;
 
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta;
  double llv = loglike_Normal_ClaytonS(bekk,theta_candidate,r);
  double best_val = llv;
  //set the seed
  
  
  // Generating random values theta
  while(l<500){
    theta_candidate[0]= arma::randn()*0.05+theta_mu[0];
    theta_candidate[1]= arma::randn()*0.05+theta_mu[1];
    theta_candidate[2]= arma::randn()*8+theta_mu[2];
    llv = loglike_Normal_ClaytonS(bekk,theta_candidate,r);
    
      l++;
    
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      theta_mu=theta_optim;
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_frank(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=1;
  int j=1;
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta;
  double llv = loglike_Normal_Frank(bekk,theta_candidate,r);
  double best_val = llv;
  //set the seed
  
  
  // Generating random values theta
  while(l<100){
    theta_candidate[0]= arma::randn()*0.05+theta_mu[0];
    theta_candidate[1]= arma::randn()*0.05+theta_mu[1];
    theta_candidate[2]= arma::randn()*8+theta_mu[2];
    llv = loglike_Normal_Frank(bekk,theta_candidate,r);
    
      l++;
    
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      
        theta_mu=theta_optim;
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_clayton_gumbel90(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=0;
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  double best_val =loglike_Clayton_Gumbel90(bekk,theta_candidate,r);
  //set the seed
  
  
  // Generating random values theta
  while(l<(100/nc)){
    for(int i=0; i<std::pow(n,2)-n +2;i++){
      theta_candidate[i]= arma::randn()*0.04+theta_mu[i];
    }
    double llv = loglike_Clayton_Gumbel90(bekk,theta_candidate,r);
    if(llv > -1e25){
      l++;
    }
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      
      
        theta_mu=theta_optim;
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}



// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_gumbel_gumbelsurvival(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r) {
    int n =r.n_cols;
    int l=0;
    int j=1;
    
    arma::vec theta_mu;
    theta_mu=theta;
    arma::vec theta_optim=theta;
    arma::vec theta_candidate=theta; 
    arma::vec theta_sd={0.06,0.06,0.06,0.06,0.06,0.06,10,10};
    double best_val = loglike_Normal_Gumbel_GumbelSurvival(bekk,theta,r);
    
    
    // Generating random values theta
    while(l<200 && j<15){
      
      for(int i=0; i< theta.n_rows;i++){
        theta_candidate[i]= arma::randn()*theta_sd[i]+theta_mu[i];
      }
      double llv = loglike_Normal_Gumbel_GumbelSurvival(bekk,theta_candidate,r);
     
  
  
  
    l++;
    
  
      if(llv>best_val){
        Rcpp::Rcout << llv;
        best_val=llv;
        theta_optim=theta_candidate;
        j++;
        if(j==3 || j==6 || j>8){
        theta_mu=theta_optim;
        }
      }
    }
    return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                              Rcpp::Named("best_val") = best_val);
    
  }

// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_gumbelS_claytonS(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r) {
  int n =r.n_cols;
  int l=0;
  int j=1;
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  arma::vec theta_sd={0.06,0.06,0.06,0.06,0.06,0.06,10,10};
  double best_val = loglike_Normal_GumbelS_ClaytonS(bekk,theta,r);
  
  
  // Generating random values theta
  while(l<200 && j<15){
    
    for(int i=0; i< theta.n_rows;i++){
      theta_candidate[i]= arma::randn()*theta_sd[i]+theta_mu[i];
    }
    double llv = loglike_Normal_GumbelS_ClaytonS(bekk,theta_candidate,r);
    
    
    
    
    l++;
    
    
    if(llv>best_val){
      Rcpp::Rcout << llv;
      best_val=llv;
      theta_optim=theta_candidate;
      j++;
      if(j==3 || j==6 || j>8){
        theta_mu=theta_optim;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_gumbel_clayton(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r) {
  int n =r.n_cols;
  int l=0;
  int j=1;
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  arma::vec theta_sd={0.06,0.06,0.06,0.06,0.06,0.06,10,10};
  double best_val = loglike_Normal_Gumbel_Clayton(bekk,theta,r);
  
  
  // Generating random values theta
  while(l<200 && j<15){
    
    for(int i=0; i< theta.n_rows;i++){
      theta_candidate[i]= arma::randn()*theta_sd[i]+theta_mu[i];
    }
    double llv = loglike_Normal_Gumbel_Clayton(bekk,theta_candidate,r);
    
    
    
    
    l++;
    
    
    if(llv>best_val){
      Rcpp::Rcout << llv;
      best_val=llv;
      theta_optim=theta_candidate;
      j++;
      if(j==3 || j==6 || j>8){
        theta_mu=theta_optim;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}


// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_clayton_claytonsurvival(const arma::vec& bekk,const arma::vec& theta,const arma::mat& r) {
  int n =r.n_cols;
  int l=0;
  int j=1;
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  arma::vec theta_sd={0.06,0.06,0.06,0.06,0.06,0.06,10,10};
  double best_val =  loglike_Normal_Clayton_ClaytonSurvival(bekk,theta,r);
  
  // Generating random values theta
  while(l<200){
   
    for(int i=0; i< theta.n_rows;i++){
      theta_candidate[i]= arma::randn()*theta_sd[i]*(1/j)+theta_mu[i];
    }
    double llv = loglike_Normal_Clayton_ClaytonSurvival(bekk,theta_candidate,r);
    
    
  
      l++;
  
      
      if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      j++;
      if(j==3 || j==6 || j>8){
        theta_mu=theta_optim;
      }
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
arma::mat FilterProbs_normal_clayton_claytonS(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  double par_gumbel =theta[6];
  double par_gumbelSurvival =theta[7];
  
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"clayton",0,par_gumbel));
  arma::mat cor_gumbel=cor_Clayton(par_gumbel);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double gumbel_det = arma::det(cor_gumbel_chol);
  
  // set_seed(3592);
  // arma::mat cor_gumbelSurvival=arma::cor(rbicop_cpp(100000,"clayton",180,par_gumbelSurvival));
  arma::mat cor_gumbelSurvival=cor_Clayton(par_gumbelSurvival);
  arma::mat cor_gumbelSurvival_chol=arma::chol(cor_gumbelSurvival).t();
  
  double gumbelSurvival_det = arma::det(cor_gumbelSurvival_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  arma::mat filterProbs = arma::zeros(NoOBs,3);
  filterProbs(0,0)=1;
  double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(0).t());
  
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(i).t());
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*H_det );
    double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*claytonPDF(p_norm_1(0),p_norm_1(1),par_gumbel)*H_det* gumbel_det);
    double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*claytonPDF(1-p_norm_2(0),1-p_norm_2(1),par_gumbelSurvival)*H_det*gumbelSurvival_det);
    
    // double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*dbicop_cpp(p_norm_1,"clayton",0,par_gumbel)*H_det* gumbel_det);
    // double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*dbicop_cpp(p_norm_2,"clayton",180,par_gumbelSurvival)*H_det*gumbelSurvival_det);
    
    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2+llv_temp3);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = llv_temp3/llv_temp;
    
    filterProbs(i,0)=p1t;
    filterProbs(i,1)=p2t;
    filterProbs(i,2)=p3t;
    
    
  }
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat  FilterProbs_normal_gumbel_clayton(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  double par_gumbel =theta[6];
  double par_clayton =theta[7];
  
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",0,par_gumbel));
  
   arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double gumbel_det = arma::det(cor_gumbel_chol);
  
  // set_seed(3592);
  // arma::mat cor_clayton=arma::cor(rbicop_cpp(100000,"clayton",0,par_clayton));
  arma::mat cor_clayton=cor_Clayton(par_clayton);
  arma::mat cor_clayton_chol=arma::chol(cor_clayton).t();
  
  double clayton_det = arma::det(cor_clayton_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double H_det=arma::det(H_eigen_inv);
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  
  arma::mat filterProbs = arma::zeros(NoOBs,3);
  filterProbs(0,0)=1;
  
  double llv=log(arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*H_det);;
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t());
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*H_det);
    //double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*dbicop_cpp(p_norm_1,"gumbel",0,par_gumbel)*H_det* gumbel_det);
    //double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*dbicop_cpp(p_norm_2,"clayton",0,par_clayton)*H_det*clayton_det);
    double llv_temp2 = (p1t*p12+p2t*p22+p3t*p32)*d_norm_1*gumbelPDF(p_norm_1(0),p_norm_1(1),par_gumbel)*H_det* gumbel_det;
    double llv_temp3 = (p1t*p13+p2t*p23+p3t*p33)*d_norm_2*claytonPDF(p_norm_2(0),p_norm_2(1),par_clayton)*H_det*clayton_det;
    
    double llv_temp = llv_temp1+llv_temp2+llv_temp3;
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = 1-p1t-p2t;
    
    filterProbs(i,0)=p1t;
    filterProbs(i,1)=p2t;
    filterProbs(i,2)=p3t;
    
    
    
  }
  
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat FilterProbs_normal_gumbelS_claytonS(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  double par_gumbel =theta[6];
  double par_clayton =theta[7];
  
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",180,par_gumbel));
  //   
  arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double gumbel_det = arma::det(cor_gumbel_chol);
  
  // set_seed(3592);
  // arma::mat cor_clayton=arma::cor(rbicop_cpp(100000,"clayton",180,par_clayton));
  arma::mat cor_clayton=cor_Clayton(par_clayton);
  arma::mat cor_clayton_chol=arma::chol(cor_clayton).t();
  
  double clayton_det = arma::det(cor_clayton_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double H_det=arma::det(H_eigen_inv);
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  
  arma::mat filterProbs = arma::zeros(NoOBs,3);
  filterProbs(0,0)=1;
  
  double llv=log(arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()))*H_det);;
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t());
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*H_det);
    //double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*dbicop_cpp(p_norm_1,"gumbel",180,par_gumbel)*H_det* gumbel_det);
    //double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*dbicop_cpp(p_norm_2,"clayton",180,par_clayton)*H_det*clayton_det);
    double llv_temp2 = (p1t*p12+p2t*p22+p3t*p32)*d_norm_1*gumbelPDF(1-p_norm_1(0),1-p_norm_1(1),par_gumbel)*H_det* gumbel_det;
    double llv_temp3 = (p1t*p13+p2t*p23+p3t*p33)*d_norm_2*claytonPDF(1-p_norm_2(0),1-p_norm_2(1),par_clayton)*H_det*clayton_det;
    
    double llv_temp = llv_temp1+llv_temp2+llv_temp3;
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = 1-p1t-p2t;
    
    filterProbs(i,0)=p1t;
    filterProbs(i,1)=p2t;
    filterProbs(i,2)=p3t;
    
    
    
    
  }
  
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat FilterProbs_normal_gumbel_gumbelS(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  double par_gumbel =theta[6];
  double par_gumbelSurvival =theta[7];
  
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",0,par_gumbel));
  arma::mat cor_gumbel=cor_Gumbel(par_gumbel);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double gumbel_det = arma::det(cor_gumbel_chol);
  
  // set_seed(3592);
  // arma::mat cor_gumbelSurvival=arma::cor(rbicop_cpp(100000,"gumbel",180,par_gumbelSurvival));
  arma::mat cor_gumbelSurvival=cor_Gumbel(par_gumbelSurvival);
  arma::mat cor_gumbelSurvival_chol=arma::chol(cor_gumbelSurvival).t();
  
  double gumbelSurvival_det = arma::det(cor_gumbelSurvival_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  arma::mat filterProbs = arma::zeros(NoOBs,3);
  filterProbs(0,0)=1;
  double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(0).t());
  
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(i).t());
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()))*H_det );
    // double llv_temp2 = ((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*dbicop_cpp(p_norm_1,"gumbel",0,par_gumbel)*H_det* gumbel_det);
    // double llv_temp3 = ((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",180,par_gumbelSurvival)*H_det*gumbelSurvival_det);
    
    double llv_temp2 = ((p1t*p12+p2t*p22+p3t*p32)*d_norm_1*gumbelPDF(p_norm_1(0),p_norm_1(1),par_gumbel)*H_det* gumbel_det);
    double llv_temp3 = ((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*gumbelPDF(1-p_norm_2(0),1-p_norm_2(1),par_gumbelSurvival)*H_det*gumbelSurvival_det);
    
    double llv_temp = llv_temp1+llv_temp2+llv_temp3;
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = llv_temp3/llv_temp;
    
    filterProbs(i,0)=p1t;
    filterProbs(i,1)=p2t;
    filterProbs(i,2)=p3t;
    
    
  }
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat FilterProbs_Clayton_Gumbel90(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  arma::mat filterProbs = arma::zeros(NoOBs,1);
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_clayton =theta[2];
  double par_gumbel =theta[3];
  //check for identification
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  set_seed(3592);
  arma::mat cor_clayton=arma::cor(BiCopSim_cpp(100000,3,par_clayton));
  arma::mat cor_clayton_chol=arma::chol(cor_clayton).t();
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,24,par_gumbel));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double clay_det=arma::det(cor_clayton_chol);
  double gumb_det=arma::det(cor_gumbel_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  filterProbs(0,0)=1;
  
  double d_norm_1=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*clay_det*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),3,par_clayton)));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double H_det =arma::det(H_eigen_inv);
    
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_clayton_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*clay_det*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),3,par_clayton)+(p1t*(1-p)+(1-p1t)*(q))*H_det*gumb_det*d_norm_2*BiCopPDF_cpp(p_norm_2(0),p_norm_2(1),24,par_gumbel));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*clay_det*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),3,par_clayton))/llv_temp;
    
    filterProbs(i,0)=p1t;
    
    
  }
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat FilterProbs_Normal_Gumbel(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  arma::mat filterProbs = arma::zeros(NoOBs,1);
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }

  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_copula =theta[2];

  //check for identification


  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;

  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();


  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"gumbel",0,par_copula));
  arma::mat cor_gumbel=cor_Gumbel(par_copula);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();


  double gumb_det=arma::det(cor_gumbel_chol);

  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));



  double p1t = 1;

  filterProbs(0,0)=1;

  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));

  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());

  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));


  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;

    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double H_det =arma::det(H_eigen_inv);


    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));

     arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());

    double llv_temp =((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*H_det*gumb_det*d_norm_2*gumbelPDF(p_norm_2(0),p_norm_2(1),par_copula));
    p1t=((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1)/llv_temp;

    filterProbs(i,0)=p1t;


  }
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat FilterProbs_Normal_Clayton(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  arma::mat filterProbs = arma::zeros(NoOBs,1);
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_copula =theta[2];
  
  //check for identification
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_gumbel=arma::cor(rbicop_cpp(100000,"clayton",0,par_copula));
  arma::mat cor_gumbel=cor_Clayton(par_copula);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  
  double gumb_det=arma::det(cor_gumbel_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  filterProbs(0,0)=1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double H_det =arma::det(H_eigen_inv);
    
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*H_det*gumb_det*d_norm_2*claytonPDF(p_norm_2(0),p_norm_2(1),par_copula));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1)/llv_temp;
    
    filterProbs(i,0)=p1t;
    
    
  }
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat FilterProbs_Normal_GumbelS(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  arma::mat filterProbs = arma::zeros(NoOBs,1);
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_copula =theta[2];
  
  //check for identification
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  //set_seed(3592);
  //arma::mat cor_gumbel = arma::cor(rbicop_cpp(100000,"gumbel",180,par_copula));
  arma::mat cor_gumbel=cor_Gumbel(par_copula);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  
  double gumb_det=arma::det(cor_gumbel_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  filterProbs(0,0)=1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double H_det =arma::det(H_eigen_inv);
    
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp =((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*H_det*gumb_det*d_norm_2*gumbelPDF(1-p_norm_2(0),1-p_norm_2(1),par_copula));
    //double llv_temp =((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*H_det*gumb_det*d_norm_2*dbicop_cpp(p_norm_2,"gumbel",180,par_copula));
    
    p1t=((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1)/llv_temp;
    
    filterProbs(i,0)=p1t;
    
    
  }
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat FilterProbs_Normal_ClaytonS(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  arma::mat filterProbs = arma::zeros(NoOBs,1);
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_copula =theta[2];
  
  //check for identification
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  // set_seed(3592);
  // arma::mat cor_gumbel = arma::cor(rbicop_cpp(100000,"clayton",180,par_copula));
  arma::mat cor_gumbel=cor_Clayton(par_copula);
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  
  double gumb_det=arma::det(cor_gumbel_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  filterProbs(0,0)=1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double H_det =arma::det(H_eigen_inv);
    
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*H_det*gumb_det*d_norm_2*claytonPDF(1-p_norm_2(0),1-p_norm_2(1),par_copula));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1)/llv_temp;
    
    filterProbs(i,0)=p1t;
    
    
  }
  return filterProbs;
}
// [[Rcpp::export]]
arma::mat FilterProbs_Normal_Frank(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  arma::mat filterProbs = arma::zeros(NoOBs,1);
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_copula =theta[2];
  
  //check for identification
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  set_seed(3592);
  arma::mat cor_gumbel = arma::cor(rbicop_cpp(100000,"frank",0,par_copula));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  
  double gumb_det=arma::det(cor_gumbel_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  filterProbs(0,0)=1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    double H_det =arma::det(H_eigen_inv);
    
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp =arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1+(p1t*(1-p)+(1-p1t)*(q))*H_det*gumb_det*d_norm_2*dbicop_cpp(p_norm_2,"frank",0,par_copula));
    p1t=arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1)/llv_temp;
    
    filterProbs(i,0)=p1t;
    
    
  }
  return filterProbs;
}















// [[Rcpp::export]]
double loglike_LL_Gumbel(const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(theta.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(theta.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double par_gumbel = theta[11];
  //check for identification
  if(valid_bekk(C,A,G)==false || par_gumbel<= 1|| par_gumbel>=17){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  set_seed(3592);
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,4,par_gumbel));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double det_gumbel = arma::det(cor_gumbel_chol);
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),4,par_gumbel)));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    llv+=log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),4,par_gumbel)));
    
    
    
    
  }
  return llv;
}


// [[Rcpp::export]]
double loglike_LL_GumbelS(const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(theta.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(theta.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double par_gumbel = theta[11];
  //check for identification
  if(valid_bekk(C,A,G)==false || par_gumbel<= 1|| par_gumbel>=17){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  set_seed(3592);
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,14,par_gumbel));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double det_gumbel = arma::det(cor_gumbel_chol);
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),14,par_gumbel)));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    llv+=log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),14,par_gumbel)));
    
    
    
    
  }
  return llv;
}
// [[Rcpp::export]]
double loglike_LL_Clayton(const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(theta.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(theta.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double par_clayton = theta[11];
  //check for identification
  if(valid_bekk(C,A,G)==false || par_clayton<=0 || par_clayton>=17){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  set_seed(3592);
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,3,par_clayton));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double det_gumbel = arma::det(cor_gumbel_chol);
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),3,par_clayton)));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    llv+=log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),3,par_clayton)));
    
    
    
    
  }
  return llv;
}

// [[Rcpp::export]]
double loglike_LL_ClaytonS(const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(theta.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(theta.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double par_clayton = theta[11];
  //check for identification
  if(valid_bekk(C,A,G)==false || par_clayton<=0 || par_clayton>=17){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  set_seed(3592);
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,13,par_clayton));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double det_gumbel = arma::det(cor_gumbel_chol);
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),13,par_clayton)));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    llv+=log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),13,par_clayton)));
    
    
    
    
  }
  return llv;
}

// [[Rcpp::export]]
double loglike_LL_Frank(const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(theta.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat G = arma::reshape(theta.subvec((index +  std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double par_frank = theta[11];
  //check for identification
  if(valid_bekk(C,A,G)==false || par_frank<=-35 || par_frank ==0|| par_frank>=35){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  set_seed(3592);
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(100000,5,par_frank));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double det_gumbel = arma::det(cor_gumbel_chol);
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),5,par_frank)));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    
    
    double d_norm=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    
    llv+=log(arma::as_scalar(arma::det(H_eigen_inv)*det_gumbel*d_norm*BiCopPDF_cpp(p_norm(0),p_norm(1),5,par_frank)));
    
    
    
    
  }
  return llv;
}
// std::function<double(double)> f(arma::mat BEKK, arma::mat r) {
//   auto ax2 = (double x) {  return (-loglike_Normal_Gumbel_GumbelSurvival(BEKK,x,r)); };
//   return ax2;
// }
// 
// // [[Rcpp::export]]
// Rcpp::List optim_3state(arma::vec init_val,arma::mat BEKK, arma::mat r)
// {
//   double fopt;
//   int res = optim_lbfgs(f(BEKK,r), init_val, fopt);
//   return Rcpp::List::create(
//     Rcpp::Named("xopt") = init_val,
//     Rcpp::Named("fopt") = fopt,
//     Rcpp::Named("status") = res
//   );
// }
// 
