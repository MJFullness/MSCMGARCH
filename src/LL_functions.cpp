
#include <RcppArmadillo.h>
#include "LL_2dim.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]




// [[Rcpp::export]]
int indicatorFunction(arma::mat r,arma::mat signs){
  r=r.t();
  
  int indicator=1;
  int n=r.n_rows;
  for (int i=0; i<n; i++){
    if(arma::as_scalar(signs.row(i)) * arma::as_scalar(r.row(i)) < 0){
      indicator = 0;
    }
  }
  return indicator;
}

// [[Rcpp::export]]
arma::mat comp_bekk_forecast(const arma::vec& bekk, const arma::mat& r) {
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
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  for (int i = 1; i <= NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
  }
  return H;
}

// [[Rcpp::export]]
arma::mat comp_asymm_bekk_forecast(const arma::vec& bekk, const arma::mat& r, arma::mat signs) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 3 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }

  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat B = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), index +  2*std::pow(n, 2)-1).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  2*std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;

  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Bt  = B.t();
  arma::mat Gt  = G.t();


  for (int i = 1; i <= NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + indicatorFunction(r.row(i - 1), signs) * Bt * r.row(i - 1).t() * r.row(i - 1) * B + Gt * H * G;

  }
  return H;
}


// [[Rcpp::export]]
bool valid_asymm_bekk(arma::mat& C,arma::mat& A, arma::mat& B ,arma::mat& G, arma::mat r, arma::mat signs){
  int n =C.n_cols;
  int N =r.n_rows;
  double exp_indicator_value = 0;
  for (int i=0; i<N;i++){
    exp_indicator_value+=indicatorFunction(r.row(i),signs);
    
  }
  exp_indicator_value=exp_indicator_value/N;
  
  arma::mat prod = kron(A,A)+exp_indicator_value*kron(B,B)+kron(G,G);
  
  arma::vec eigvals;
  eigvals= abs(arma::eig_gen(prod));
  double max=0;
  for (int i=0; i< eigvals.n_elem; i++){
    if(eigvals[i]>max){
      max=eigvals[i];
    }
  }
  if(max >= 1){
    return false;
  }
  
  for (int i=0; i<n;i++){
    if(C(i,i)<=0){
      return false;
    }
  }
  if(A(0,0)<=0 || B(0,0)<=0 || G(0,0)<=0) {
    return false;
  }
  
  else{
    return true;
  }
}


// [[Rcpp::export]]
arma::vec copula_cor(arma::vec& theta, arma::vec& type){
  int n = type.n_rows;
  arma::vec cor_vec = type;
  for (int i=0; i<n; i++){
    if(type[i]==3 ||type[i]==13 ){
      cor_vec[i]=cor_Clayton(theta[i])(0,1);
    }
    else{
      cor_vec[i]=cor_Gumbel(theta[i])(0,1);
    }
  }
  return cor_vec;
}

// [[Rcpp::export]]
arma::mat cor_mat(arma::vec& theta, arma::vec& type){
  int n = type.n_rows;
  arma::mat cor_mat = arma::eye(n,n);
  arma::vec cor_vec = type;
  for (int i=0; i<n; i++){
    if(type[i]==3 || type[i]==13){
      cor_vec[i]=cor_Clayton(theta[i])(0,1);
    }
    else{
      cor_vec[i]=cor_Gumbel(theta[i])(0,1);
    }
  }
  // works only for dim=3 yet
  cor_mat(0,1)=cor_vec[0]*sqrt(1-pow(cor_vec[1],2))*sqrt(1-pow(cor_vec[2],2))+cor_vec[1]*cor_vec[2];
  cor_mat(1,0)=cor_mat(0,1);
  cor_mat(0,2)=cor_vec[1];
  cor_mat(2,0)=cor_vec[1];
  cor_mat(1,2)=cor_vec[2];
  cor_mat(2,1)=cor_vec[2];
  return cor_mat;
}

// [[Rcpp::export]]
double loglike_Normal_Copula_3(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r,arma::vec& type) {
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
 
  
  //check for identification
  if( valid_bekk(C,A,G)==false || p<=0 || q<=0 || p>=1 || q>=1 || valid_copula(theta.subvec(2,4),type)==false){
    return -1e25;
  }
  
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
  arma::vec copula_par = theta.subvec(2,4);
  arma::mat cor_cop=cor_mat(copula_par,type);
  
  if(arma::det(cor_cop)<=0){
    return -1e25;
  }
  
  arma::mat cor_cop_chol=arma::chol(cor_cop).t();
  
  double cop_chol_det= arma::det(cor_cop_chol);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  //arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  //arma::vec p_norm_2=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t());
    
    double llv_1 = arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1);
    double llv_2 = arma::as_scalar((p1t*(1-p)+(1-p1t)*q)*VineCopula(p_norm_2,copula_par,type)*d_norm_2*H_det*cop_chol_det);
    
    double llv_temp = llv_1+llv_2;
      p1t=llv_1/llv_temp;
    
    llv+=log(llv_temp);
    
    
  }
  return llv;
}

// [[Rcpp::export]]
double loglike_Normal_Copula_3_asymm(const arma::vec& bekk, arma::vec signs, const arma::vec& theta, const arma::mat& r,arma::vec& type) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 3 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat B = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), index +  2*std::pow(n, 2)-1).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  2*std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  
  
  //check for identification
  if(valid_asymm_bekk(C,A,B,G,r,signs)==false || p<=0 || q<=0 || p>=1 || q>=1 || valid_copula(theta.subvec(2,4),type)==false){
    return -1e25;
  }
  
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Bt  = B.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
  arma::vec copula_par = theta.subvec(2,4);
 
  arma::mat cor_cop=cor_mat(copula_par,type);

  if(arma::det(cor_cop)<=0){
    return -1e26;
  }
  
  arma::mat cor_cop_chol=arma::chol(cor_cop).t();
  
  double cop_chol_det= arma::det(cor_cop_chol);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  
 
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + indicatorFunction(r.row(i - 1), signs) * Bt * r.row(i - 1).t() * r.row(i - 1) * B + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t());
    
    double llv_1 = arma::as_scalar((p1t*p+(1-p1t)*(1-q))*H_det*d_norm_1);
    double llv_2 = arma::as_scalar((p1t*(1-p)+(1-p1t)*q)*VineCopula(p_norm_2,copula_par,type)*d_norm_2*H_det*cop_chol_det);
    
    double llv_temp = llv_1+llv_2;
    p1t=llv_1/llv_temp;
    
    llv+=log(llv_temp);
    
    
  }
  return llv;
}


// [[Rcpp::export]]
double loglike_Normal_Copula_Copula_3(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r,arma::vec& type) {
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
  arma::vec copula_par = theta.subvec(6,11);
  //check for identification
  
  //check for identification
  if(valid_bekk(C,A,G)==false  || p13 >= 1||   p13 <= 0|| p23 >= 1|| p21+p22>=1 || p32+p33>=1 || p11+p12>=1 || p23 <= 0||  p31 >= 1||  p31 <= 0 || p11<=0 || p12<=0 || p11>=1 || p12>=1 || p21<=0 || p22<=0 || p21>=1 || p22>=1  || p33<=0 || p32<=0 || p33>=1 || p32>=1|| valid_copula(copula_par,type)==false){
    return -1e25;
  }
  arma::vec copula_par_1 = theta.subvec(6,8);
  arma::vec copula_par_2 = theta.subvec(9,11);
  arma::vec type1 = type.subvec(0,2);
  arma::vec type2 = type.subvec(3,5);
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
 
  
  arma::mat cor_cop_1=cor_mat(copula_par_1,type1);
  
  arma::mat cor_cop_2=cor_mat(copula_par_2,type2);
  
  
  if(arma::det(cor_cop_1)<=0 || arma::det(cor_cop_2)<=0){
    return -1e25;
  }
  
  arma::mat cor_cop_chol_1=arma::chol(cor_cop_1).t();
  
  double cop_chol_det_1= arma::det(cor_cop_chol_1);
  
  arma::mat cor_cop_chol_2=arma::chol(cor_cop_2).t();
  
  double cop_chol_det_2= arma::det(cor_cop_chol_2);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  //arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_cop_chol_1*H_eigen_inv*r.row(i).t()));
    double d_norm_3=arma::prod(dnorm_cpp(cor_cop_chol_2*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_cop_chol_1*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_3=pnorm_cpp(cor_cop_chol_2*H_eigen_inv*r.row(i).t());
    
    double llv_1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*H_det*d_norm_1);
    double llv_2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*VineCopula(p_norm_2,copula_par_1,type1)*d_norm_2*H_det*cop_chol_det_1);
    double llv_3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*VineCopula(p_norm_3,copula_par_2,type2)*d_norm_3*H_det*cop_chol_det_2);
    
    double llv_temp = llv_1+llv_2+llv_3;
    p1t=llv_1/llv_temp;
    p2t=llv_2/llv_temp;
    p3t=llv_3/llv_temp;
    
    llv+=log(llv_temp);
    
    
  }
  return llv;
}

// [[Rcpp::export]]
double loglike_Normal_Copula_Copula_3_asymm(const arma::vec& bekk, arma::vec signs, const arma::vec& theta, const arma::mat& r,arma::vec& type) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 3 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat B = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), index +  2*std::pow(n, 2)-1).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  2*std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  arma::vec copula_par = theta.subvec(6,11);
  //check for identification
  
  //check for identification
  if(valid_asymm_bekk(C,A,B,G,r,signs)==false  || p13 >= 1||   p13 <= 0|| p23 >= 1|| p21+p22>=1 || p32+p33>=1 || p11+p12>=1 || p23 <= 0||  p31 >= 1||  p31 <= 0 || p11<=0 || p12<=0 || p11>=1 || p12>=1 || p21<=0 || p22<=0 || p21>=1 || p22>=1  || p33<=0 || p32<=0 || p33>=1 || p32>=1|| valid_copula(copula_par,type)==false){
    return -1e25;
  }
  arma::vec copula_par_1 = theta.subvec(6,8);
  arma::vec copula_par_2 = theta.subvec(9,11);
  arma::vec type1 = type.subvec(0,2);
  arma::vec type2 = type.subvec(3,5);
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Bt  = B.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
  
  
  arma::mat cor_cop_1=cor_mat(copula_par_1,type1);
  
  arma::mat cor_cop_2=cor_mat(copula_par_2,type2);
  
  
  if(arma::det(cor_cop_1)<=0 || arma::det(cor_cop_2)<=0){
    return -1e25;
  }
  
  arma::mat cor_cop_chol_1=arma::chol(cor_cop_1).t();
  
  double cop_chol_det_1= arma::det(cor_cop_chol_1);
  
  arma::mat cor_cop_chol_2=arma::chol(cor_cop_2).t();
  
  double cop_chol_det_2= arma::det(cor_cop_chol_2);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  //arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + indicatorFunction(r.row(i - 1), signs) * Bt * r.row(i - 1).t() * r.row(i - 1) * B + Gt * H * G;
    
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_cop_chol_1*H_eigen_inv*r.row(i).t()));
    double d_norm_3=arma::prod(dnorm_cpp(cor_cop_chol_2*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_cop_chol_1*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_3=pnorm_cpp(cor_cop_chol_2*H_eigen_inv*r.row(i).t());
    
    double llv_1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*H_det*d_norm_1);
    double llv_2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*VineCopula(p_norm_2,copula_par_1,type1)*d_norm_2*H_det*cop_chol_det_1);
    double llv_3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*VineCopula(p_norm_3,copula_par_2,type2)*d_norm_3*H_det*cop_chol_det_2);
    
    double llv_temp = llv_1+llv_2+llv_3;
    p1t=llv_1/llv_temp;
    p2t=llv_2/llv_temp;
    p3t=llv_3/llv_temp;
    
    llv+=log(llv_temp);
    
    
  }
  return llv;
}



// [[Rcpp::export]]
double LL_loglike_Copula_3(arma::vec& theta, const arma::mat& r,arma::vec& type) {
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
  
  arma::mat A = arma::reshape(theta.subvec(index, (index + std::pow(n, 2)) - 1 ), n, n);
  arma::mat G = arma::reshape(theta.subvec((index +  std::pow(n, 2)), numb_of_vars-1), n, n);
 
  
  //check for identification
  if(valid_bekk(C,A,G)==false || valid_copula(theta.subvec(numb_of_vars,numb_of_vars+2),type)==false){
    return -1e25;
  }
  
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
  arma::vec copula_par = theta.subvec(numb_of_vars,numb_of_vars+2);
  arma::mat cor_cop = cor_mat(copula_par,type);
  
  
  if(arma::det(cor_cop)<=0){
    return -1e25;
  }
  
  arma::mat cor_cop_chol=arma::chol(cor_cop).t();
  
  double cop_chol_det= arma::det(cor_cop_chol);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  double H_det=arma::det(H_eigen_inv);
  
  
  double d_norm_1=arma::prod(dnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  
  arma::vec p_norm_1=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(VineCopula(p_norm_1,copula_par,type)*d_norm_1*H_det*cop_chol_det));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
      double d_norm_1=arma::prod(dnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp = log(arma::as_scalar(VineCopula(p_norm_1,copula_par,type)*d_norm_1*H_det*cop_chol_det));
    
    
    llv+=llv_temp;
    
    
    
  }
  return llv;
}


// [[Rcpp::export]]
double LL_loglike_Copula_3_asymm(const arma::vec& theta, arma::vec signs, const arma::mat& r,arma::vec& type) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 3 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(theta.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat B = arma::reshape(theta.subvec((index +  std::pow(n, 2)), index +  2*std::pow(n, 2)-1).t(), n, n);
  arma::mat G = arma::reshape(theta.subvec((index +  2*std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  
  
  
  //check for identification
  if(valid_asymm_bekk(C,A,B,G,r,signs)==false || valid_copula(theta.subvec(numb_of_vars,numb_of_vars+2),type)==false){
    return -1e25;
  }
  
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Bt  = B.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
  arma::vec copula_par = theta.subvec(numb_of_vars,numb_of_vars+2);
  arma::mat cor_cop = cor_mat(copula_par,type);
  
  
  if(arma::det(cor_cop)<=0){
    return -1e25;
  }
  
  arma::mat cor_cop_chol=arma::chol(cor_cop).t();
  
  double cop_chol_det= arma::det(cor_cop_chol);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  double H_det=arma::det(H_eigen_inv);
  
  
  double d_norm_1=arma::prod(pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  
  arma::vec p_norm_1=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(VineCopula(p_norm_1,copula_par,type)*d_norm_1*H_det*cop_chol_det));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + indicatorFunction(r.row(i - 1), signs) * Bt * r.row(i - 1).t() * r.row(i - 1) * B + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t());
    
    double llv_temp = log(arma::as_scalar(VineCopula(p_norm_1,copula_par,type)*d_norm_1*H_det*cop_chol_det));
    
    
    llv+=llv_temp;
    
    
  }
  return llv;
}


// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_copula_3(const arma::vec& bekk,const arma::vec& theta,arma::vec& type,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=0;
  
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  double best_val = loglike_Normal_Copula_3(bekk,theta,r,type);
  arma::vec alternative_theta = {0.9,0.05,15,15,15};
  double best_val_alt = loglike_Normal_Copula_3(bekk,alternative_theta,r,type);
  
  
  
  //set the seed
  if(best_val_alt>best_val){
    best_val=best_val_alt;
    theta_mu=alternative_theta;
    theta_candidate=alternative_theta;
  }
  //set the seed
  
  
  // Generating random values theta
  while(l<200){
    for(int i=0; i<2;i++){
      theta_candidate[i]= arma::randn()*0.001+theta_mu[i];
    }
    for(int i=2; i<5;i++){
      theta_candidate[i]= arma::randn()+theta_mu[i];
    }
    double llv = loglike_Normal_Copula_3(bekk,theta_candidate,r,type);
    
    if(llv > -1e25){
      l++;
    }
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      if(l>20){
        theta_mu=theta_optim;
        
      }
        
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_copula_copula_3(const arma::vec& bekk,const arma::vec& theta,arma::vec& type,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=0;
  
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  
  
  //set the seed
  double best_val = loglike_Normal_Copula_Copula_3(bekk,theta,r,type);
  
  
  // Generating random values theta
  while(l<300){
    
    theta_candidate[0]=arma::randu(); 
    theta_candidate[1]=arma::randu()*(1-theta_candidate[0]); 
    theta_candidate[2]=arma::randu();
    theta_candidate[3]=arma::randu()*(1-theta_candidate[2]); 
    theta_candidate[4]=arma::randu();
    theta_candidate[5]=arma::randu()*(1-theta_candidate[4]); 
    
      
    
    for(int i=6; i<12;i++){
      theta_candidate[i]= arma::randn()+theta_mu[i];
    }
    double llv = loglike_Normal_Copula_Copula_3(bekk,theta_candidate,r,type);
   
    if(llv > -1e25){
      l++;
    }
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      
      
      if(l>20){
        theta_mu=theta_optim;
      }
      
      
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}

// [[Rcpp::export]]
Rcpp::List random_grid_search_normal_copula_3_asymm(const arma::vec& bekk,arma::vec signs,const arma::vec& theta,arma::vec& type,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=0;
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  arma::vec alternative_theta = {0.9,0.1,15,15,15};
  double best_val_alt = loglike_Normal_Copula_3_asymm(bekk,signs,alternative_theta,r,type);
  
  
  double best_val = loglike_Normal_Copula_3_asymm(bekk,signs,theta,r,type);
  //set the seed
  if(best_val_alt>best_val){
    best_val=best_val_alt;
    theta_mu=alternative_theta;
    theta_candidate=alternative_theta;
  }
  
  // Generating random values theta
  while(l<200){
    for(int i=0; i<2;i++){
      theta_candidate[i]= arma::randn()*0.005+theta_mu[i];
    }
    for(int i=2; i<5;i++){
      theta_candidate[i]= arma::randn()+theta_mu[i];
    }
    double llv = loglike_Normal_Copula_3_asymm(bekk,signs,theta_candidate,r,type);
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
Rcpp::List random_grid_search_normal_copula_copula_3_asymm(const arma::vec& bekk,arma::vec signs,const arma::vec& theta,arma::vec& type,const arma::mat& r, const int nc) {
  int n =r.n_cols;
  int l=0;
  
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  arma::vec theta_alternative = theta;
 
  double best_val = loglike_Normal_Copula_Copula_3_asymm(bekk,signs,theta,r,type);
  //set the seed
  
  
  
  // Generating random values theta
  while(l<300){
    theta_candidate[0]=arma::randu(); 
    theta_candidate[1]=arma::randu()*(1-theta_candidate[0]); 
    theta_candidate[2]=arma::randu();
    theta_candidate[3]=arma::randu()*(1-theta_candidate[2]); 
    theta_candidate[4]=arma::randu();
    theta_candidate[5]=arma::randu()*(1-theta_candidate[4]); 
    for(int i=6; i<12;i++){
      theta_candidate[i]= arma::randn()+theta_mu[i];
    }
    double llv = loglike_Normal_Copula_Copula_3_asymm(bekk,signs,theta_candidate,r,type);
    if(llv > -1e25){
      l++;
    }
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      if(l>20){
        theta_mu=theta_optim;
        
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}


// [[Rcpp::export]]
arma::mat FilterProbs_Normal_Copula_3(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r,arma::vec& type) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  arma::mat filterProbs = arma::zeros(NoOBs,1);
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
  
  
  
  
  
  filterProbs(0,0)=1;
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
  arma::vec copula_par = theta.subvec(2,4);
 
  arma::mat cor_cop=cor_mat(copula_par,type);
  
  arma::mat cor_cop_chol=arma::chol(cor_cop).t();
  
  double cop_chol_det= arma::det(cor_cop_chol);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t());
    
    double llv_1 = arma::as_scalar((p1t*p+(1-p1t)*(1-q))*arma::det(H_eigen_inv)*d_norm_1);
    double llv_2 = arma::as_scalar((p1t*(1-p)+(1-p1t)*q)*VineCopula(p_norm_2,copula_par,type)*d_norm_2*H_det*cop_chol_det);
    
    double llv_temp = llv_1+llv_2;
    p1t=llv_1/llv_temp;
    
    llv+=log(llv_temp);
    filterProbs(i,0)=p1t;
    
    
  }
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat FilterProbs_Normal_Copula_3_asymm(const arma::vec& bekk, arma::vec signs, const arma::vec& theta, const arma::mat& r,arma::vec& type) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 3 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  arma::mat filterProbs = arma::zeros(NoOBs,1);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat B = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), index +  2*std::pow(n, 2)-1).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  2*std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p = theta[0];
  double q = theta[1];
  double par_c12 =theta[2];
  double par_c13 =theta[3];
  double par_c23 =theta[4];
  
  
  filterProbs(0,0)=1;
  
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Bt  = B.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
  arma::vec copula_par = theta.subvec(2,4);

  arma::mat cor_cop=cor_mat(copula_par,type);
  
  
  arma::mat cor_cop_chol=arma::chol(cor_cop).t();
  
  double cop_chol_det= arma::det(cor_cop_chol);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + indicatorFunction(r.row(i - 1), signs) * Bt * r.row(i - 1).t() * r.row(i - 1) * B + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_cop_chol*H_eigen_inv*r.row(i).t());
    
    double llv_1 = arma::as_scalar((p1t*p+(1-p1t)*(1-q))*arma::det(H_eigen_inv)*d_norm_1);
    double llv_2 = arma::as_scalar((p1t*(1-p)+(1-p1t)*q)*VineCopula(p_norm_2,copula_par,type)*d_norm_2*H_det*cop_chol_det);
    
    double llv_temp = llv_1+llv_2;
    p1t=llv_1/llv_temp;
    
    llv+=log(llv_temp);
    
    filterProbs(i,0)=p1t;
  }
  return filterProbs;
}


// [[Rcpp::export]]
arma::mat FilterProbs_Normal_Copula_Copula_3(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r,arma::vec& type) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  arma::mat filterProbs = arma::zeros(NoOBs,3);
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
  arma::vec copula_par = theta.subvec(6,11);
  //check for identification
  
  
  arma::vec copula_par_1 = theta.subvec(6,8);
  arma::vec copula_par_2 = theta.subvec(9,11);
  arma::vec type1 = type.subvec(0,2);
  arma::vec type2 = type.subvec(3,5);
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
  
  
  arma::mat cor_cop_1=cor_mat(copula_par_1,type1);

  
  
  arma::mat cor_cop_2=cor_mat(copula_par_2,type2);
  
  
  
  arma::mat cor_cop_chol_1=arma::chol(cor_cop_1).t();
  
  double cop_chol_det_1= arma::det(cor_cop_chol_1);
  
  arma::mat cor_cop_chol_2=arma::chol(cor_cop_2).t();
  
  double cop_chol_det_2= arma::det(cor_cop_chol_2);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  
  filterProbs(0,0)=p1t;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_cop_chol_1*H_eigen_inv*r.row(i).t()));
    double d_norm_3=arma::prod(dnorm_cpp(cor_cop_chol_2*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_cop_chol_1*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_3=pnorm_cpp(cor_cop_chol_2*H_eigen_inv*r.row(i).t());
    
    double llv_1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::det(H_eigen_inv)*d_norm_1);
    double llv_2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*VineCopula(p_norm_2,copula_par_1,type1)*d_norm_2*H_det*cop_chol_det_1);
    double llv_3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*VineCopula(p_norm_3,copula_par_2,type2)*d_norm_2*H_det*cop_chol_det_2);
    
    double llv_temp = llv_1+llv_2+llv_3;
    p1t=llv_1/llv_temp;
    p2t=llv_2/llv_temp;
    p3t=llv_3/llv_temp;
    
    llv+=log(llv_temp);
    filterProbs(i,0)=p1t;
    filterProbs(i,1)=p2t;
    filterProbs(i,2)=p3t;
    
  }
  return filterProbs;
}

// [[Rcpp::export]]
arma::mat FilterProbs_Normal_Copula_Copula_3_Asymm(const arma::vec& bekk, arma::vec signs, const arma::vec& theta, const arma::mat& r,arma::vec& type) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 3 * std::pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  
  int index = 0;
  
  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = bekk(index);
      index += 1;
    }
  }
  
  arma::mat A = arma::reshape(bekk.subvec(index, (index + std::pow(n, 2)) - 1 ).t(), n, n);
  arma::mat B = arma::reshape(bekk.subvec((index +  std::pow(n, 2)), index +  2*std::pow(n, 2)-1).t(), n, n);
  arma::mat G = arma::reshape(bekk.subvec((index +  2*std::pow(n, 2)), numb_of_vars-1).t(), n, n);
  double p11 = theta[0];
  double p12 = theta[1];
  double p13 =1 -p11 - p12;
  double p21 = theta[2];
  double p22 = theta[3];
  double p23 =1 -p21 - p22;
  double p32 = theta[4];
  double p33 = theta[5];
  double p31 =1 -p33 - p32;
  arma::vec copula_par = theta.subvec(6,11);
  
  arma::mat filterProbs = arma::zeros(NoOBs,3);
  //check for identification
  
  
  arma::vec copula_par_1 = theta.subvec(6,8);
  arma::vec copula_par_2 = theta.subvec(9,11);
  arma::vec type1 = type.subvec(0,2);
  arma::vec type2 = type.subvec(3,5);
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Bt  = B.t();
  arma::mat Gt  = G.t();
  
  //compute covariance matrix
  
 
  arma::mat cor_cop_1=cor_mat(copula_par_1,type1);
 
  
  
  arma::mat cor_cop_2=cor_mat(copula_par_2,type2);
  
  
 
  
  arma::mat cor_cop_chol_1=arma::chol(cor_cop_1).t();
  
  double cop_chol_det_1= arma::det(cor_cop_chol_1);
  
  arma::mat cor_cop_chol_2=arma::chol(cor_cop_2).t();
  
  double cop_chol_det_2= arma::det(cor_cop_chol_2);
  
  arma::mat H_eigen_inv=arma::inv(eigen_value_decomposition(H));
  
  
  
  double p1t = 1;
  double p2t = 0;
  double p3t = 0;
  filterProbs(0,0)=p1t;
  
  double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(H_eigen_inv*r.row(0).t());
  
  
  double llv = log(arma::as_scalar(arma::det(H_eigen_inv)*d_norm_1));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + indicatorFunction(r.row(i - 1), signs) * Bt * r.row(i - 1).t() * r.row(i - 1) * B + Gt * H * G;
    
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_cop_chol_1*H_eigen_inv*r.row(i).t()));
    double d_norm_3=arma::prod(dnorm_cpp(cor_cop_chol_2*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_2=pnorm_cpp(cor_cop_chol_1*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_3=pnorm_cpp(cor_cop_chol_2*H_eigen_inv*r.row(i).t());
    
    double llv_1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::det(H_eigen_inv)*d_norm_1);
    double llv_2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p32)*VineCopula(p_norm_2,copula_par_1,type1)*d_norm_2*H_det*cop_chol_det_1);
    double llv_3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*VineCopula(p_norm_3,copula_par_2,type2)*d_norm_2*H_det*cop_chol_det_2);
    
    double llv_temp = llv_1+llv_2+llv_3;
    p1t=llv_1/llv_temp;
    p2t=llv_2/llv_temp;
    p3t=llv_3/llv_temp;
    
    llv+=log(llv_temp);
    
    filterProbs(i,0)=p1t;
    filterProbs(i,1)=p2t;
    filterProbs(i,2)=p3t;
  }
  return filterProbs;
}



// 
// // [[Rcpp::export]]
// Rcpp::List random_grid_search_copula_3(const arma::vec& bekk,const arma::vec& theta,arma::vec& type,const arma::mat& r, const int nc) {
//   int n =r.n_cols;
//   int l=0;
//   
//   
//   arma::vec theta_mu;
//   theta_mu=theta;
//   arma::vec theta_optim=theta;
//   arma::vec theta_candidate=theta; 
//   double best_val = loglike_Normal_Copula_3(bekk,theta,r,type);
//   //set the seed
//   
//   
//   // Generating random values theta
//   while(l<(100/nc)){
//     for(int i=0; i<3;i++){
//       theta_candidate[i]= arma::randn()*0.04+theta_mu[i];
//     }
//     double llv = loglike_Clayton_Gumbel(bekk,theta_candidate,r);
//     if(llv > -1e25){
//       l++;
//     }
//     if(llv>best_val){
//       best_val=llv;
//       theta_optim=theta_candidate;
//       theta_mu=theta_optim;
//       
//     }
//   }
//   return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
//                             Rcpp::Named("best_val") = best_val);
//   
// }
// 
// // [[Rcpp::export]]
// Rcpp::List random_grid_search_copula_3_asymm(const arma::vec& bekk,arma::vec signs,const arma::vec& theta,arma::vec& type,const arma::mat& r, const int nc) {
//   int n =r.n_cols;
//   int l=0;
//   
//   
//   arma::vec theta_mu;
//   theta_mu=theta;
//   arma::vec theta_optim=theta;
//   arma::vec theta_candidate=theta; 
//   double best_val = loglike_Normal_Copula_3_asymm(bekk,signs,theta,r,type);
//   //set the seed
//   
//   
//   // Generating random values theta
//   while(l<(100/nc)){
//     for(int i=0; i<3;i++){
//       theta_candidate[i]= arma::randn()*0.04+theta_mu[i];
//     }
//     double llv = loglike_Clayton_Gumbel(bekk,theta_candidate,r);
//     if(llv > -1e25){
//       l++;
//     }
//     if(llv>best_val){
//       best_val=llv;
//       theta_optim=theta_candidate;
//       theta_mu=theta_optim;
//       
//     }
//   }
//   return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
//                             Rcpp::Named("best_val") = best_val);
//   
// }
// 

// [[Rcpp::export]]
Rcpp::List random_grid_search_LL_copula_3(arma::vec& theta,arma::mat& r, arma::vec& type, const int nc) {
  int n =r.n_cols;
  int l=0;
  int nOf_BEKK_Par = 2*pow(n,2)+n*(n+1)/2;
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  double best_val = LL_loglike_Copula_3(theta,r,type);
  //set the seed
  
  
  // Generating random values theta
  while(l<(100/nc)){
    for(int i=nOf_BEKK_Par; i<(nOf_BEKK_Par+3);i++){
      theta_candidate[i]= arma::randn()*0.01+theta_mu[i];
    }
    
    double llv = LL_loglike_Copula_3(theta_candidate,r,type);
    
    if(llv > -1e25){
      l++;
    }
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}


// [[Rcpp::export]]
Rcpp::List random_grid_search_LL_copula_3_asymm(arma::vec& theta, arma::vec& signs, arma::mat& r,arma::vec& type, int nc) {
  int n =r.n_cols;
  int l=0;
  int nOf_BEKK_Par = 3*pow(n,2)+n*(n+1)/2;
  
  arma::vec theta_mu;
  theta_mu=theta;
  arma::vec theta_optim=theta;
  arma::vec theta_candidate=theta; 
  double best_val = LL_loglike_Copula_3_asymm(theta,signs,r,type);
  //set the seed
  
  
  // Generating random values theta
  while(l<(100/nc)){
    for(int i=nOf_BEKK_Par; i<(nOf_BEKK_Par+3);i++){
      theta_candidate[i]= arma::randn()*0.01+theta_mu[i];
    }
    
    double llv = LL_loglike_Copula_3_asymm(theta_candidate,signs,r,type);
    
    if(llv > -1e25){
      l++;
    }
    if(llv>best_val){
      best_val=llv;
      theta_optim=theta_candidate;
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("theta_optim") = theta_optim,
                            Rcpp::Named("best_val") = best_val);
  
}


