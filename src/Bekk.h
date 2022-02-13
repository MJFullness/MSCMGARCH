
#ifndef Bekk_H
#define Bekk_H





inline bool valid_bekk(arma::mat& C,arma::mat& A,arma::mat& G){
  int n =C.n_cols;
  arma::mat prod = kron(A,A)+kron(G,G);
  
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
  if(A(0,0)<=0 || G(0,0)<=0) {
    return false;
  }
  else{
    return true;
  }
}




inline int indicatorFunction(arma::mat r, arma::mat signs){
  r = r.t();
  
  int indicator = 1;
  int n = r.n_rows;
  for (int i = 0; i<n; i++){
    if(arma::as_scalar(signs.row(i)) * arma::as_scalar(r.row(i)) < 0){
      indicator = 0;
    }
  }
  return indicator;
}


inline bool valid_asymm_bekk(arma::mat& C,arma::mat& A, arma::mat& B ,arma::mat& G, arma::mat r, arma::mat signs){
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



#endif




