#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::vec pnorm_cpp(arma::vec r){
  
  // Obtaining namespace of Matrix package
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("stats");
  // Picking up Matrix() function from Matrix package
  Rcpp::Function f = pkg["pnorm"];
  Rcpp::NumericVector res =f(r);
  // Executing Matrix( m, sparse = TRIE )
  return res;
}
// [[Rcpp::export]]
arma::vec dnorm_cpp(arma::vec r){
  
  // Obtaining namespace of stats package
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("stats");
  
  
  Rcpp::Function f = pkg["dnorm"];
  Rcpp::NumericVector res =f(r);
  
  // if(res(0)==0.0){
  //   res(0)=1e-324;
  // }
  // if(res(1)==0.0){
  //   res(1)=1e-324;
  // }
  return res;
}

// [[Rcpp::export]]
bool valid_bekk(arma::mat& C,arma::mat& A,arma::mat& G){
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

// [[Rcpp::export]]
double gumbelCDF(double u1, double u2, double theta){
  double res = std::exp(-std::pow(std::pow(-std::log(u1),theta) + std::pow(-std::log(u2),theta),1.0/theta));
  return res;
}

// [[Rcpp::export]]
double gumbelPDF(double u1,double u2, double theta){
  
  if(theta==1.0){
    return 1.0;
  }
  if(u1<1e-10){
    u1=1e-10;

  }
  if(u2<1e-10){
    u2=1e-10;

  }
  if(u1>1.0-1e-10){
    u1=1.0-1e-10;

  }
  if(u2>1.0-1e-10){
    u2=1.0-1e-10;

  }
 //double res=gumbelCDF(u1,u2,theta) *std::pow(u1*u2,-1.0)*std::pow(log(u1)*log(u2),theta-1.0)*(std::pow(std::pow(-log(u1),theta) +std::pow(-log(u2),theta),-2.0+(2.0/theta))) *(1+ (theta-1)*std::pow(std::pow(-log(u1),theta) +std::pow(-log(u2),theta),-1.0/theta));
 // double thetha1 = 1.0 / theta;
 //  double t1 = std::pow(-std::log(u1), theta) + std::pow(-std::log(u2), theta);
 //  double temp = -std::pow(t1, thetha1) + (2 * thetha1 - 2.0) * std::log(t1) +
 //    (theta - 1.0) * std::log(std::log(u1) * std::log(u2)) -
 //    std::log(u1 * u2) +   std::log(1+(theta - 1.0) * std::pow(t1, -thetha1));
 //  double res = std::exp(temp);
  
  double res=std::log(gumbelCDF(u1,u2,theta))-std::log(u1*u2)+(theta-1.0)*std::log(log(u1)*log(u2))+(-2.0+(2.0/theta))*std::log(std::pow(-log(u1),theta) +std::pow(-log(u2),theta)) +std::log(1+ (theta-1)*std::pow(std::pow(-log(u1),theta) +std::pow(-log(u2),theta),-1.0/theta));
  
   res= std::exp(res);

  
  if((std::isnan(res)|| res== 0.0 || res==R_PosInf) && u1>0.5 && u2>0.5 && u1==u2){
    double res =std::nextafter(R_PosInf,1.0);
    return res;
  }
  else if(res==R_PosInf && u1==u2){
    res = std::nextafter(R_PosInf,1.0);
    return res;
  }
  
  else if((std::isnan(res) || res==R_PosInf) && u1==u2 && u1<0.5 && u2<0.5){
    return std::nextafter(R_PosInf,1.0);
  } 
  else if(std::isnan(res) || res==R_PosInf || res==0.0){
    double res = std::nextafter(0.0,1.0);
    
    return res;
    
  }
  return res;
  
}
// [[Rcpp::export]]
double gumbelPDF_raw(double u1,double u2, double theta){
  if(u1<1e-10){
    u1=1e-10;
    
  }
  if(u2<1e-10){
    u2=1e-10;
    
  }
  if(u1>1.0-1e-10){
    u1=1.0-1e-10;
    
  }
  if(u2>1.0-1e-10){
    u2=1.0-1e-10;
    
  }
  if(theta==1){
    return 1;
  }
  double res= gumbelCDF(u1,u2,theta) *std::pow(u1*u2,-1.0)*std::pow(log(u1)*log(u2),theta-1.0)*(std::pow(std::pow(-log(u1),theta) +std::pow(-log(u2),theta),-2.0+(2.0/theta))) *(1+ (theta-1)*std::pow(std::pow(-log(u1),theta) +std::pow(-log(u2),theta),-1.0/theta));
if(std::isnan(res) ){
  res = std::nextafter(res,1.0);
}
return res;
  }
// [[Rcpp::export]]
double gumbelH1(double u1,double u2, double theta){
  
  
  if(u1<1e-10){
    u1=1e-10;
    
  }
  if(u2<1e-10){
    u2=1e-10;
    
  }
  if(u1>1.0-1e-10){
    u1=1.0-1e-10;
    
  }
  if(u2>1.0-1e-10){
    u2=1.0-1e-10;
    
  }
  
  double x= pow(-log(u1),theta) + pow(-log(u2),theta);
  double res=gumbelCDF(u1,u2,theta)*pow(-log(u1),-1+theta)*pow(x,-1+1/theta)/u1;
  //double res=std::log(gumbelCDF(u1,u2,theta))-std::log(u1*u2)+(theta-1.0)*std::log(log(u1)*log(u2))+(-2.0+(2.0/theta))*std::log(std::pow(-log(u1),theta) +std::pow(-log(u2),theta)) +std::log(1+ (theta-1)*std::pow(std::pow(-log(u1),theta) +std::pow(-log(u2),theta),-1.0/theta));
  // double thetha1 = 1.0 / theta;
  // double t1 = std::pow(-std::log(u1), theta) + std::pow(-std::log(u2), theta);
  // double res = -std::pow(t1, thetha1) + (2 * thetha1 - 2.0) * std::log(t1) +
  //    (theta - 1.0) * std::log(std::log(u1) * std::log(u2)) -
  //    std::log(u1 * u2) +
  //    std::log((theta) * std::pow(t1, -thetha1));
  //res= std::exp(res);
  if((std::isnan(res)|| res== 0.0 ||res==R_PosInf) && u1>0.5 && u2>0.5 && u1==u2){
    double res =std::nextafter(res,0.0);
    return res;
  }
  
  return res;
  
  
}
// [[Rcpp::export]]
double gumbelH2(double u1,double u2, double theta){
  
  
  if(u1<1e-10){
    u1=1e-10;
    
  }
  if(u2<1e-10){
    u2=1e-10;
    
  }
  if(u1>1.0-1e-10){
    u1=1.0-1e-10;
    
  }
  if(u2>1.0-1e-10){
    u2=1.0-1e-10;
    
  }
  double x= pow(-log(u1),theta) + pow(-log(u2),theta);
  double res=gumbelCDF(u1,u2,theta)*pow(-log(u2),-1+theta)*pow(x,-1+1/theta)/u2;
  //double res=std::log(gumbelCDF(u1,u2,theta))-std::log(u1*u2)+(theta-1.0)*std::log(log(u1)*log(u2))+(-2.0+(2.0/theta))*std::log(std::pow(-log(u1),theta) +std::pow(-log(u2),theta)) +std::log(1+ (theta-1)*std::pow(std::pow(-log(u1),theta) +std::pow(-log(u2),theta),-1.0/theta));
  // double thetha1 = 1.0 / theta;
  // double t1 = std::pow(-std::log(u1), theta) + std::pow(-std::log(u2), theta);
  // double res = -std::pow(t1, thetha1) + (2 * thetha1 - 2.0) * std::log(t1) +
  //    (theta - 1.0) * std::log(std::log(u1) * std::log(u2)) -
  //    std::log(u1 * u2) +
  //    std::log((theta) * std::pow(t1, -thetha1));
  //res= std::exp(res);
  if((std::isnan(res)|| res== 0.0 ||res==R_PosInf) && u1>0.5 && u2>0.5 && u1==u2){
    double res =std::nextafter(res,1.0);
    return res;
  }
  return res;
  
  
}


// [[Rcpp::export]]
arma::mat cor_Gumbel(double theta){
  arma::mat cor_mat = arma::zeros(2,2);
  cor_mat(0,0)=1;
  cor_mat(1,1)=1;
  cor_mat(1,0)=sin(M_PI/2 * (1-1/theta));
  cor_mat(0,1)=sin(M_PI/2 * (1-1/theta));
  return cor_mat;
}


// [[Rcpp::export]]
double claytonCDF(double u1, double u2, double theta){
  double res = std::pow(std::pow(u1,-theta)+std::pow(u2,-theta)-1,-1/theta);
  return res;
}
// [[Rcpp::export]]
double claytonPDF(double u1,double u2, double theta){
  if(u1<1e-10){
    u1=1e-10;
    
  }
  if(u2<1e-10){
    u2=1e-10;
    
  }
  if(u1>1.0-1e-10){
    u1=1.0-1e-10;
    
  }
  if(u2>1.0-1e-10){
    u2=1.0-1e-10;
    
  }
  if(theta==0){
    return 1;
  }
  double temp = std::log(theta+1) - (1.0 + theta) * std::log(u1 * u2);
  temp = temp - (2.0 + 1.0 / (theta)) *  std::log(std::pow(u1, -theta) + std::pow(u2, -theta) - 1.0);
  
  double res= std::exp(temp);
  
  if((std::isnan(res) || res==0.0 || res==R_PosInf ||  res==R_NegInf) && u1<0.5 && u2<0.5 && u1==u2){
    double res = std::nextafter(R_PosInf,1.0);
    return res;
  }
  else if(std::isnan(res)){
    double res = std::nextafter(0.0,1.0);
    return res;
  }
  
  if(res==R_PosInf || res==R_NegInf ){
    res=std::nextafter(0.0,1.0);
    return res;
  }
  if(res==0.0){
    res=std::nextafter(0.0,1.0);
  }
  
  return res;
  
}

// [[Rcpp::export]]
double claytonPDF_raw(double u1,double u2, double theta){
  
  if(u1<1e-10){
    u1=1e-10;
    
  }
  if(u2<1e-10){
    u2=1e-10;
    
  }
  if(u1>1.0-1e-10){
    u1=1.0-1e-10;
    
  }
  if(u2>1.0-1e-10){
    u2=1.0-1e-10;
    
  }
  if(theta==0){
    return 1;
  }

  double temp = std::log(theta+1) - (1.0 + theta) * std::log(u1 * u2);
  temp = temp - (2.0 + 1.0 / (theta)) *  std::log(std::pow(u1, -theta) + std::pow(u2, -theta) - 1.0);
  
  double res= std::exp(temp);
  
  if(std::isnan(res)){
    res = std::nextafter(res,1.0);
  }
 
  return res;
  
}

// [[Rcpp::export]]
double claytonH1(double u1,double u2, double theta){
  if(u1<1e-10){
    u1=1e-10;
    
  }
  if(u2<1e-10){
    u2=1e-10;
    
  }
  if(u1>1.0-1e-10){
    u1=1.0-1e-10;
    
  }
  if(u2>1.0-1e-10){
    u2=1.0-1e-10;
    
  }
  
  
  double res = std::pow(std::pow(u1,-theta)+std::pow(u2,-theta)-1,-1/theta-1)*pow(u1,-theta-1);
  
  
  if(std::isnan(res)|| res== 0.0 ||res==R_PosInf){
    double res =std::nextafter(res,0.0);
    return res;
  }
  
  return res;
  
}


// [[Rcpp::export]]
double claytonH2(double u1,double u2, double theta){
  if(u1<1e-10){
    u1=1e-10;
    
  }
  if(u2<1e-10){
    u2=1e-10;
    
  }
  if(u1>1.0-1e-10){
    u1=1.0-1e-10;
    
  }
  if(u2>1.0-1e-10){
    u2=1.0-1e-10;
    
  }
 
  double res = std::pow(std::pow(u1,-theta)+std::pow(u2,-theta)-1,-1/theta-1)*pow(u2,-theta-1);
  
  if(std::isnan(res)|| res== 0.0 ||res==R_PosInf){
    double res =std::nextafter(res,0.0);
    return res;
  }
  return res;
  
}

// [[Rcpp::export]]
double claytonPDF_log(double u1,double u2, double theta){
  
  double res= std::log(1.0+theta)+std::log(std::pow(u1,-theta)+std::pow(u2,-theta)-1.0)*(-1.0/theta -2.0) -std::log(u1*u2)*(theta+1);
  
  if((std::isnan(res) || res==R_NegInf) && u1<0.5 && u2<0.5 && u1==u2){
    double res = std::nextafter(R_PosInf,0.0);
    return res;
  }
  else if(std::isnan(res)){
    double res = std::nextafter(res,1.0);
    return res;
  }
  
  if(res==R_PosInf){
    res=std::nextafter(res,0.0);
    return res;
  }
  if(res==R_NegInf){
    res = std::nextafter(R_NegInf,0.0);
    return res;
  }
  
  return res;
  
}

// [[Rcpp::export]]
arma::mat cor_Clayton(double theta){
  arma::mat cor_mat = arma::zeros(2,2);
  cor_mat(0,0)=1;
  cor_mat(1,1)=1;
  cor_mat(1,0)=sin(M_PI/2 * (theta/(theta+2)));
  cor_mat(0,1)=sin(M_PI/2 * (theta/(theta+2)));
  return cor_mat;
}



// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::export]]
arma::mat eigen_value_decomposition(arma::mat& A){
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym( eigval, eigvec, A );
  
  
  arma::mat diag_mat_eigval = arma::diagmat(sqrt(eigval));
  return eigvec*diag_mat_eigval*eigvec.t();
  
}

// [[Rcpp::export]]
arma::mat BiCopSim_cpp(int NoOSim, int family, double par){
  
  Rcpp::Environment pkg =  Rcpp::Environment::namespace_env("VineCopula");
  
  
  Rcpp::Function f = pkg["BiCopSim"];
  Rcpp::NumericMatrix res= f(NoOSim, family, par);
  
  
  return  Rcpp::as<arma::mat>(res);
}

// [[Rcpp::export]]
arma::mat rbicop_cpp(int NoOSim, 	const char* family , int rot , double par){
  
  Rcpp::Environment pkg =  Rcpp::Environment::namespace_env("rvinecopulib");
  
  
  Rcpp::Function f = pkg["rbicop"];
  Rcpp::NumericMatrix res= f(NoOSim, family,rot, par);
  
  
  return  Rcpp::as<arma::mat>(res);
}

// [[Rcpp::export]]
double dmvnorm_cpp(arma::mat r, arma::vec& mean, arma::mat& sigma){
  
  // Obtaining namespace of Matrix package
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("mvtnorm");
  
  // Picking up Matrix() function from Matrix package
  Rcpp::Function f = pkg["dmvnorm"];
  double res =Rcpp::as<double>(f(r, mean,  sigma));
  // Executing Matrix( m, sparse = TRIE )
  return res;
}
//rotation matrices

// choose 2 int from N (as index)
arma::imat choose2_fast(int &N)
{
  arma::imat Out(N * (N - 1) / 2, 2);
  int c = 0;
  for (int i = 0; i < (N - 1); i++)
  {
    for (int j = i + 1; j < N; j++)
    {
      Out(c, 0) = i;
      Out(c, 1) = j;
      c++;
    }
  }
  return Out;
}

// Givens rotation matrix
arma::mat rot_mat(arma::vec &thetas, int K)
{
  arma::mat Out = arma::eye(K, K);
  arma::imat Cmat = choose2_fast(K);
  for (int i = 0; i < K * (K - 1) / 2; i++)
  {
    arma::mat temp = arma::eye(K, K);
    temp(Cmat(i, 0), Cmat(i, 0)) = cos(thetas(i));
    temp(Cmat(i, 1), Cmat(i, 1)) = cos(thetas(i));
    temp(Cmat(i, 0), Cmat(i, 1)) = -sin(thetas(i));
    temp(Cmat(i, 1), Cmat(i, 0)) = sin(thetas(i));
    
    Out = Out * temp;
  }
  return Out.t();
}  



// [[Rcpp::export]]
double BiCopPDF_cpp(double returns1, double returns2 , int family, double par){
  
  // Obtaining namespace of Matrix package
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("VineCopula");
  
  
  // Picking up Matrix() function from Matrix package
  Rcpp::Function f = pkg["BiCopPDF"];
  Rcpp::NumericVector res = f(returns1, returns2, family, par);
  // Executing Matrix( m, sparse = TRIE )
  return Rcpp::as<double>(res);
}
// [[Rcpp::export]]
double dbicop_cpp(arma::mat returns, 	const char* family , int rot, double par){
  
  // Obtaining namespace of Matrix package
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("rvinecopulib");
  
  
  // Picking up Matrix() function from Matrix package
  Rcpp::Function f = pkg["dbicop"];
  Rcpp::NumericVector res = f(returns, family,rot, par);
  // Executing Matrix( m, sparse = TRIE )
  return Rcpp::as<double>(res);
}

// [[Rcpp::export]]
double copulaPDF(double u1,double u2, double theta, int type){
  if(type==3){
    return claytonPDF(u1,u2,theta);
  }
  if(type==13){
    return claytonPDF(1-u1,1-u2,theta);
  }
  if(type==4){
    return gumbelPDF(u1,u2,theta); 
  }
  if(type==14){
    return gumbelPDF(1-u1,1-u2,theta); 
  }
}

// [[Rcpp::export]]
double copulaPDF_raw(double u1,double u2, double theta, int type){
  if(type==3){
    return claytonPDF_raw(u1,u2,theta);
  }
  if(type==13){
    return claytonPDF_raw(1-u1,1-u2,theta);
  }
  if(type==4){
    return gumbelPDF_raw(u1,u2,theta); 
  }
  if(type==14){
    return gumbelPDF_raw(1-u1,1-u2,theta); 
  }
}

// [[Rcpp::export]]
double copulaH1(double u1,double u2, double theta, int type){
  if(type==3){
    return claytonH1(u1,u2,theta);
  }
  if(type==13){
    return 1-claytonH1(1-u1,1-u2,theta);
  }
  if(type==4){
    return gumbelH1(u1,u2,theta); 
  }
  if(type==14){
    return 1-gumbelH1(1-u1,1-u2,theta); 
  }
}

// [[Rcpp::export]]
double copulaH2(double u1,double u2, double theta, int type){
  if(type==3){
    return claytonH2(u1,u2,theta);
  }
  if(type==13){
    return 1-claytonH2(1-u1,1-u2,theta);
  }
  if(type==4){
    return gumbelH2(u1,u2,theta); 
  }
  if(type==14){
    return 1-gumbelH2(1-u1,1-u2,theta); 
  }
}

// [[Rcpp::export]]
bool valid_copula(arma::vec theta, arma::vec type){
  for (int i=0; i< type.n_rows; i++){
  if(type[i]==3 || type[i]==13){
    if(theta[i]<0 || theta[i] >= 30.7){
      return false;
    }
  }
  if(type[i]==4 || type[i]==14){
    if(theta[i]<1 ||  theta[i] > 31.5){
      return false;
    }
  }
  
  }
  return true;
  
  
}
// [[Rcpp::export]]
double VineCopula2(arma::vec u, arma::vec copula_par, arma::vec type){
  double res=  copulaPDF(copulaH2(u[0],u[2],copula_par[1],type[1]),copulaH2(u[1],u[2],copula_par[2],type[2]),copula_par[0],type[0])*copulaPDF(u[0],u[2],copula_par[1],type[1])*copulaPDF(u[1],u[2],copula_par[2],type[2]);
if(res == R_PosInf ||  std::isnan(res) || res == 0.0){
  return nextafter(res, 1.0);
}
return res;
}

// [[Rcpp::export]]
double VineH2(double u1,double u2,double u3, arma::vec copula_par, arma::vec type){
  double res=  copulaH2(copulaH2(u1,u3,copula_par[1],type[1]),copulaH2(u2,u3,copula_par[2],type[2]),copula_par[0],type[0])*copulaPDF(u2,u3,copula_par[2],type[2]);
  if(res == R_PosInf ){
    return nextafter(res, 1.0);
  }
  // if(std::isnan(res)){
  //   res=0.0;
  // }
  
  return res;
}
// [[Rcpp::export]]
double VineCopula(arma::vec u, arma::vec copula_par, arma::vec type){
  double res = copulaPDF(copulaH2(u[0],u[2],copula_par[1],type[1]),copulaH2(u[1],u[2],copula_par[2],type[2]),copula_par[0],type[0])*copulaPDF(u[0],u[2],copula_par[1],type[1])*copulaPDF(u[1],u[2],copula_par[2],type[2]);
  if(res == R_PosInf){
    return std::nextafter(R_PosInf,1.0);
  }
  if(std::isnan(res) || res==0.0){
    return std::nextafter(res,1.0);
  }
  
  return res;
}



    
