gumbelCDF<-function(u,v,theta){
return(exp(-((-log(u))^theta+(-log(v))^theta)^((1/theta))))
}

gumbelCDF_derivative_u<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
  res=gumbelCDF(u,v,theta)*(-log(u))^(-1+theta)*x^(-1+1/theta)/u
  return(res)
}
gumbelCDF_derivative_v<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
  res=(gumbelCDF(u,v,theta)*x^(-1+1/theta)*(-log(v))^(-1+theta))/v
  return(res)
}

gumbelCDF_derivative_u_theta<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
  res=(1/(theta^2*u*log(u)))*gumbelCDF(u,v,theta)*(-log(u))^theta*x^(-2+1/theta)*(-theta*(-(-log(u))^theta*(-1+x^(1/theta))+theta*(-log(v))^theta)*log(-log(u))-(-1+x^(1/theta))*x*log((-log(u))^theta+(-log(v))^theta)+theta*(-1+theta+x^(1/theta))*(-log(v))^theta*log(-log(v)))
  return(res)
  }

gumbelCDF_derivative_v_theta<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
 res=(1/(theta^2*v*log(v)))*gumbelCDF(u,v,theta)*x^(-2+1/theta)*(-log(v))^theta*(theta*(-log(u))^theta*(-1+theta+x^(1/theta))*log(-log(u))-(-1+x^(1/theta))*x*log((-log(u))^theta+(-log(v))^theta)-theta*(theta*(-log(u))^theta-(-1+x^(1/theta))*(-log(v))^theta)*log(-log(v)))
 return(res)
 }

gumbelCDF_derivative_u_u<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
 res=(gumbelCDF(u,v,theta)*(-log(u))^(-2+theta)*x^(-2+1/theta)*((-log(u))^theta*(log(u)+x^(1/theta))+(1-theta+log(u))*(-log(v))^theta))/u^2
 return(res)
 }
gumbelCDF_derivative_v_v<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
  res=(gumbelCDF(u,v,theta)*x^(-2+1/theta)*(-log(v))^(-2+theta)*((-log(u))^theta*(1-theta+log(v))+(-log(v))^theta*(x^(1/theta)+log(v))))/v^2
  return(res)
  }

gumbelCDF_derivative_u_u_u<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
  res=(1/(u^3*log(u)^3))*gumbelCDF(u,v,theta)*(-log(u))^theta*x^(-3+1/theta)*(-(-log(u))^(2*theta)*(log(u)+x^(1/theta))*(2*log(u)+x^(1/theta))+(-log(u))^theta*(-1+theta^2-3*log(u)+3*theta*log(u)-4*log(u)^2+3*(-1+theta-log(u))*x^(1/theta))*(-log(v))^theta-((-2+theta)*(-1+theta)+log(u)*(3-3*theta+2*log(u)))*(-log(v))^(2*theta))
  return(res)
}

gumbelCDF_derivative_v_v_v<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
  res=(1/(v^3*log(v)^3))*gumbelCDF(u,v,theta)*x^(-3+1/theta)*(-log(v))^theta*(-(-log(v))^(2*theta)*(x^(1/theta)+log(v))*(x^(1/theta)+2*log(v))+(-log(u))^theta*(-log(v))^theta*(-1+theta^2+3*x^(1/theta)*(-1+theta-log(v))-3*log(v)+3*theta*log(v)-4*log(v)^2)-(-log(u))^(2*theta)*((-2+theta)*(-1+theta)+log(v)*(3-3*theta+2*log(v))))
  return(res)
}

gumbelCDF_derivative_u_theta_theta<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
res=(1/(theta^4*u*log(u)))*gumbelCDF(u,v,theta)*(-log(u))^theta*x^(-3+1/theta)*(theta^2*(-(-log(u))^(2*theta)*(1-3*x^(1/theta)+x^(2/theta))+theta*((-log(u))^theta*(-3+theta+3*x^(1/theta))-theta*(-log(v))^theta)*(-log(v))^theta)*log(-log(u))^2-(1-3*x^(1/theta)+x^(2/theta))*x^2*log((-log(u))^theta+(-log(v))^theta)^2+theta^2*(-log(v))^theta*log(-log(v))*(-2*(-1+x^(1/theta))*x+(theta*(-log(u))^theta*(-1+theta+x^(1/theta))-((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))*(-log(v))^theta)*log(-log(v)))+2*theta*x*log((-log(u))^theta+(-log(v))^theta)*((-1+x^(1/theta))*x+(1-theta+(-3+theta)*x^(1/theta)+x^(2/theta))*(-log(v))^theta*log(-log(v)))+2*theta*log(-log(u))*(-theta*(-log(u))^theta*(-1+x^(1/theta))*x+x*((-log(u))^theta*(1-3*x^(1/theta)+x^(2/theta))-theta*(-1+x^(1/theta))*(-log(v))^theta)*log((-log(u))^theta+(-log(v))^theta)+theta*(-(-log(u))^theta*((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))+theta*(-1+theta+x^(1/theta))*(-log(v))^theta)*(-log(v))^theta*log(-log(v))))
return(res)
}
gumbelCDF_derivative_v_theta_theta<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
res=(1/(theta^4*v*log(v)))*gumbelCDF(u,v,theta)*x^(-3+1/theta)*(-log(v))^theta*(theta^2*(-log(u))^theta*(-(-log(u))^theta*((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))+theta*(-1+theta+x^(1/theta))*(-log(v))^theta)*log(-log(u))^2-(1-3*x^(1/theta)+x^(2/theta))*x^2*log((-log(u))^theta+(-log(v))^theta)^2+2*theta*(-log(u))^theta*log(-log(u))*(-theta*(-1+x^(1/theta))*x+(1-theta+(-3+theta)*x^(1/theta)+x^(2/theta))*x*log((-log(u))^theta+(-log(v))^theta)+theta*(theta*(-log(u))^theta*(-1+theta+x^(1/theta))-((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))*(-log(v))^theta)*log(-log(v)))+theta^2*log(-log(v))*(-2*(-1+x^(1/theta))*x*(-log(v))^theta+(-theta^2*(-log(u))^(2*theta)+theta*(-log(u))^theta*(-3+theta+3*x^(1/theta))*(-log(v))^theta-(1-3*x^(1/theta)+x^(2/theta))*(-log(v))^(2*theta))*log(-log(v)))-2*theta*x*log((-log(u))^theta+(-log(v))^theta)*((-log(u))^theta*(-1+x^(1/theta))*(-1+theta*log(-log(v)))-(-log(v))^theta*(-1+x^(1/theta)*(1-3*log(-log(v)))+log(-log(v))+x^(2/theta)*log(-log(v)))))
return(res)
}

gumbelPDF<-function(u,v,theta){
return((gumbelCDF(u,v,theta)*(-log(u))^(-1+theta)* (-1+theta+x^(1/theta)) *x^(-2+1/theta)* (-log(v))^(-1+theta))/(u*v))

}

gumbelPDF_derivative_u<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta

  res=-(1/(u^2*v))*gumbelCDF(u,v,theta)*(-log(u))^(-2 +theta)* x^(-3 + 1/theta)* (-(-log(u))^theta * ((-1 + theta) * (theta +log(u)) + (-2 + 2 * theta +log(u)) * x^(1/theta) + x^(2/theta)) + (-1 + theta -log(u)) * (-1 + theta + x^(1/theta)) * (-log(v))^theta)* (-log(v))^(-1 + theta)
  return(res)
}

gumbelPDF_derivative_v<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta

  res=-(1/(u*v^2))*gumbelCDF(u,v,theta)*(-log(u))^(-1+theta) * x^(-3+1/theta) * (-log(v))^(-2+theta)* ((-log(u))^theta * (-1+theta+x^(1/theta)) * (-1+theta-log(v))-(-log(v))^theta * (x^(2/theta)+(-1+theta) * (theta+log(v))+(x)^(1/theta)*(-2+2* theta+log(v))))
  return(res)
}


gumbelPDF_derivative_theta<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
  res=(1/(theta^2 * u * v))*gumbelCDF(u,v,theta)* (-log(u))^(-1+theta) * x^(-2+1/theta) * (-log(v))^(-1+theta) * (theta^2+(1/x)*theta*(-(-log(u))^theta*((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))+theta*(-1+theta+x^(1/theta)) * (-log(v))^theta) *log(-log(u))+(1-theta+(-3+theta)*x^(1/theta)+x^(2/theta))*logx+(1/x)*theta*(theta*(-log(u))^theta*(-1+theta+x^(1/theta))-((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))*(-log(v))^theta)*log(-log(v)))
  return(res)
}

gumbelPDF_derivative_theta_theta<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
  res=(1/(theta^4*u*v)) *gumbelCDF(u,v,theta)* (-log(u))^(-1+theta)*x^(-4+1/theta)*(-log(v))^(-1+theta)*(theta^2*((-log(u))^(2*theta)*((-1+theta)^3+(7+3*(-3+theta)*theta)*x^(1/theta)+3*(-2+theta)*x^(2/theta)+x^(3/theta))-theta*(-log(u))^theta*(3-7*theta+4*theta^2+(-9+7*theta)*x^(1/theta)+3*x^(2/theta))*(-log(v))^theta+theta^2*(-1+theta+x^(1/theta))*(-log(v))^(2*theta))*log(-log(u))^2+(-1+theta+(7-3*theta)*x^(1/theta)+(-6+theta)*x^(2/theta)+x^(3/theta))*x^2*logx^2+2*theta*x*logx*(-((1-3*x^(1/theta)+x^(2/theta))*x)+(theta*(-log(u))^theta*(1-theta+(-3+theta)*x^(1/theta)+x^(2/theta))-(-(-1+theta)^2+(7+(-6+theta)*theta)*x^(1/theta)+2*(-3+theta)*x^(2/theta)+x^(3/theta))*(-log(v))^theta)*log(-log(v)))+theta^2*log(-log(v))*(theta^2*(-log(u))^(2*theta)*(2+(-1+theta+x^(1/theta))*log(-log(v)))+(-log(u))^theta*(-log(v))^theta*(2-3*theta*log(-log(v))+7*theta^2*log(-log(v))-4*theta^3*log(-log(v))+x^(2/theta)*(2-3*theta*log(-log(v)))+x^(1/theta)*(-6+(9-7*theta)*theta*log(-log(v))))+(-log(v))^(2*theta)*(2-2*theta^2+(-1+theta)^3*log(-log(v))+x^(3/theta)*log(-log(v))+x^(2/theta)*(2+3*(-2+theta)*log(-log(v)))+x^(1/theta)*(-6+(7+3*(-3+theta)*theta)*log(-log(v)))))-2*theta*log(-log(u))*(-(x*(-(-log(u))^theta*(-(-1+theta)^2+(7+(-6+theta)*theta)*x^(1/theta)+2*(-3+theta)*x^(2/theta)+x^(3/theta))+theta*(1-theta+(-3+theta)*x^(1/theta)+x^(2/theta))*(-log(v))^theta)*logx)+theta*(theta*(-log(v))^(2*theta)*(-theta+((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))*log(-log(v)))+(-log(u))^(2*theta)*(-1+theta^2+(-1+theta)^2*theta*log(-log(v))+x^(2/theta)*(-1+theta*log(-log(v)))+x^(1/theta)*(3+theta*(-3+2*theta)*log(-log(v))))-(-log(u))^theta*(-log(v))^theta*(1-log(-log(v))+4*theta*log(-log(v))-7*theta^2*log(-log(v))+4*theta^3*log(-log(v))+x^(3/theta)*log(-log(v))+x^(2/theta)*(1+(-6+4*theta)*log(-log(v)))+x^(1/theta)*(-3+(7+theta*(-12+7*theta))*log(-log(v)))))))
return(res)
}

gumbelPDF_derivative_theta_u<-function(u,v,theta){
  x<-(-log(u))^theta + (-log(v))^theta
  res=(1/(theta^4*u*v))*gumbelCDF(u,v,theta)*(-log(u))^(-1+theta)*x^(-4+1/theta)*(-log(v))^(-1+theta)*(theta^2*((-log(u))^(2*theta)*((-1+theta)^3+(7+3*(-3+theta)*theta)*x^(1/theta)+3*(-2+theta)*x^(2/theta)+x^(3/theta))-theta*(-log(u))^theta*(3-7*theta+4*theta^2+(-9+7*theta)*x^(1/theta)+3*x^(2/theta))*(-log(v))^theta+theta^2*(-1+theta+x^(1/theta))*(-log(v))^(2*theta))*log(-log(u))^2+(-1+theta+(7-3*theta)*x^(1/theta)+(-6+theta)*x^(2/theta)+x^(3/theta))*x^2*logx^2+2*theta*x*logx*(-((1-3*x^(1/theta)+x^(2/theta))*x)+(theta*(-log(u))^theta*(1-theta+(-3+theta)*x^(1/theta)+x^(2/theta))-(-(-1+theta)^2+(7+(-6+theta)*theta)*x^(1/theta)+2*(-3+theta)*x^(2/theta)+x^(3/theta))*(-log(v))^theta)*log(-log(v)))+theta^2*log(-log(v))*(theta^2*(-log(u))^(2*theta)*(2+(-1+theta+x^(1/theta))*log(-log(v)))+(-log(u))^theta*(-log(v))^theta*(2-3*theta*log(-log(v))+7*theta^2*log(-log(v))-4*theta^3*log(-log(v))+x^(2/theta)*(2-3*theta*log(-log(v)))+x^(1/theta)*(-6+(9-7*theta)*theta*log(-log(v))))+(-log(v))^(2*theta)*(2-2*theta^2+(-1+theta)^3*log(-log(v))+x^(3/theta)*log(-log(v))+x^(2/theta)*(2+3*(-2+theta)*log(-log(v)))+x^(1/theta)*(-6+(7+3*(-3+theta)*theta)*log(-log(v)))))-2*theta*log(-log(u))*(-(x*(-(-log(u))^theta*(-(-1+theta)^2+(7+(-6+theta)*theta)*x^(1/theta)+2*(-3+theta)*x^(2/theta)+x^(3/theta))+theta*(1-theta+(-3+theta)*x^(1/theta)+x^(2/theta))*(-log(v))^theta)*logx)+theta*(theta*(-log(v))^(2*theta)*(-theta+((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))*log(-log(v)))+(-log(u))^(2*theta)*(-1+theta^2+(-1+theta)^2*theta*log(-log(v))+x^(2/theta)*(-1+theta*log(-log(v)))+x^(1/theta)*(3+theta*(-3+2*theta)*log(-log(v))))-(-log(u))^theta*(-log(v))^theta*(1-log(-log(v))+4*theta*log(-log(v))-7*theta^2*log(-log(v))+4*theta^3*log(-log(v))+x^(3/theta)*log(-log(v))+x^(2/theta)*(1+(-6+4*theta)*log(-log(v)))+x^(1/theta)*(-3+(7+theta*(-12+7*theta))*log(-log(v)))))))
}

gumbelPDF_derivative_theta_v<-function(u,v,theta){
x<-(-log(u))^theta + (-log(v))^theta
  res=-(1/(theta^2*u*v^2))*gumbelCDF(u,v,theta)*(-log(u))^(-1+theta)*x^(-4+1/theta)*(-log(v))^(-2+theta)*(theta*(-((-1+theta)*(-log(u))^(2*theta)*((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta)))+(-log(u))^theta*((-1+theta)*theta*(-3+4*theta)+(-1+theta)*(-4+7*theta)*x^(1/theta)+(-5+4*theta)*x^(2/theta)+x^(3/theta))*(-log(v))^theta-theta*((-1+theta)*theta+2*(-1+theta)*x^(1/theta)+x^(2/theta))*(-log(v))^(2*theta))*log(-log(u))-(-log(u))^theta*(-log(v))^theta*((1-theta-(-1+theta)*x^(1/theta)+(-4+theta)*x^(2/theta)+x^(3/theta))*logx+theta*(theta*(1+x^(1/theta))+(-1+4*theta-7*theta^2+4*theta^3+(-1+theta)*(-3+7*theta)*x^(1/theta)+(-1+3*theta)*x^(2/theta))*log(-log(v))))+(-log(u))^(2*theta)*((-1+theta)*(1-theta+(-3+theta)*x^(1/theta)+x^(2/theta))*logx+theta^2*(x^(1/theta)*(1+(-1+theta)*log(-log(v)))+(-1+theta)*(2+(-1+theta)*log(-log(v)))))+(-log(v))^(2*theta)*(x^(3/theta)*(-logx+theta*log(-log(v)))+x^(2/theta)*((5-2*theta)*logx+theta*(-5+3*theta)*log(-log(v)))+theta*(theta-2*theta^2+(-1+theta)*(logx+(-1+theta)*theta*log(-log(v))))+x^(1/theta)*(-2*theta^2+(-1+theta)*(-((-4+theta)*logx)+theta*(-4+3*theta)*log(-log(v)))))+x*log(v)*(theta*((-log(u))^theta*((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))-theta*(-1+theta+x^(1/theta))*(-log(v))^theta)*log(-log(u))-(1-theta+(-3+theta)*x^(1/theta)+x^(2/theta))*x*logx+theta*(-theta*(-log(u))^theta*(1+(-1+theta+x^(1/theta))*log(-log(v)))+(-log(v))^theta*(-theta+((-1+theta)^2+(-3+2*theta)*x^(1/theta)+x^(2/theta))*log(-log(v))))))

}

gumbelPDF_derivative_u_u<-function(u,v,theta){
x<-(-log(u))^theta + (-log(v))^theta
  res=(1/(u^3*v))*gumbelCDF(u,v,theta)*(-log(u))^(-3+theta)*x^(-4+1/theta)*((-log(u))^(2*theta)*((-1+theta)*(theta+theta^2+3*theta*log(u)+2*log(u)^2)+(3*(-1+theta)*theta+2*log(u)*(-3+3*theta+log(u)))*x^(1/theta)+3*(-1+theta+log(u))*x^(2/theta)+x^(3/theta))+(-log(u))^theta*(-((-1+theta)*((-1+theta)*(1+4*theta)-log(u)*(3+4*log(u))))+(-5+12*theta-7*theta^2+3*(-1+theta)*log(u)+4*log(u)^2)*x^(1/theta)+3*(1-theta+log(u))*x^(2/theta))*(-log(v))^theta+((-2+theta)*(-1+theta)+log(u)*(3-3*theta+2*log(u)))*(-1+theta+x^(1/theta))*(-log(v))^(2*theta))*(-log(v))^(-1+theta)
return(res)
}

gumbelPDF_derivative_v_v<-function(u,v,theta){
x<-(-log(u))^theta + (-log(v))^theta
  res=(1/(u*v^3))*gumbelCDF(u,v,theta)*(-log(u))^(-1+theta)*x^(-4+1/theta)*(-log(v))^(-3+theta)*((-log(u))^(2*theta)*(-1+theta+x^(1/theta))*((-2+theta)*(-1+theta)+log(v)*(3-3*theta+2*log(v)))+(-log(v))^(2*theta)*(x^(3/theta)+3*x^(2/theta)*(-1+theta+log(v))+(-1+theta)*(theta+theta^2+3*theta*log(v)+2*log(v)^2)+x^(1/theta)*(3*(-1+theta)*theta+2*log(v)*(-3+3*theta+log(v))))+(-log(u))^theta*(-log(v))^theta*(-3*x^(2/theta)*(-1+theta-log(v))+x^(1/theta)*(-5+12*theta-7*theta^2+3*(-1+theta)*log(v)+4*log(v)^2)-(-1+theta)*((-1+theta)*(1+4*theta)-log(v)*(3+4*log(v)))))
  return(res)
}

gumbelPDF_derivative_u_v<-function(u,v,theta){
 x<-(-log(u))^theta + (-log(v))^theta
  res=(1/(u^2*v^2))*gumbelCDF(u,v,theta)*(-log(u))^(-2+theta)*x^(-4+1/theta)*(-log(v))^(-2+theta)*(-(-log(u))^(2*theta)*((-1+theta)*(theta+log(u))+(-2+2*theta+log(u))*x^(1/theta)+x^(2/theta))*(-1+theta-log(v))+(1-theta+log(u))*(-log(v))^(2*theta)*(x^(2/theta)+(-1+theta)*(theta+log(v))+x^(1/theta)*(-2+2*theta+log(v)))+(-log(u))^theta*(-log(v))^theta*(x^(3/theta)+x^(2/theta)*(-4+4*theta+log(u)+log(v))+(-1+theta)*(1+theta*(-3+4*theta)+log(v)+log(u)*(1+2*log(v)))+x^(1/theta)*((-1+theta)*(-3+7*theta+log(v))+log(u)*(-1+theta+2*log(v)))))
  return(res)
}

