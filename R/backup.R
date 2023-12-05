
#'@examples
#'a = c(0,1)
#'print(a)
#'print(compute_z(a))
compute_z<-function(beta){
  n = length(beta)
  z <- vector(mode="numeric", length=n)
  for (i in 1:n){
    if(beta[i]==0){
      z[i]<-0
    }else{
      z[i]<-1
    }
  }
  return (z)
}
sum_var <-function(X){
  temp <-  apply(X, 2, var)/nrow(X)
  return (sum(temp)/length(temp))
}

#'@examples
#'
#'

conditional_posterior<- function(R2, q, Y,U,X,sigma2,phi,beta,z,a,b,A,B){
  if(R2==0){
    return(0)
  }
  z <- compute_z(beta)
  k <- ncol(X)
  v_x <- sum_var(X)
  s_z <- sum(z)
  res1 <- exp( -1/(2*sigma2)*(k*v_x*q*(1-R2))/R2 * t(beta) %*%  diag(z) %*% beta)
  if (is.na(res1)) {
    print("res1**")
    print(v_x)
    print(sigma2)
    print(k)
    print(q)
    print(R2)
    print(t(beta) %*%  diag(z) %*% beta)
  }
  res2 <- q^(s_z + 0.5*s_z + a -1)
  if (is.na(res2)) {
    print("res2")
  }
  res3 <-(1-q)^(k-s_z + b +1 )
  if (is.na(res3)) {
    print("res3")
  }
  res4 <- R2^(A-1-0.5*s_z)
  if (is.na(res4)) {
    print("res4")
  }
  res5 <- (1-R2)^(0.5*s_z+B-1)
  if (is.na(res5)) {
    print("res5")
  }
  
  return (res1*res2*res3*res4*res5)
}

#'@examples
#'beta=1
#'X=1
#'draw_conditional_posterior_R2_q(Y,U,X,sigma2,phi,beta,z,1,1,1,1)







draw_conditional_posterior_R2_q<-function(Y,U,X,sigma2,phi,beta,z,a,b,A,B){
  m<-1
  grid <- c(seq(0,0.1,0.001),seq(0.11,0.9,0.01),seq(0.901,1,0.001))
  p = length(grid)
  list_indexes <- vector(mode="numeric", length=p^2)
  list_prob <- vector(mode="numeric", length=p^2)
  for (i in 1:p){
    for(j in 1:p){
      list_prob[i*p+j] <-conditional_posterior(grid[i],grid[j], Y,U,X,sigma2,phi,beta,z,a,b,A,B)
      list_indexes[i*p+j] <- i*p+j
    }
  }
  empty_matrix <- matrix(0, 4,  6)
  list_samples <- sample(list_indexes,m,replace =TRUE ,prob=list_prob)
  print(2)
  matrix_samples <-matrix(0,m,2)
  for (i in 1:m){
    matrix_samples[i,1]<- grid[list_samples[m]%/%p]
    matrix_samples[i,2]<-grid[list_samples[m]%%p]
  }
  return (matrix_samples)
}

#'@examples
#'draw_conditional_posterior_R2_q(Y,U,X,sigma2,phi,beta,z,1,1,1,1)
#'#test des fonctions ci dessus 
#' Y<-c(1, 0.5,1, 0.5)
#' U<-matrix(c(1, 0.5, 0.5, 1, 0.5, 1, 0.5, 1), nrow = 4, byrow = TRUE)
#' X<-matrix(c(1, 0.2, 0.32, 1, 0.2, 0.45, 0, 1.3, 0.6, 0.45, 0, 2, 0.6, 45, 0, 0.2, 0.6, 0.45, 0, 2, 6, 0.45, 0, 0.2), nrow = 4, byrow = TRUE)
#' phi<-c(1.2,0.2)
#' beta<-c(0,0,1,0,0,1)
#' z<-c(0,0,1,0,0,1)
#' sigma2<-1
#' R2<-0.9
#' q<-0.15

#' draw_conditional_posterior_R2_q(Y,U,X,phi,beta,z,1,1,1,1)



## **Question 2 - (II)**


draw_conditional_posterior_phi<-function(Y,U,X,sigma2,beta){
  mean_vector<-solve(t(U) %*% U)%*%t(U)%*%(Y-X%*%beta)
  covariance_matrix<- solve(t(U) %*% U)*sigma2
  multivariate_data <- mvrnorm(n = 1, mu = mean_vector, Sigma = covariance_matrix)
  return (multivariate_data)
}
#'@examples
#'Y<-c(1, 0.5,1, 0.5)
#'X<-matrix(c(1, 0.5, 0.5, 1, 0.5, 1, 0.5, 1), nrow = 4, byrow = TRUE)
#'U<-matrix(c(1, 0.2, 0.32, 1, 126, 0.45, 0, 1.3, 6, 45, 0, 2), nrow = 4, byrow = TRUE)
#'sigma2<-1
#'beta<-c(1, 0.5)
#'draw_conditional_posterior_phi(Y,U,X,sigma2,beta)


## **Question 2 - (III)**
def_X_tilde<-function(X,z){
  X_tilde <- diag(0, nrow = nrow(X), ncol = sum(z))
  cpt<-0
  for (i in 1:ncol(X)) {
    if (z[i]==1) {
      cpt<-cpt+1
      X_tilde[,cpt]<-X[, i]
    }
  }
  return (X_tilde)
}


sample_conditional_posterior_zi<-function(Y,U,X,sigma2,phi,gamma2,z,i){
  Y_tilde<- Y - U%*%phi
  T<-length(Y)
  k<-dim(X)[2]
  
  #cas zi=0
  z[i]<-0
  s_z<-sum(z)
  if (s_z==0){
    return(1)
  }else{
    I_sz <- diag(1, nrow = sum(z), ncol = sum(z))
    X_tilde<-def_X_tilde(X,z)
    W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2
    beta_tilde_hat<-solve(W_tilde)%*%t(X_tilde)%*%Y_tilde
    proba_0<-((t(Y_tilde)%*%Y_tilde - t(beta_tilde_hat)%*%W_tilde%*%beta_tilde_hat)/2 )
    proba_0<-proba_0^(-T/2)
    proba_0<-proba_0*det(W_tilde)^0.5*gamma(T/2)*q^s_z*(1-q)^(k-s_z)*(1/gamma2)^(s_z/2)
  }
  
  
  #cas zi=1
  z[i]<-1
  s_z<-sum(z)
  I_sz <- diag(1, nrow = sum(z), ncol = sum(z))
  X_tilde<-def_X_tilde(X,z)
  W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2
  beta_tilde_hat<-solve(W_tilde)%*%t(X_tilde)%*%Y_tilde
  proba_1<- (t(Y_tilde)%*%Y_tilde - t(beta_tilde_hat)%*%W_tilde%*%beta_tilde_hat)/2 
  proba_1<-proba_1^(-T/2)
  proba_1<-proba_1*det(W_tilde)^0.5*gamma(T/2)*q^s_z*(1-q)^(k-s_z)*(1/gamma2)^(s_z/2)
  
  p<-proba_1/(proba_1+proba_0)
  #print('proba')
  #print(proba_0)
  #print(proba_1)
  #print(p)
  return (rbinom(1, size = 1, prob = p))
}

sample_conditional_posterior_z<-function(Y,U,X,phi,R2,q,z0){
  k<-dim(X)[2]
  v_X <- sum_var(X)
  gamma2<-R2/(k*q*v_X*(1-R2))
  for (j in 1:100) {
    for (i in 1:length(z0)) {
      z0[i]<-sample_conditional_posterior_zi(Y,U,X,sigma2,phi,gamma2,z0,i)
    }
  }
  return (z0)
}

#'@examples
#'#test des fonctions ci dessus 
#'Y<-c(1, 0.5,1, 0.5)
#'U<-matrix(c(1, 0.5, 0.5, 1, 0.5, 1, 0.5, 1), nrow = 4, byrow = TRUE)
#'X<-matrix(c(1, 0.2, 0.32, 1, 0.2, 0.45, 0, 1.3, 0.6, 0.45, 0, 2, 0.6, 45, 0, 0.2, 0.6, 0.45, 0, 2, 6, 0.45, 0, 0.2), nrow = 4, byrow = TRUE)
#'phi<-c(1.2,0.2)
#'z0<-c(0,0,1,0,0,1)
#'R2<-0.9
#'q<-0.15
#'
#'#sample_conditional_posterior_z(Y,U,X,phi,R2,q,z0)



## **Question 2 - (IV)**

sample_conditional_posterior_sigma2<-function(Y,U,X,phi,R2, q, z){
  shape<-length(Y)/2 
  I_sz <- diag(1, nrow = sum(z), ncol = sum(z))
  Y_tilde<- Y - U%*%phi
  X_tilde<-def_X_tilde(X,z)
  k<-dim(X)[2]
  v_X <- sum_var(X)
  gamma2<-R2/(k*q*v_X*(1-R2))  
  W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2
  beta_tilde_hat <- solve(W_tilde)%*%t(X_tilde)%*%Y_tilde
  scale<-t(beta_tilde_hat)%*%W_tilde%*%beta_tilde_hat
  scale<- t(Y)%*%Y - scale
  scale<-scale/2
  echantillon_inverse_gamma <- rinvgamma(1, shape = shape, rate = 1/scale)
  return (echantillon_inverse_gamma)
}

#test des fonctions ci dessus 
#'@examples
#'Y<-c(1, 0.5,1, 0.5)
#'U<-matrix(c(1, 0.5, 0.5, 1, 0.5, 1, 0.5, 1), nrow = 4, byrow = TRUE)
#'X<-matrix(c(1, 0.2, 0.32, 1, 0.2, 0.45, 0, 1.3, 0.6, 0.45, 0, 2, 0.6, 45, 0, 0.2, 0.6, 0.45, 0, 2, 6, 0.45, 0, 0.2), nrow = 4, byrow = TRUE)
#'phi<-c(1.2,0.2)
#'z<-c(0,0,1,0,0,1)
#'sigma2<-1
#'R2<-0.9
#'q<-0.15

#'draw_conditional_posterior_R2_q(Y,U,X,sigma2,phi,beta,z,1,1,1,1)



## **Question 2 - (V)**

sample_conditional_posterior_beta_tilde<-function(Y,U,X,phi,R2,q,sigma2,z){
  I_sz <- diag(1, nrow = sum(z), ncol = sum(z))
  Y_tilde<- Y - U%*%phi
  X_tilde<-def_X_tilde(X,z)
  k<-dim(X)[2]
  v_X <- sum_var(X)
  gamma2<-R2/(k*q*v_X*(1-R2))  
  W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2
  
  mean_vector<-solve(W_tilde)%*%t(X_tilde)%*%Y_tilde
  covariance_matrix<-sigma2*solve(W_tilde)
  multivariate_data <- mvrnorm(n = 1, mu = mean_vector, Sigma = covariance_matrix)
  return (multivariate_data)
}
#test des fonctions ci dessus 
#'@examples
#'Y<-c(1, 0.5,1, 0.5)
#'U<-matrix(c(1, 0.5, 0.5, 1, 0.5, 1, 0.5, 1), nrow = 4, byrow = TRUE)
#'X<-matrix(c(1, 0.2, 0.32, 1, 0.2, 0.45, 0, 1.3, 0.6, 0.45, 0, 2, 0.6, 45, 0, 0.2, 0.6, 0.45, 0, 2, 6, 0.45, 0, 0.2), nrow = 4, byrow = TRUE)
#'phi<-c(1.2,0.2)
#'z<-c(0,0,1,0,0,1)
#'R2<-0.9
#'q<-0.15
#'sigma2<-1


#'sample_conditional_posterior_beta_tilde(Y,U,X,phi,R2,q,sigma2,z)

sample_data0<-function(T,k,rho,s,Ry){#tu recalcules corr_matrix chaque fois
  corr_matrix <- toeplitz(rho^(0:(k-1)))
  mu <- rep(0, times = k)
  X<-mvrnorm(n = T, mu = mu, Sigma = corr_matrix)
  beta<-mvrnorm(n = k, mu = 0, Sigma = 1)
  indices_zeros <- sample(length(beta), size = s)
  beta[indices_zeros] <- 0
  sigma2<-sum(sapply(X%*%beta, function(x) x^2))*(1/Ry-1)/T
  epsilon<-mvrnorm(n = T, mu = 0, Sigma = sigma2)
  Y<-X%*%beta + epsilon
  return(list(X,Y,beta))
}
res<-sample_data0(5,3,0.75,1,0.25)
X<-res[1]
Y<-res[2]
beta<-res[3]

res<-sample_data0(200,100,0.75,5,0.25)
X<-res[1]
Y<-res[2]
beta<-res[3]

for(i in 1:10){
  
}



