compute_z<-function(beta){
  z=1*((beta!=0)&!is.na(beta))
}

compute_beta<-function(z, beta_tilde){
  z[is.na(z)]<-0
  z[z==1]<-beta_tilde[!is.na(beta_tilde)]
  z}

sum_var <-function(X){
  temp <-  apply(X, 2, var)
  return (sum(temp)/length(temp))
}

conditional_posterior<- function(R2, q, Y,U,X,sigma2,phi,beta,z,a,b,A,B){
  if(R2==0){
    return(0)
  }
  #z <- compute_z(beta)
  k <- ncol(X)
  v_X <- sum_var(X)
  s_z <- sum(z,na.rm=TRUE)
  res1 <- exp( -1/(2*sigma2)*(k*v_X*q*(1-R2))/R2 * t(beta) %*%  diag(z) %*% beta)
  res2 <- q^(s_z + 0.5*s_z + a -1)
  res3 <-(1-q)^(k-s_z + b -1 )
  res4 <- R2^(A-1-0.5*s_z)
  res5 <- (1-R2)^(0.5*s_z+B-1)
  return (res1*res2*res3*res4*res5)
}


#'@examples
#'a=default_a
#'b=default_b
#'A=aa=default_aa
#'B=bb=default_bb
#'k=default_k
#'s=50
#'r_y=default_r_y[3]
#'X=xx=generate_multiple_x(number_of_datasets = 1)[[1]]
#'Y=y=generate_single_y(s=s,r_y=r_y,xx = xx)
#'sigma2=12;
#'sigma_epsilon=sqrt(sigma2)
#'phi=0
#'z=rep(1:0,each=50)
#'beta<-(`[<-`)(z,1:50,rnorm(50))
#'m=1000
#'aug=draw_conditional_posterior_R2_q(m,Y,U,X,sigma2,phi,beta,z,a,b,A,B)
#'U=u=generate_u();
#'barvx=barvx_f(xx)
#'tbetabeta=t(beta)%*%beta
#'s_z=sum(z,na.rm=TRUE)
#'r2_q_grid=r2_q_grid_generate()
#'dan=sample_r2_q_cond_y_u_x_theta_z(sigma_epsilon,barvx,k,a,b,aa,bb,tbetabeta,s_z,r2_q_grid,m=m)
#'plot(sort(dan[,1]),sort(aug[,1]));abline(0,1,col='red')
#'plot(sort(dan[,2]),sort(aug[,2]));abline(0,1,col='red')



likelihood_conditional_posterior_R2_q<-
  function(R2,q,sigma2,beta,z,a,b,A,B,v_X,k,s_z,correction=FALSE){
    (exp( -1/(2*sigma2)*(k*v_X*q*(1-R2))/R2 * sum(beta[z==1]^2)))*
      (q^(s_z + 0.5*s_z + a -1))*
      ((1-q)^(k-s_z + b -1 ))*
      (R2^(A-1-0.5*s_z))*
      ((1-R2)^(0.5*s_z+B-1))*
      ((1-correction) + (correction*(.001+.009*(0.11 <= q & q <= 0.9))*(.001+.009*(0.11 <= R2 & R2 <= 0.9))))}

dist_conditional_posterior_R2_q<-function(Y,U,X,sigma2,phi,beta,z,a,b,A,B){
   grid <- c(seq(0,0.1,0.001),seq(0.11,0.9,0.01),seq(0.901,1,0.001))
  p = length(grid)
  list_indexes <- vector(mode="numeric", length=p^2)
  list_prob <- vector(mode="numeric", length=p^2)
  list_prob_df<-data.frame(q=numeric(0),R2=numeric(0),prob=numeric(0))
  k <- ncol(X)
  v_X <- sum_var(X)
  s_z <- sum(z,na.rm=TRUE)
  calc_mat<- t(beta) %*%  diag(z) %*% beta
  for (i in 1:p){
    R2<-grid[i]
    res4 <- R2^(A-1-0.5*s_z)
    res5 <- (1-R2)^(0.5*s_z+B-1)
    print(paste0("-",i))
    for(j in 1:p){
      print(j)
      q<-grid[j]
      if(R2==0){
        res<-0
      }else{
        res1 <- exp( -1/(2*sigma2)*(k*v_X*q*(1-R2))/R2 * calc_mat)
        res2 <- q^(s_z + 0.5*s_z + a -1)
        res3 <-(1-q)^(k-s_z + b -1 )
        res<- (res1*res2*res3*res4*res5)
        
        if(0.11 <= q & q <= 0.9){
          res<-res*0.01
        }else{
          res<-res*0.001
        }
        
        
        if(0.11 <= R2 & R2 <= 0.9){
          res<-res*0.01
        }else{
          res<-res*0.001
        }
        
      }
      
      list_prob[i*p+j] <-res
      list_indexes[i*p+j] <- i*p+j
      list_prob_df<-rbind(list_prob_df,data.frame(q=q,R2=R2,prob=res))
    }
  }
  
  
  list_prob_df
}  


grid_a=  expand.grid(
  q=c(seq(0,0.1,0.001),seq(0.11,0.9,0.01),seq(0.901,1,0.001)),
  R2=c(seq(0,0.1,0.001),seq(0.11,0.9,0.01),seq(0.901,1,0.001)))

draw_conditional_posterior_R2_q<-function(sigma2,beta,z,a,b,A,B,v_X,k,s_z,grid_a,m=1){
  grid_a|>
    dplyr::mutate(prob=likelihood_conditional_posterior_R2_q(R2,q,sigma2,beta,z,a,b,A,B,v_X,k,s_z,correction=TRUE),
                  prob=ifelse(is.na(prob),0,prob))|>
    dplyr::slice(sample(dplyr::n(),size=m,prob=prob+(all(prob==0)),replace=TRUE))|>
    (`[`)(c("R2","q"))|>
    as.matrix()|>
    (`[`)(1:m,)
}


def_X_tilde<-function(X,z){
  X[,is.na(z)|(z==1)]
}


#'@examples
#'a=default_a
#'b=default_b
#'A=aa=default_aa
#'B=bb=default_bb
#'k=default_k
#'s=50
#'r_y=default_r_y[3]
#'X=xx=generate_multiple_x(number_of_datasets = 1)[[1]]
#'Y=y=generate_single_y(s=s,r_y=r_y,xx = xx)
#'sigma2=12;
#'sigma_epsilon=sqrt(sigma2)
#'phi=0
#'z=rep(1:0,each=50)
#'beta<-(`[<-`)(z,1:50,rnorm(50))
#'m=1000
#'R2=r2=.8
#'v_X=sum_var(X)
#'U=u=generate_u();
#'barvx=barvx_f(xx)
#'tbetabeta=t(beta)%*%beta
#'s_z=sum(z,na.rm=TRUE)
#'r2_q_grid=r2_q_grid_generate()
#'q=mean(z)
#'gamma2=gamma2_f(r2 = r2,k = k,q = q,barvx = barvx)
#'gamma2_aug=R2/(k*q*v_X*(1-R2))
#'tilde_y=y
#'ttildeytildey=t(tilde_y)%*%tilde_y
#'tt=default_tt
#'aug=sample_conditional_posterior_zi(Y,U,X,sigma2,phi,gamma2_aug,q,z,i=1,m=m)
#'dan=sample_zi_cond_zj_y_u_x_phi_gamma(z,i=1,tilde_y,ttildeytildey,xx,q,tt,k,gamma2,m=m)
#'mean(aug);mean(dan)


sample_conditional_posterior_zi<-function(Y,U,X,sigma2,phi,gamma2,q,z,i,m=1){
  Y_tilde<- Y
  T<-length(Y)
  k<-dim(X)[2]
  
  
  z[i]<-0
  z[is.na(z)]<-0
  s_z<-sum(z,na.rm=TRUE)
  if (s_z==0){
    return(1)
  }else{
  #cas zi=0
  I_sz_0 <- diag(1, nrow = s_z, ncol = s_z)
  X_tilde_0<-def_X_tilde(X,z)
  W_tilde_0<-t(X_tilde_0)%*%X_tilde_0 + I_sz_0/gamma2
  beta_tilde_hat_0<-try(solve(W_tilde_0)%*%t(X_tilde_0)%*%Y_tilde)
  if(is.element("try-error",class(beta_tilde_hat_0))){beta_tilde_hat_0=rep(0,nrow(Y_tilde));
  proba_0<-sum(Y_tilde^2)
  }else{
  proba_0<-sum(Y_tilde^2) - t(beta_tilde_hat_0)%*%W_tilde_0%*%beta_tilde_hat_0 
  }
  
  #cas zi=1
  z[i]<-1
  z[is.na(z)]<-0
  s_z<-sum(z,na.rm=TRUE)
  I_sz_1 <- diag(1, nrow = s_z, ncol = s_z)
  X_tilde_1<-def_X_tilde(X,z)
  W_tilde_1<-t(X_tilde_1)%*%X_tilde_1 + I_sz_1/gamma2
  beta_tilde_hat_1<-try(solve(W_tilde_1)%*%t(X_tilde_1)%*%Y_tilde)
  if(is.element("try-error",class(beta_tilde_hat_1))){beta_tilde_hat_1=rep(0,nrow(Y_tilde))}
  proba_1<- (t(Y_tilde)%*%Y_tilde - t(beta_tilde_hat_1)%*%W_tilde_1%*%beta_tilde_hat_1)
  
  p<-(proba_0/proba_1)^(-T/2)
  p<- det(W_tilde_1)^(-1/2) /( det(W_tilde_1)^(-1/2)  + ((1-q)/q)*gamma2^(1/2)* det(W_tilde_0)^(-1/2) * p )
  
  #print('proba')
  #print(proba_0)
  #print(proba_1)
  #print(p)
  #print('fin')
  return (rbinom(m, size = 1, prob = p))}
}
#'@examples
#'a=default_a
#'b=default_b
#'A=aa=default_aa
#'B=bb=default_bb
#'k=default_k
#'s=50
#'r_y=default_r_y[3]
#'x=X=xx=generate_multiple_x(number_of_datasets = 1)[[1]]
#'Y=y=generate_single_y(s=s,r_y=r_y,xx = xx)
#'sigma2=12;
#'sigma_epsilon=sqrt(sigma2)
#'phi=0
#'z=rep(1:0,each=50)
#'beta<-(`[<-`)(z,1:50,rnorm(50))
#'m=1000
#'R2=r2=.8
#'v_X=sum_var(X)
#'U=u=generate_u();
#'barvx=barvx_f(xx)
#'tbetabeta=t(beta)%*%beta
#'s_z=sum(z,na.rm=TRUE)
#'r2_q_grid=r2_q_grid_generate()
#'q=mean(z)
#'gamma2=gamma2_f(r2 = r2,k = k,q = q,barvx = barvx)
#'gamma2_aug=R2/(k*q*v_X*(1-R2))
#'tilde_y=y
#'ttildeytildey=t(tilde_y)%*%tilde_y
#'tt=default_tt
#'aug=plyr::raply(100,sample_conditional_posterior_z(Y,U,X,phi,R2,q,z))
#'dan=plyr::raply(100,sample_z_cond_y_u_x_phi_gamma(z,tilde_y,ttildeytildey,x,q,tt,k,gamma2))
#'plot(aug|>plyr::aaply(2,mean),dan|>plyr::aaply(2,mean));abline(0,1,col="red")
#'plot(sort(aug|>plyr::aaply(1,mean)),sort(dan|>plyr::aaply(1,mean)));abline(0,1,col="red")

sample_conditional_posterior_z<-function(Y,U,X,phi,R2,q,z0){
  k<-dim(X)[2]
  v_X <- sum_var(X)
  gamma2<-R2/(k*q*v_X*(1-R2))
  for (i in 1:length(z0)) {
    z0[i]<-sample_conditional_posterior_zi(Y,U,X,sigma2,phi,gamma2,q,z0,i)
  }
  return (z0)
}

#'@examples
#'a=default_a
#'b=default_b
#'A=aa=default_aa
#'B=bb=default_bb
#'k=default_k
#'s=50
#'r_y=default_r_y[3]
#'x=X=xx=generate_multiple_x(number_of_datasets = 1)[[1]]
#'Y=y=generate_single_y(s=s,r_y=r_y,xx = xx)
#'sigma2=12;
#'sigma_epsilon=sqrt(sigma2)
#'phi=0
#'z=rep(1:0,each=50)
#'beta<-(`[<-`)(z,1:50,rnorm(50))
#'m=1000
#'R2=r2=.8
#'v_X=sum_var(X)
#'U=u=generate_u();
#'barvx=barvx_f(xx)
#'tbetabeta=t(beta)%*%beta
#'s_z=sum(z,na.rm=TRUE)
#'r2_q_grid=r2_q_grid_generate()
#'q=mean(z)
#'gamma2=gamma2_f(r2 = r2,k = k,q = q,barvx = barvx)
#'gamma2_aug=R2/(k*q*v_X*(1-R2))
#'tilde_y=y
#'Y_tilde=Y
#'ttildeytildey=t(tilde_y)%*%tilde_y
#'tt=default_tt
#'tilde_x=x[,z==1]
#'X_tilde<-def_X_tilde(X,z)
#'k_aug<-dim(X)[2]
#'  I_sz <- diag(1, nrow = sum(z,na.rm=TRUE), ncol = sum(z,na.rm=TRUE))
#'v_X <- sum_var(X)
#'gamma2<-R2/(k*q*v_X*(1-R2))
#'W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2
#'tilde_w=(t(tilde_x)%*%tilde_x+diag(s_z)/gamma2)
#'inv_tilde_w=if(s_z==0){0}else{solve(tilde_w)}
#'hat_tilde_beta=inv_tilde_w%*%t(tilde_x)%*%(tilde_y)
#'range(tilde_w-W_tilde)
#'range(tilde_x-X_tilde)
#'range(tilde_y-Y_tilde)
#'shape<-length(Y)/2 
#'gamma2<-R2/(k*q*v_X*(1-R2))
#'W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2_aug
#'beta_tilde_hat <- solve(W_tilde)%*%t(X_tilde)%*%Y_tilde
#'range(hat_tilde_beta-beta_tilde_hat)
#'scale<-t(beta_tilde_hat)%*%W_tilde%*%beta_tilde_hat
#'scale<- t(Y_tilde)%*%Y_tilde - scale
#'scale_aug<-scale/2
#'rate_dan=(ttildeytildey-t(hat_tilde_beta)%*%tilde_w%*%hat_tilde_beta)/2
#'aug=sample_conditional_posterior_sigma2(Y,U,X,phi,R2, q, z,m=m)|>sort()
#'dan=sample_sigmaepsilon_cond_y_u_x_phi_r2_q_z(z,tt,ttildeytildey,tilde_w,hat_tilde_beta,m=m)|>(`^`)(2)|>sort()
#'plot(aug,dan);abline(0,1,col='red')

sample_conditional_posterior_sigma2<-function(Y,U,X,phi,R2, q, z,m=1){
  shape<-length(Y)/2 
  I_sz <- diag(1, nrow = sum(z,na.rm=TRUE), ncol = sum(z,na.rm=TRUE))
  Y_tilde<- Y
  z[is.na(z)]<-0
  X_tilde<-def_X_tilde(X,z)
  k<-dim(X)[2]
  v_X <- sum_var(X)
  gamma2<-R2/(k*q*v_X*(1-R2))
  W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2
  beta_tilde_hat <- solve(W_tilde)%*%t(X_tilde)%*%Y_tilde
  scale<-t(beta_tilde_hat)%*%W_tilde%*%beta_tilde_hat
  scale<- t(Y_tilde)%*%Y_tilde - scale
  scale<-scale/2
  echantillon_inverse_gamma <- invgamma::rinvgamma(m, shape = shape, rate = scale)
  return (echantillon_inverse_gamma)
}



## **Question 2 - (V)**
#'@examples
#'a=default_a
#'b=default_b
#'A=aa=default_aa
#'B=bb=default_bb
#'k=default_k
#'s=50
#'r_y=default_r_y[3]
#'x=X=xx=generate_multiple_x(number_of_datasets = 1)[[1]]
#'Y=y=generate_single_y(s=s,r_y=r_y,xx = xx)
#'sigma2=12;
#'sigma_epsilon=sqrt(sigma2)
#'phi=0
#'z=rep(1:0,each=50)
#'beta<-(`[<-`)(z,1:50,rnorm(50))
#'m=1000
#'R2=r2=.8
#'v_X=sum_var(X)
#'U=u=generate_u();
#'barvx=barvx_f(xx)
#'tbetabeta=t(beta)%*%beta
#'s_z=sum(z,na.rm=TRUE)
#'r2_q_grid=r2_q_grid_generate()
#'q=mean(z)
#'gamma2=gamma2_f(r2 = r2,k = k,q = q,barvx = barvx)
#'gamma2_aug=R2/(k*q*v_X*(1-R2))
#'tilde_y=y
#'Y_tilde=Y
#'ttildeytildey=t(tilde_y)%*%tilde_y
#'tt=default_tt
#'tilde_x=x[,z==1]
#'X_tilde<-def_X_tilde(X,z)
#'k_aug<-dim(X)[2]
#'  I_sz <- diag(1, nrow = sum(z,na.rm=TRUE), ncol = sum(z,na.rm=TRUE))
#'v_X <- sum_var(X)
#'gamma2<-R2/(k*q*v_X*(1-R2))
#'W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2
#'tilde_w=(t(tilde_x)%*%tilde_x+diag(s_z)/gamma2)
#'inv_tilde_w=if(s_z==0){0}else{solve(tilde_w)}
#'hat_tilde_beta=inv_tilde_w%*%t(tilde_x)%*%(tilde_y)
#'range(tilde_w-W_tilde)
#'range(tilde_x-X_tilde)
#'range(tilde_y-Y_tilde)
#'shape<-length(Y)/2 
#'gamma2<-R2/(k*q*v_X*(1-R2))
#'W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2_aug
#'beta_tilde_hat <- solve(W_tilde)%*%t(X_tilde)%*%Y_tilde
#'range(hat_tilde_beta-beta_tilde_hat)
#'scale<-t(beta_tilde_hat)%*%W_tilde%*%beta_tilde_hat
#'scale<- t(Y_tilde)%*%Y_tilde - scale
#'scale_aug<-scale/2
#'rate_dan=(ttildeytildey-t(hat_tilde_beta)%*%tilde_w%*%hat_tilde_beta)/2
#'mean_vector<-solve(W_tilde)%*%t(X_tilde)%*%Y_tilde
#'covariance_matrix<-sigma2*solve(W_tilde)
#'mu = inv_tilde_w%*%t(tilde_x)%*%(tilde_y)
#'Sigma = (sigma_epsilon^2)*inv_tilde_w
#'range(mu-mean_vector)
#'range(Sigma-covariance_matrix)
#'aug=sample_conditional_posterior_beta_tilde(Y,U,X,phi,R2,q,sigma2,z,m=m)
#'dan=sample_tildebeta_cond_y_u_x_phi_r2_q_z_sigma2(sigma_epsilon,tilde_y,tilde_x,tilde_w,m=m)
#'plot(aug[1,],dan[1,]);abline(0,1,col='red')
#'plot(aug[,1]|>sort(),dan[,1]|>sort());abline(0,1,col='red')


sample_conditional_posterior_beta_tilde<-function(Y,U,X,phi,R2,q,sigma2,z,m=1){
  z[is.na(z)]<-0
  I_sz <- diag(1, nrow = sum(z,na.rm=TRUE), ncol = sum(z,na.rm=TRUE))
  Y_tilde<- Y
  X_tilde<-def_X_tilde(X,z)
  k<-dim(X)[2]
  v_X <- sum_var(X)
  gamma2<-R2/(k*q*v_X*(1-R2))  
  W_tilde<-t(X_tilde)%*%X_tilde + I_sz/gamma2
  
  mean_vector<-solve(W_tilde)%*%t(X_tilde)%*%Y_tilde
  covariance_matrix<-sigma2*solve(W_tilde)
  multivariate_data <- MASS::mvrnorm(n = m, mu = mean_vector, Sigma = covariance_matrix)
  return (multivariate_data)
}

sample_data0<-function(TT,k,rho,s,Ry){
  require(glmnet)
  corr_matrix <- toeplitz(rho^(0:(k-1)))
  mu <- rep(0, times = k)
  X<-mvrnorm(n = TT, mu = mu, Sigma = corr_matrix)
  X <- scale(X)
  beta<-mvrnorm(n = k, mu = 0, Sigma = 1)
  indices_zeros <- sample(length(beta), size = k-s)
  beta[indices_zeros] <- 0
  sigma2<-sum(sapply(X%*%beta, function(x) x^2))*(1/Ry-1)/TT
  epsilon<-mvrnorm(n = TT, mu = 0, Sigma = sigma2)
  Y<-X%*%beta + epsilon
  return(list(X=X,Y=Y))}

init_betasigma2<-function(X,Y){
    lasso_model <- glmnet::cv.glmnet(X, Y, alpha = 1, intercept = FALSE)
  beta<-(coef(lasso_model, s = "lambda.min")[-1, 1])
  
  sigma2<-var(Y-X%*%beta)
  #sigma2<-((y-x%*%beta)^2)|>sum()|>(`/`)(nrow(x)-sum(beta!=0))
  
  return(list(beta=beta,sigma2=sigma2))
}
Gibbs<-function(N,a,A,b,B,k,U,phi,X,Y){
  v_X <- sum_var(X)
  
  init<-init_betasigma2(X,Y)
  beta<-init$beta
  sigma2<-init$sigma2
  phi<-0
  z <- compute_z(beta)
  s_z<-sum(z,na.rm=TRUE)
  list_R2 <- vector(mode="numeric", length=N)
  list_q <- vector(mode="numeric", length=N)
  list_s_z <- vector(mode="numeric", length=N)
  list_sigma_epsilon <- vector(mode="numeric", length=N)
  for(i in 1:N){
    if (i %% 50 == 0) {
      print(i)
    }
    mat_R2q<-draw_conditional_posterior_R2_q(sigma2,beta,z,a,b,A,B,v_X,k,s_z,grid_a,m=1)
    R2<-mat_R2q["R2"]
    q<-mat_R2q["q"]
    z<-sample_conditional_posterior_z(Y,U,X,phi,R2,q,z)
    s_z<-sum(z,na.rm=TRUE)
    sigma2<-sample_conditional_posterior_sigma2(Y,U,X,phi,R2, q, z)
    beta_tilde<-sample_conditional_posterior_beta_tilde(Y,U,X,phi,R2,q,sigma2,z)
    beta<-compute_beta(z,beta_tilde)
    list_R2[i]<-R2
    list_q[i]<-q
    list_s_z[i]<-s_z
    list_sigma_epsilon[i]=sqrt(sigma2)
  }
  cbind(q=list_q,
        R2=list_R2,
        s_x=list_s_z,
        sigma_epsilon=list_sigma_epsilon)
}
if(F){
  
  
  
  N<-6000
  
  a<-1
  A<-1
  b<-1
  B<-1
  
  TT<-200
  k<-100
  rho<-0.75
  U<-0
  phi<-0
  
  Ry<-0.25
  s<-10
  res<-sample_data0(T,k,rho,s,Ry)
  gibbs_sample<-Gibbs(N,a,A,b,B,TT,k,rho,U,phi,Ry,s,X,Y)
    

hist(gibbs_sample$R2, main = "Histogramme de R2", xlab = "Valeurs", ylab = "Fréquence", col = "skyblue", border = "black")
  hist(gibbs_sample$q, main = "Histogramme de q", xlab = "Valeurs", ylab = "Fréquence", col = "skyblue", border = "black")
}
