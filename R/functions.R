default_k=100
default_tt=200

default_a=1
default_b=1
default_aa=1
default_bb=1

default_rho=.75
default_number_of_datasets=100
default_s=c(5,10,100)
default_r_y=c(.02,.25,.5)

default_nrep=110000
default_burning=10000

#'@description
#' Randomly generates multiple dataset version of the matrix X.
#'@examples
#'k=default_k
#'tt=default_tt
#'rho=default_rho
#'number_of_datasets=default_number_of_datasets
#'x<-generate_multiple_x(k=k,tt=tt,rho=rho,number_of_datasets=number_of_datasets)
generate_multiple_x<-
  function(k=default_k,
           tt=default_tt,
           rho=default_rho,
           number_of_datasets=3){
    rho=.75
    big_sigma=stats::toeplitz(.75^(0:(k-1)))
    #toeplitz=rho^(outer(1:.T,1:.T,`-`)|>abs())
    plyr::rlply(number_of_datasets,
                {mvrnorm(n = tt,mu = rep(0,k),Sigma=big_sigma)|>
                    plyr::aaply(2,function(x){(x-mean(x))/sd(x)})})}

#'@examples
#'k=default_k
#'s=5
#'draw_beta(s)
draw_beta<-function(s){c(rnorm(s))}


#'@examples
#'k=default_k
#'s=default_s
#'r_y=default_r_y
#'rho=default_rho
#'number_of_datasets=default_number_of_datasets
#'x<-generate_y_sample_q(k=k,tt=tt,rho=rho,number_of_datasets=number_of_datasets)
#'generate_single_y()
generate_single_y<-function(
  s=default_s[1],
  r_y=default_r_y[1],
  xx=generate_multiple_x(k=default_k,
                        tt=default_tt,
                        rho=default_rho,
                        number_of_datasets=1)[[1]]){
                  beta=draw_beta(s)
                  xbeta=xx[,1:s]%*%beta
                  #xbeta=xx[,sample(ncol(xx),size = d$s,replace=FALSE)]%*%beta#we could take the n first one would do the same
                  sigma_epsilon<-sqrt((1/r_y-1)/T)*abs(xbeta)
                  y=xbeta +rnorm(nrow(xx),sd=sigma_epsilon)
                  y
                  }
  
  generate_y_sample_q<-function(
  s=default_s,
  r_y=default_r_y,
  a=default_a,
  b=default_b,
  aa=default_aa,
  bb=default_bb,
  nrep=default_nrep,
  burning=default_burning,
  x=generate_multiple_x(k=default_k,
                        tt=default_tt,
                        rho=default_rho,
                        number_of_datasets=default_number_of_datasets)){
  expand.grid(s=s,r_y=r_y,i=1:length(xx))|>
      plyr::daply(~s+r_y+i,
                function(d){
                  xx<-x[[i]]
                  beta=draw_beta(d$s)
                  xbeta=xx[,1:d$s]%*%beta
                  #xbeta=xx[,sample(ncol(xx),size = d$s,replace=FALSE)]%*%beta#we could take the n first one would do the same
                  sigma_epsilon<-sqrt((1/d$r_y-1)/T)*abs(xbeta)
                  y=xbeta +rnorm(nrow(xx),sd=sigma_epsilon)
                  q<-sample_q(xx,y,a=a,
                              b=b,
                              aa=aa,
                              bb=bb,
                              nrep=nrep,
                              burning=burning)})
  
  
}
  

## **Question 2 - (I)**

#'@examples
#'empty_matrix <- matrix(0, 4,  6)
#'print(empty_matrix)
#'empty_matrix[1,1:2] <- c(0,1)
#'print(empty_matrix)

#'@examples
#'r2_q_grid_generate()|>nrow()|>sqrt()
#'
r2_q_grid_generate<-function(){
  base_grid=c(seq(0,.1,by=.001),
             seq(.11,.9,by=.01),
             seq(.901,1,by=.001))
  dbase_grid=c(.0005,rep(.001,99),.0055,
              rep(.01,79),.0055,
              rep(.001,99),.0005)
  cbind(
  expand.grid(q=base_grid,
              r2=base_grid),
  expand.grid(dq=dbase_grid,
              dr2=dbase_grid))|>
    dplyr::mutate(dqdr2=dq*dr2)|>
    dplyr::select(q,r2,dqdr2)
}




#'@examples
#'r2=.5;q=.5;u=0;x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,195));sigma=1;beta=rnorm(z)*z;y=generate_single_y(xx=x)
#'loglikelihood_r2_q_cond_y_u_x_theta_z_a_b_aa_bb(r2,q,y,u,x,sigma,beta,z)

loglikelihood_r2_q_cond_y_u_x_theta_z_a_b_aa_bb<-
  function(r2,q,y,u,x,sigma,beta,z,barvx=1,
           k=ncol(x),  
           a=default_a,
           b=default_b,
           aa=default_aa,
           bb=default_bb,
           s_z=sum(z)){
    -(2*sigma^2)^(-1)*k*barvx*q*(1-r2)*r2^(-1)*c(t(beta)%*%diag(z)%*%beta)+
      log(q^(s_z*(3/2)+a-1))+
      log((1-q)^(k-s_z+b-1))+
      log(r2^(aa-1-(s_z/2)))+
      log((1-r2)^((s_z/2)+bb-1))
  }




#'@examples
#'u=0;x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,195));sigma=1;beta=z*rnorm(z);y=generate_single_y(xx=x)
#'sample_r2_q_cond_y_u_x_theta_z(r2,q,y,u,x,sigma,beta,z)

sample_r2_q_cond_y_u_x_theta_z<-
  function(y,u,x,sigma,beta,z,
           k=default_k,
           r2_q_grid=r2_q_grid_generate()){
  r2_q_grid|>
    dplyr::mutate(logpi=loglikelihood_r2_q_cond_y_u_x_theta_z_a_b_aa_bb(
      r2,q,y,u,x,sigma,beta,z,barvx=1,
                             k=ncol(x),  
                             a=default_a,
                             b=default_b,
                             aa=default_aa,
                             bb=default_bb,
                             s_z=sum(z)),
      pi=dqdr2*#we multiply the probabilities by dq dr2 as the grid is not regular
        exp(logpi-max(logpi,na.rm=TRUE))|>
        (function(x){ifelse(is.na(x),0,x)})())|>
      dplyr::slice(sample(dplyr::n(),size=1,prob=pi))|>
      (`[`)(c("r2","q"))|>
      unlist()
}

## III.

#'@examples
#'r2=.5;q=.5;u=0;x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,195));sigma=1;beta=rnorm(z)*z;y=generate_single_y(xx=x)
#'loglikelihood_r2_q_cond_y_u_x_theta_z_a_b_aa_bb(r2,q,y,u,x,sigma,beta,z)

loglikelihood_z_cond_y_u_x_phi_gamma<-
  function(y,u,x,phi,gamma,tt=nrow(tt),s_z=sum(z),k=ncol(x)){
    tilde_w=
    (q^s_z)*((1-q)^(k-s_z))*
      gamma^(-s_z)*
      det(tilde_w)^(-1/2)*
      gamma(tt/2)
    
  }





### IV.
