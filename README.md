Assignment 2
================
Daniel Bonnery, Max Kramkimel, Augustin Poissonnier
Dec 5, 2023

# Bayesian Statistics, Assignment 2

## Results

### Plot 1

### Plot 2

## How to run

### Install the following packages:

dplyr, ggplot2, glmnet, MASS, plyr, stats

### Copy and execute the project:

-   Clone the project to a local repository via:

``` r
system("git clone https://github.com/Augustin202/Stat_bayes.git")
```

-   open the .rproj file in rstudio
-   run:

``` r
source(main.R)
```

## R code

The R code consists a file R/functions.R that contains all necessary
functions to run the main code, a main.R file that will run the main
program. The main file generates the data, draw samples of q, creates
plots and compiles the Rmd file to create this markdown file, and an
Rproj file, all available on github there:.

### main.R

### R/functions.R

``` r
#fixed parameters
default_l=0
default_k=100
default_tt=200
default_rho=.75

#hyperpriors
default_a=1
default_b=1
default_aa=1
default_bb=1

#Simulation parameters
default_number_of_datasets=100
default_s=c(5,10,100)
default_r_y=c(.02,.25,.5)

#Sampling parameters
default_nrep=5000
default_burning=1000

#'@description
#' Compute gamma2 as a function of k, q, barvx
#'@examples
#'r2=.5;k=default_k,q=.1,barvx=1
#'gamma2_f(r2,k,q,barvx)
gamma2_f<-function(r2,k,q,barvx){r2/(k*q*barvx*(1-r2))}

#'@description
#' Randomly generates multiple dataset version of the matrix X.
#'@examples
#'k=4
#'tt=10
#'rho=default_rho
#'number_of_datasets=3
#'generate_multiple_x(k=k,tt=tt,rho=rho,number_of_datasets=number_of_datasets)
generate_multiple_x<-
  function(k=default_k,
           tt=default_tt,
           rho=default_rho,
           number_of_datasets=3){
    big_sigma=stats::toeplitz(.75^(0:(k-1)))
    #toeplitz=rho^(outer(1:.T,1:.T,`-`)|>abs())
    plyr::rlply(number_of_datasets,
                {MASS::mvrnorm(n = tt,mu = rep(0,k),Sigma=big_sigma)|>
                    plyr::aaply(2,function(x){(x-mean(x))/sd(x)})|>
                    aperm(2:1)})}

#'@examples
#'xx=generate_multiple_x(k=default_k,tt=default_tt,rho=default_rho,number_of_datasets=1)[[1]]
#'barvx_f(xx)
barvx_f<-function(x){x|>plyr::aaply(2,var)|>mean()}

#'@examples
#'s=5
#'draw_beta(s)
draw_tildebeta<-rnorm


#'@examples
#'k=default_k
#'s=default_s
#'r_y=default_r_y
#'rho=default_rho
#'number_of_datasets=default_number_of_datasets
#'xx=generate_multiple_x(k=default_k,tt=default_tt,rho=default_rho,number_of_datasets=1)[[1]]
#'generate_single_y(s[1],r_y[1],xx)
generate_single_y<-function(
  s=default_s[1],
  r_y=default_r_y[1],
  xx=generate_multiple_x(k=default_k,
                        tt=default_tt,
                        rho=default_rho,
                        number_of_datasets=1)[[1]]){
                  beta=draw_tildebeta(s)
                  xbeta=xx[,1:s]%*%beta
                  #xbeta=xx[,sample(ncol(xx),size = d$s,replace=FALSE)]%*%beta#we could take the n first one would do the same
                  sigma_epsilon<-sqrt((1/r_y-1)*mean(xbeta^2))
                  y=xbeta +rnorm(nrow(xx),sd=sigma_epsilon)
                  y}
#'@examples
#'x=generate_multiple_x(k=default_k,tt=default_tt,rho=default_rho,number_of_datasets=3)
#'nrep=30
#'burning=10
#'s=default_s
#'r_y=default_r_y
#'rho=default_rho
#'a=default_a
#'b=default_b
#'aa=default_aa
#'bb=default_bb
#'r2_q_grid=r2_q_grid_generate()
#'generate_y_sample_q(x,s,r_y,a,b,aa,bb,r2_q_grid,nrep,burning)

generate_y_sample_q<-function(
  x,
  u,
  s,
  r_y,
  a,
  b,
  aa,
  bb,
  r2_q_grid,
  nrep,
  burning){
  expand.grid(s=s,r_y=r_y,i=1:length(x))|>
      plyr::daply(~s+r_y+i,
                function(d){
                  xx<-x[[d$i]]
                  tt=nrow(xx)
                  k=ncol(xx)
                  beta=draw_tildebeta(d$s)
                  xbeta=xx[,1:d$s]%*%beta
                  #xbeta=xx[,sample(ncol(xx),size = d$s,replace=FALSE)]%*%beta#we could take the n first one would do the same
                  sigma_epsilon<-sqrt((1/r_y-1)*mean(xbeta^2))
                  y=xbeta +rnorm(nrow(xx),sd=sigma_epsilon)
                  barvx=barvx_f(xx)
                  q<-Gibbs_q(x = xx,y = y,u = u,barvx = barvx,
                             tt = tt,k = k,phi = phi,r2_q_grid = r2_q_grid,a=a,
                              b=b,
                              aa=aa,
                              bb=bb,
                              nrep=nrep,
                              burning=burning)},.progress="text")|>
    (function(x){names(dimnames(x))<-c("s","r_y","dataset","j");x})()

}
  

## **Question 2 - (I)**

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


generate_u<-function(tt=default_tt){u=matrix(0,tt,1)}

#'@examples
#'r2=.5;q=.5;u=generate_u();x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,95));sigma_epsilon=1;beta=rnorm(z)*z;y=generate_single_y(xx=x)
#'tbetabeta=t(beta)%*%beta
#'loglikelihood_r2_q_cond_y_u_x_theta_z_a_b_aa_bb(r2,q,sigma_epsilon,barvx,k,a,b,aa,bb,tbetabeta,s_z)

loglikelihood_r2_q_cond_y_u_x_theta_z_a_b_aa_bb<-
  function(r2,q,sigma_epsilon,barvx,k,a,b,aa,bb,tbetabeta,s_z){
    -(2*sigma_epsilon^2)^(-1)*k*barvx*q*(1-r2)*(r2^(-1))*c(tbetabeta)+
      (s_z*(3/2)+a-1)*log(q)+
      (k-s_z+b-1)*log(1-q)+
      (aa-1-(s_z/2))*log(r2)+
      ((s_z/2)+bb-1)*log(1-r2)
  }


sample_phi<-function(){0}

#'@examples
#'u=generate_u();x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,95));sigma_epsilon=1;beta=z*rnorm(z);y=generate_single_y(xx=x)
#'sample_r2_q_cond_y_u_x_theta_z(r2,q,sigma_epsilon,barvx,k,a,b,aa,bb,tbetabeta,s_z,r2_q_grid)
sample_r2_q_cond_y_u_x_theta_z<-
  function(sigma_epsilon,barvx,k,a,b,aa,bb,tbetabeta,s_z,r2_q_grid){
  r2_q_grid|>
    dplyr::mutate(logpi=loglikelihood_r2_q_cond_y_u_x_theta_z_a_b_aa_bb(
      r2=r2,q=q,sigma_epsilon=sigma_epsilon,barvx=barvx,k=k,a=a,b=b,aa=aa,bb=bb,tbetabeta=tbetabeta,s_z=s_z),
      pi=dqdr2*#we multiply the probabilities by dq dr2 as the grid is not regular
        exp(logpi-max(logpi,na.rm=TRUE))|>
        (function(x){ifelse(is.na(x),0,x)})())|>
      dplyr::slice(sample(dplyr::n(),size=1,prob=pi))|>
      (`[`)(c("r2","q"))|>
      unlist()
}
## **Question 2 - (III)**

#'@examples
#'r2=.5;q=.5;u=matrix(0,nrow(x),1);x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,95));sigma_epsilon=1;beta=rnorm(z)*z;y=generate_single_y(xx=x)
#'likelihood_z_cond_y_u_x_phi_gamma(z,y,u,x,phi=0,r2,q,tt=nrow(xx),s_z=sum(z),k=ncol(x),barvx=1,gamma2=gamma2_f(r2,k,q,barvx))
loglikelihood_z_cond_y_u_x_phi_gamma<-
  function(z,tilde_y,
           ttildeytildey,
           x,q,tt,k,
           gamma2){
    s_z=sum(z)
    if(s_z==0){-Inf}else{{
    tilde_x=x[,z==1]
    tilde_w=(t(tilde_x)%*%tilde_x+diag(s_z)/gamma2);
    hat_tilde_beta=solve(tilde_w)%*%t(tilde_x)%*%tilde_y
    log((q^s_z)*((1-q)^(k-s_z)))+
      (-s_z/2)*log(gamma2)+
      (-1/2)*log(det(tilde_w))+
      (-tt/2)*log((ttildeytildey-t(hat_tilde_beta)%*%tilde_w%*%hat_tilde_beta)/2)}}
      #log(gamma(tt/2))
  }

#'@examples
#'r2=.5;q=.5;u=generate_u();x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,95));sigma_epsilon=1;beta=rnorm(z)*z;y=generate_single_y(xx=x)
#'sample_zi_cond_zj_y_u_x_phi_gamma(z,2,y,u=generate_u(),x,phi=0,r2,q,tt=nrow(x),k=ncol(x),barvx=1,gamma2=gamma2_f(r2 = r2,k =k,q = q,barvx = barvx))

sample_zi_cond_zj_y_u_x_phi_gamma<-
  function(z,i,tilde_y,ttildeytildey,
           x,q,tt,k,
           gamma2){
    logprob=c(loglikelihood_z_cond_y_u_x_phi_gamma((`[<-`)(z,i,0),tilde_y=tilde_y,ttildeytildey=ttildeytildey,
                                                x=x,q=q,tt=tt,k=k,
                                                gamma2=gamma2),
           loglikelihood_z_cond_y_u_x_phi_gamma((`[<-`)(z,i,1),tilde_y=tilde_y,ttildeytildey=ttildeytildey,
                                                x=x,q=q,tt=tt,k=k,
                                                gamma2=gamma2))|>
             (function(x){x-max(x)})()
           
    if(max(logprob,na.rm=TRUE)>-Inf){1}else{
    sample(0:1,
           size=1,
           prob=logprob|>exp())}
    }

#'@examples
#'r2=.5;q=.5;u=generate_u();x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,95));sigma_epsilon=1;beta=rnorm(z)*z;y=generate_single_y(xx=x)
#'sample_z_cond_y_u_x_phi_gamma(z,y,u=generate_u(),x,phi=0,r2,q,tt=nrow(x),k=ncol(x),barvx=1,gamma2=gamma2_f(r2 = r2,k =k,q = q,barvx = barvx))

sample_z_cond_y_u_x_phi_gamma<-
  function(z,tilde_y,ttildeytildey,
           x,q,tt,k,
           gamma2){
    for(i in 1:k){
      z[i]=
      sample_zi_cond_zj_y_u_x_phi_gamma(z=z,i=i,tilde_y=tilde_y,ttildeytildey=ttildeytildey,
                                        x=x,q=q,tt=tt,k=k,
                                        gamma2=gamma2)}
    z}



### IV.

#'@examples
#'r2=.5;q=.5;u=generate_u();x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,95));sigma_epsilon=1;beta=rnorm(z)*z;y=generate_single_y(xx=x)
#'loglikelihood_r2_q_cond_y_u_x_theta_z_a_b_aa_bb(r2,q,y,u,x,sigma_epsilon,beta,z)

sample_sigmaepsilon_cond_y_u_x_phi_r2_q_z<-
  function(z,tt,ttildeytildey,tilde_w,hat_tilde_beta){
    (1/rgamma(
      1, 
      shape=tt/2,  
      rate =(ttildeytildey-t(hat_tilde_beta)%*%tilde_w%*%hat_tilde_beta)/2 ))|>sqrt()
  }



### V.

#'@examples
#'r2=.5;q=.5;u=generate_u();x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'z=rep(c(1,0),c(5,95));sigma_epsilon=1;beta=rnorm(z)*z;y=generate_single_y(xx=x)
#'loglikelihood_r2_q_cond_y_u_x_theta_z_a_b_aa_bb(r2,q,y,u,x,sigma_epsilon,beta,z)
sample_tildebeta_cond_y_u_x_phi_r2_q_z_sigma2<-
  function(sigma_epsilon,
           tilde_y,
           tilde_x,
           tilde_w){
    inv_tilde_w<-solve(tilde_w)
    MASS::mvrnorm(n = 1, 
            mu = inv_tilde_w%*%t(tilde_x)%*%(tilde_y), 
            Sigma = (sigma_epsilon^2)*inv_tilde_w)
    }

#0. Initial values
initial_values_f<-function(y,x){
  lasso_model <- glmnet::glmnet(x, y, lambda = .1, intercept = FALSE,family="gaussian")
  beta<-lasso_model$beta|>as.vector()
  tilde_x=x[,beta!=0]
  ((y-x%*%beta)^2)|>sum()|>(`/`)(nrow(x)-sum(beta!=0))|>sqrt()->sigma_epsilon
  #(lm(y~tilde_x+0)|>summary())$sigma
  return(list(beta=beta,sigma_epsilon=sigma_epsilon))}

#'@examples
#'tt=default_tt
#'k=default_k
#'x=generate_multiple_x(number_of_datasets = 1)[[1]]
#'y=generate_single_y(s=50,r_y=.5,xx=x)
#'u=generate_u()
#'a=default_a
#'b=default_a
#'aa=default_b
#'bb=default_bb
#'barvx=default_barvx
#'phi=0
#'r2_q_grid=r2_q_grid_generate()
#'nrep=default_nrep
#'burning=default_burning
Gibbs_q<-function(x,
                  barvx=barvx_f(x),
                  y,
                  u,
                  a,
                  b,
                  aa,
                  bb,
                  tt=nrow(x),
                  k=ncol(x),
                  phi,
                  r2_q_grid,
                  nrep,
                  burning){
  the_sample<-matrix(ncol=4)
  #Initialise
  initial_values_f(x=x,y=y)->initial_values
  beta=initial_values$beta
  z=beta!=0
  s_z=sum(z)
  sigma_epsilon=  initial_values$sigma_epsilon
  tilde_beta=beta[z==1]
  tbetabeta=t(tilde_beta)%*%tilde_beta
  tilde_y=y-u%*%phi
  ttildeytildey=t(tilde_y)%*%tilde_y
  #I.
  for(i in 1:nrep){
  r2_q=sample_r2_q_cond_y_u_x_theta_z(sigma_epsilon=sigma_epsilon,barvx=barvx,k=k,a=a,b=b,aa=aa,bb=bb,tbetabeta=tbetabeta,s_z=s_z,r2_q_grid=r2_q_grid)
  r2=r2_q["r2"]
  q=r2_q["q"]
  gamma2=gamma2_f(r2,k,q,barvx)
#II.
  if(FALSE){
  phi=sample_phi()
  tilde_y=y-u%*%phi
  ttildeytildey=t(tilde_y)%*%tilde_y}
#III.
  z=sample_z_cond_y_u_x_phi_gamma(z=z,tilde_y=tilde_y,ttildeytildey = ttildeytildey,
                                  x=x,q=q,tt=tt,k=k,
                                  gamma2=gamma2)
  s_z=sum(z)
  tilde_x=x[,z==1]
  tilde_w=(t(tilde_x)%*%tilde_x+diag(s_z)/gamma2)
  inv_tilde_w=if(s_z==0){0}else{solve(tilde_w)}
  hat_tilde_beta=inv_tilde_w%*%t(tilde_x)%*%(tilde_y)
  #IV.
  sigma_epsilon=sample_sigmaepsilon_cond_y_u_x_phi_r2_q_z(
    tt=tt,ttildeytildey=ttildeytildey,tilde_w=tilde_w,
    hat_tilde_beta=hat_tilde_beta)    
  #V.
  tilde_beta=sample_tildebeta_cond_y_u_x_phi_r2_q_z_sigma2(
    sigma_epsilon=sigma_epsilon,
    tilde_y=tilde_y,
    tilde_x=tilde_x,
    tilde_w=tilde_w)

  if(i%%1000==0){print(paste0(Sys.time(),i))}
  the_sample<-rbind(the_sample,c(q,s_z,r2,sigma_epsilon))
  }
  the_sample
}



plot_q_1_f<-function(q,tt){
  require(ggplot2)
  q[,,1,,drop=TRUE]|>
    as.data.frame.table(responseName = "q")|>
    dplyr::mutate(j=strtoi(j),
                  s=strtoi(levels(s)[s]))->xx
  xx[1:90,]|>
    ggplot(aes(x=j,y=q))+
    geom_line()+
    facet_grid(s~r_y)
  xx|>
    ggplot(aes(x=q,xintercept=s/tt))+
    geom_histogram()+
    facet_grid(s~r_y)+
    geom_vline(mapping = aes(xintercept=s/tt),color="red")
  
}
```
