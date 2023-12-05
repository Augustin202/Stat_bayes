library(MASS)# for multivariate normal
library(glmnet)#for lasso


# Run the R scripts in the R/ folder with our custom functions:
"R"|>list.files(full.names = TRUE)|>sapply(FUN = source)

#model parameters - constants
k     =default_k
tt    =default_tt

#model parameters - hyperpriors
a    =default_a
b    =default_b
aa   =default_aa
bb   =default_bb

#simulations tuning parameters
rho   =default_rho
number_of_datasets  =default_number_of_datasets
s     =default_s
r_y   =default_r_y
nrep            =default_nrep
burning         =default_burning

# Targets list:
    x=generate_multiple_x(k=k,tt=tt,rho = rho,number_of_datasets = number_of_datasets)
  #question 2. Generate q 
    q=generate_y_sample_q(s = s,
                          r_y = r_y,
                          a = a,
                          b = b,
                          aa = aa,
                          bb = bb,
                          nrep = nrep,
                          burning = burning,
                          x = x)
  #question 3. plots 
    plot_q_1=plot_q_1_f(q,k,burning)
    plot_q_2=plot_q_2_f(q,k,burning)

