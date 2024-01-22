cluster=file.exists("R/play_with_server.R")
#Turn to FALSE if not cluster
# Load packages required to define the pipeline:
library(targets)
library(MASS)
library(glmnet)#for lasso

options(clustermq.scheduler = "multicore")

# Run the R scripts in the R/ folder with our custom functions:
"R"|>list.files(full.names = TRUE)|>sapply(FUN = source)
try(copy_source_files_on_server(),silent = TRUE)

library(tarchetypes)
renv::dependencies("./R")$Package|>unique()->required_packages
required_packages|>writeLines("requirements.txt")
setdiff(required_packages,installed.packages()[,"Package"])|>sapply(install.packages)
library(ggplot2)
library(R2jags)

#model parameters - hyperpriors
a    =default_a
b    =default_b
aa   =default_aa
bb   =default_bb

#simulations tuning parameters
nrep            =100000
burning         =1000;default_burning
nchain=100
# Targets list:
list(
  #question 1. Generate X
  tar_target(
    name = assignment_3_data,
    command = final_project_get_data()),
  tar_target(
    name = data_dictionnary,
    command = final_project_get_dictionnary(
      names(assignment_3_data))),
  tar_target(
    name = r2_q_grid,
    command = r2_q_grid_generate()),
 tar_target(plot_1_data,plot_1_data_f(assignment_3_data,r2_q_grid)),
 tar_target(plot_1,assign3_plot1_f(plot_1_data,"output/assignment_3_fig1.pdf")),
 tar_target(
    name = sample_on_cluster,
    command = if(cluster){sendtocluster_assignment_3(dataset=assignment_3_data,
                                         r2_q_grid=r2_q_grid,
                                         nrep=nrep,
                                         burning=burning,
                                         nchain=nchain,
                                         test=FALSE)}),
 tar_target(
    name = merged_cluster_samples,
    command = if(cluster){merge_cluster_assignment_3(sample_on_cluster)}),
 tar_target(
    name = merged_samples,
    command = if(cluster){get_cluster_assignment_3(merged_cluster_samples)}else{
      final_project_run_mcmc(dataset=assignment_3_data,
                             r2_q_grid=r2_q_grid,
                             nrep=nrep,
                             burning=burning,nchain=nchain)}),
 tar_target(
   name = assign3_plot_acf,
   command = assign3_plot_acf_f(merged_samples,output = "output/assignment_3_fig2.pdf")),
 tar_target(
   name = assign3_plot_trace,
   command = assign3_plot_trace_f(merged_samples,output = "output/assignment_3_figtrace.pdf")),
 tar_target(
   name = assign3_plot_ci,
   command = assign3_plot_ci_f(merged_samples,r2_q_grid =r2_q_grid, output = "output/assignment_3_figCI.pdf")))

