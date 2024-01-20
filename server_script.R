
"R"|>list.files(full.names = TRUE)|>sapply(FUN = source)
library(targets)
tar_load(x)
copy_source_files_on_server(  source_code_dir_on_server="~/Bayes3")
sendtocluster(x,
              s=default_s,
              r_y=default_r_y,
              a=default_a,
              b=default_b,
              aa=default_aa,
              bb=default_bb,
              r2_q_grid=r2_q_grid_generate(),
              nrep=default_nrep,
              burning=default_burning,
              method="Daniel",
              test=FALSE)

sendtocluster(x[1],
              s=default_s,
              r_y=default_r_y,
              a=default_a,
              b=default_b,
              aa=default_aa,
              bb=default_bb,
              r2_q_grid=r2_q_grid_generate(),
              nrep=50000,
              burning=default_burning,
              method="alot",
              test=FALSE)






# merge_all_qs_on_server("~/Bayes3/allq")
qs<-get_data_from_server("~/Bayes3/Daniel/allq-output.rda")|>load()|>get()|>
  dplyr::select(-q)|>
  dplyr::rename(q=q.q,
                r2=q.r2, 
                s_z=q.s_z,
                sigma_epsilon= q.sigma_epsilon)

  save(qs,file="samples_of_q.rda")

  
  
  
  merge_all_qs_on_server(stamp="~/Bayes3/alot/allq",dir="~/Bayes3/alot")
  
  get_data_from_server("~/Bayes3/Daniel/allq-output.rda")|>load()|>get()|>dplyr::select(-q)|>
    dplyr::rename(q=q.q,
                  r2=q.r2, 
                  s_z=q.s_z,
                  sigma_epsilon= q.sigma_epsilon)->qs50000
  save(qs,file="samples_of_q50000.rda")
  