
Gibbs_q2<-function(x,
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
                  burning,
                  beta,
                  sigma_epsilon,
                  testgibbs=FALSE){
  if(testgibbs){the_sample<-matrix(ncol=4)}
  qq=vector()
  ((y-x%*%beta)^2)|>sum()|>(`/`)(nrow(x)-sum(beta!=0))|>sqrt()->sigma_epsilon
  z=beta!=0
  s_z=sum(z)
  tilde_beta=beta[z==1]
  tbetabeta=sum(tilde_beta^2)
  tilde_y=y
  ttildeytildey=sum(tilde_y^2)
  for(i in 1:(nrep+burning)){
    #I.
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
    z=sample_z_cond_y_u_x_phi_gamma(z=z,
                                    tilde_y=tilde_y,
                                    ttildeytildey = ttildeytildey,
                                    x=x,
                                    q=q,
                                    tt=tt,
                                    k=k,
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
    tbetabeta=sum(tilde_beta^2)
    if(i%%1000==0){print(paste0(Sys.time(),i))}
    if(testgibbs){the_sample<-rbind(the_sample,c(q,s_z,r2,sigma_epsilon))}else{qq=c(qq,q)}
    
  }
  if(testgibbs){the_sample}else{qq}}





copy_source_files_on_server<-function(
  cred =montruc::UCAMtruc(),
  source_code_dir_on_server="~/Bayes3/R"){
  "R"|>list.files(full.names = TRUE)|>sapply(FUN = source)
  list_of_files_to_sync<-
    c(file.path("R",list.files("R")))
  list_of_files_to_sync|>writeLines(con="list_of_files_to_sync.txt")
  system(paste0("sshpass -p ",montruc::UCAMtruc()$password," rsync  -av -q -e ssh  --files-from='list_of_files_to_sync.txt' . "," ",
                montruc::UCAMtruc()$username,"@camclust:",source_code_dir_on_server))}




sendtocluster<-function(x,
             s,
             r_y,
             a,
             b,
             aa,
             bb,
             r2_q_grid,
             nrep,
             burning,
             method="Daniel",
             test=FALSE){
  dd<-simulation_parameters(s,r_y,length(x))
  scriptparams=
        list(dd=dd,x=x,
             s=s,
             r_y=r_y,
             a=a,
             b=b,
             aa=aa,
             bb=bb,
             r2_q_grid=r2_q_grid,
             nrep=nrep,
             burning=burning,
             method=method)
      
      sendscripttoserver(
        .expression=expression(
          {set.seed(taskid)
            dd[taskid,]->d
            xx<-x[[d$i]]
            tt=nrow(xx)
            u=generate_u(tt)
            phi=0
            k=ncol(xx)
            tilde_beta=draw_tildebeta(d$s)
            beta=c(tilde_beta,rep(0,k-d$s))
            xbeta=xx[,1:d$s]%*%tilde_beta
            #xbeta=xx[,sample(ncol(xx),size = d$s,replace=FALSE)]%*%beta#we could take the n first one would do the same
            sigma_epsilon<-sqrt((1/r_y-1)*mean(xbeta^2))
            y=xbeta +rnorm(nrow(xx),sd=sigma_epsilon)
            barvx=barvx_f(xx)
            
            if(method=="Daniel"){
              gibbs_sample<-Gibbs_q(x = xx,
                          y = y,
                          u = u,
                          barvx = barvx,
                          tt = tt,k = k,phi = phi,r2_q_grid = r2_q_grid,a=a,
                          b=b,
                          aa=aa,
                          bb=bb,
                          nrep=nrep,
                          burning=burning)}
            
            if(method=="Augustin"){
              gibbs_sample<-Gibbs(N=6000,a=a,A=aa,b=b,B=bb,k=k,U=0,phi=0,X=xx,Y=y)}
            
            gibbs_sample
            }),
        scriptparams=scriptparams,
        stamp=paste0("~/Bayes3/",method,"/q"), 
        time="4:00:00",
        .array = if(test){1}else{nrow(dd)},
        mem=NULL,
        mem_per_cpu="3380mb",
        source_code_dir_on_server="~/Bayes3/R")
      
}

 
  
  

  
  

merge_all_qs_on_server<-function(stamp="~/Bayes3/Daniel/allq",directory_outputs="~/Bayes3/Daniel"){
  sendscripttoserver(
    .expression=expression(
      {    
        simulation_parameters(default_s,default_r_y,default_number_of_datasets)|>
          dplyr::mutate(taskid=dplyr::row_number())|>
          dplyr::mutate(file=paste0(directory_outputs,"/q",taskid,"-output.rda"))|>
          dplyr::filter(file.exists(file))|>
          plyr::ddply(~s+r_y+i,function(d){
            data.frame(q=d$file|>load()|>get())})
      }),
    scriptparams=list(directory_outputs=directory_outputs,default_s=default_s,default_r_y=default_r_y),
    stamp=stamp, 
    time="0:05:00",
    .array = 0,
    mem=NULL,
    mem_per_cpu="3380mb",
    source_code_dir_on_server="~/Bayes3/R")}









sendtocluster_assignment_3<-function(dataset,
                                     r2_q_grid,
                                     nrep,
                                     burning,
                                     nchain,
                                     test){
  copy_source_files_on_server(source_code_dir_on_server="~/Bayes3")
  xx<-dataset|>dplyr::select(-Outcome,-intercept)|>
    as.matrix()|>
    plyr::aaply(2,function(x){(x-mean(x))/sd(x)})|>
    t()
  print("")
  y=dataset[["Outcome"]]|>as.matrix()
  u=dataset|>as.matrix()|>(`[`)(,"intercept",drop=FALSE)
  tt=nrow(xx)
  k=ncol(xx)
  
  barvx=barvx_f(dataset)
  
  scriptparams=list(xx=xx,
                    y=y,u=u,tt=tt,k=k,barvx = barvx,
                    nrep=nrep,
                    burning=burning,
                    posterior="beta",
                    r2_q_grid=r2_q_grid)
  sendscripttoserver(
    .expression=expression(
      {set.seed(taskid)
        Gibbs_q(x = xx,y = y,u = u,barvx = barvx,
                tt = tt,k = k,phi = 0,r2_q_grid = r2_q_grid,a=1,b=1,aa=1,
                bb=1,nrep=nrep,burning=burning,posterior="beta")}),
    scriptparams=scriptparams,
    stamp=paste0("~/Bayes3/q"), 
    time="4:00:00",
    .array = if(test){2}else{nchain},
    mem=NULL,
    mem_per_cpu="3380mb",
    source_code_dir_on_server="~/Bayes3/R")
  
}

merge_cluster_assignment_3<-function(sample_on_cluster){

  if(is.null(sample_on_cluster)){dir="~/group_space/epidem-readonly-data/Bayes3"
  dirorfilepathonserver=file.path(dir,paste0("q",1:100,"-output.rda"))
  existingfiles=file.path(dir,list.files(dir))
  sample_on_cluster=intersect(dirorfilepathonserver,existingfiles)}
  
  scriptparams=list(sample_on_cluster=sample_on_cluster)
  sendscripttoserver(
    .expression=expression(
      {set.seed(taskid)
        sample_on_cluster|>
          lapply(load)|>
          lapply(get)|>
          c(list(along=4))|>
          do.call(abind::abind)}),
    scriptparams=scriptparams,
    stamp=paste0("~/Bayes3/qall"), 
    time="4:00:00",
    mem=NULL,
    mem_per_cpu="3380mb",
    source_code_dir_on_server="~/Bayes3/R")
  
  
  
  }

get_cluster_assignment_3<-function(merged_cluster_samples){
  if(is.null(merged_cluster_samples)){merged_cluster_samples=paste0("~/group_space/epidem-readonly-data/Bayes3/qall-output.rda")}
  get_data_from_server(merged_cluster_samples)|>load()|>get()
}














