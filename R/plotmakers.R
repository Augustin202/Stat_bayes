
plot_q_1_f<-function(qs,k=default_k,burning){
  require(ggplot2)
  qs|>
    dplyr::group_by(s,r_y,i)|>
    dplyr::filter(dplyr::row_number()>burning)|>
    dplyr::summarise(Mq=median(q,na.rm=TRUE))|>
    dplyr::ungroup()|>
    dplyr::group_by(s,r_y)|>
    dplyr::mutate(medianMq=median(Mq,na.rm=TRUE),)|>
    dplyr::ungroup()|>
    ggplot(mapping = aes(x=Mq))+
    geom_histogram()+
    #    geom_histogram(aes(y=..density..),alpha=.5)+
    #    geom_density()+
    facet_grid(s~r_y)+
    geom_vline(mapping = aes(xintercept=s/k),color="red")+
    geom_vline(mapping = aes(xintercept=medianMq),color="green")
}




plot_q_2_f<-function(q,
                     k=default_k,burning,
                     i=1,s=5,r_y=.02,a=default_a,b=default_b){
  require(ggplot2)
  q|>
    dplyr::filter(s==s,r_y==r_y,i==i,dplyr::row_number()>burning)|>
    dplyr::mutate(Eq=mean(q))|>
    ggplot(mapping = aes(x=q))+
    geom_histogram(aes(y=..density..),color="black",alpha=.5)+
    geom_line(data=expand.grid(s=default_s,r_y=default_r_y,q=seq(0,1,by=.01))|>
                dplyr::mutate(y=dbeta(q,default_a,default_b)),
              aes(x=q,y=y),linetype = "dashed")+
    geom_vline(mapping = aes(xintercept=s/k),color="red")+
    geom_vline(mapping = aes(xintercept=Eq),color="blue")}






plot_q_3_f_v_x<-function(qs,k=default_k,burning){
  require(ggplot2)
  qs|>
    dplyr::group_by(s,r_y,i)|>
    dplyr::filter(dplyr::row_number()>burning)|>
    dplyr::summarise(Eq=mean(q,na.rm=TRUE),
                     Mq=median(q,na.rm=TRUE))|>
    dplyr::ungroup()|>
    dplyr::group_by(s,r_y)|>
    dplyr::mutate(meanEq=mean(Eq,na.rm=TRUE),medianMq=median(Mq,na.rm=TRUE),)|>
    dplyr::ungroup()|>
    ggplot(mapping = aes(x=Eq))+
    geom_histogram()+
    #    geom_histogram(aes(y=..density..),alpha=.5)+
    #    geom_density()+
    facet_grid(s~r_y)+
    geom_vline(mapping = aes(xintercept=s/k),color="red")+
    geom_vline(mapping = aes(xintercept=medianMq),color="green")+
    geom_vline(mapping = aes(xintercept=meanEq),color="blue")
}



plot_best_scenario<-function(qs){
  burning=1000
  qs|>dplyr::group_by(s,r_y,i)|>
    dplyr::mutate(t=dplyr::row_number())|>
    dplyr::filter(dplyr::row_number()>burning)|>
    dplyr::ungroup()|>
    dplyr::group_by(s,r_y,i)|>
    dplyr::mutate(Mq=median(q,na.rm=TRUE),
                  Eq=mean(q,na.rm=TRUE),
                  diff=abs(s/default_k-Eq))|>
    dplyr::ungroup()|>
    dplyr::group_by(s,r_y)|>
    dplyr::filter(diff==min(diff,na.rm=TRUE))|>
    dplyr::ungroup()|>
    ggplot(mapping = aes(x=q,group=interaction(s,r_y)))+
    geom_histogram(aes(y=..density..),color="black",alpha=.5)+
    geom_density()+
    geom_line(data=expand.grid(s=default_s,r_y=default_r_y,q=seq(0,1,by=.01))|>
                dplyr::mutate(y=dbeta(q,default_a,default_b)),
              aes(x=q,y=y),linetype = "dashed")+
    geom_vline(mapping = aes(xintercept=s/default_k),color="red")+
    geom_vline(mapping = aes(xintercept=Mq),color="green")+
    geom_vline(mapping = aes(xintercept=Eq),color="blue")+
    facet_grid(s~r_y)
}



plot_random_scenario<-function(qs){
  burning=1000
  qs|>
    dplyr::filter(i==9)|>
    dplyr::group_by(s,r_y)|>
    dplyr::mutate(t=dplyr::row_number())|>
    dplyr::filter(dplyr::row_number()>burning)|>
    dplyr::mutate(Mq=median(q,na.rm=TRUE),
                  Eq=mean(q,na.rm=TRUE))|>
    dplyr::ungroup()|>
    ggplot(mapping = aes(x=q,group=interaction(s,r_y)))+
    geom_histogram(aes(y=..density..),color="black",alpha=.5)+
    geom_density()+
    geom_line(data=expand.grid(s=default_s,r_y=default_r_y,q=seq(0,1,by=.01))|>
                dplyr::mutate(y=dbeta(q,default_a,default_b)),
              aes(x=q,y=y),linetype = "dashed")+
    geom_vline(mapping = aes(xintercept=s/default_k),color="red")+
    geom_vline(mapping = aes(xintercept=Mq),color="green")+
    geom_vline(mapping = aes(xintercept=Eq),color="blue")+
    facet_grid(s~r_y)
}



plot_random_scenario_s_z<-function(qs){
  burning=1000
  qs|>
    dplyr::filter(i==25)|>
    dplyr::group_by(s,r_y)|>
    dplyr::mutate(t=dplyr::row_number())|>
    dplyr::filter(dplyr::row_number()>burning)|>
    dplyr::mutate(Mq=median(s_z,na.rm=TRUE),
                  Eq=mean(s_z,na.rm=TRUE))|>
    dplyr::ungroup()|>
    ggplot(mapping = aes(x=s_z,group=interaction(s,r_y)))+
    geom_histogram(aes(y=..density..),color="black",alpha=.5)+
    geom_density()+
    geom_line(data=expand.grid(s=default_s,r_y=default_r_y,s_z=0:default_k)|>
                dplyr::mutate(y=VGAM::dbetabinom.ab(s_z,default_k,default_a,default_b)),
              aes(x=s_z,y=y),linetype = "dashed")+
    geom_vline(mapping = aes(xintercept=s),color="red")+
    geom_vline(mapping = aes(xintercept=Mq),color="green")+
    geom_vline(mapping = aes(xintercept=Eq),color="blue")+
    facet_grid(s~r_y, scales="free")
}


plot_trace_random_scenario<-function(qs){
  burning=1000
  qs|>
    dplyr::filter(i==9)|>
    dplyr::group_by(s,r_y)|>
    dplyr::mutate(t=dplyr::row_number())|>
    dplyr::ungroup()|>
    ggplot(mapping = aes(x=t,y=q,group=interaction(s,r_y)))+
    geom_line(linewidth=.1)+
    geom_vline(mapping = aes(xintercept=burning),color="red")+
    facet_grid(s~r_y)
}
