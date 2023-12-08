
plot_q_1_f<-function(qs,k=default_k,burning){
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

plot_q_1_f_v_x<-function(qs,k=default_k,burning){
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




plot_q_2_f<-function(q,
                     k=default_k,burning,
                     i=1,s=5,r_y=.02){
  require(ggplot2)
  q|>
    dplyr::filter(s==s,r_y==r_y,i==i,dplyr::row_number()>burning)|>
    dplyr::mutate(Eq=mean(q))|>
    ggplot(mapping = aes(x=q))+
    geom_histogram(aes(y=..density..),color="black",alpha=.5)+
    geom_density()+
    geom_vline(mapping = aes(xintercept=s/k),color="red")+
    geom_vline(mapping = aes(xintercept=Eq),color="blue")}

