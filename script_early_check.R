
"R"|>list.files(full.names = TRUE)|>sapply(FUN = source)

copy_source_files_on_server(  source_code_dir_on_server="~/Bayes3")

merge_all_qs_on_server()

get_data_from_server("~/Bayes3/allq-output.rda")|>load()|>get()->qs

qs|>
  plot_q_1_f(burning=1000)

burning=1000
qs|>dplyr::group_by(s,r_y,i)|>
  dplyr::mutate(t=dplyr::row_number())|>
  dplyr::filter(dplyr::row_number()>burning)|>
  dplyr::ungroup()|>
  dplyr::group_by(s,r_y,i)|>
  dplyr::mutate(Mq=median(q,na.rm=TRUE),Eq=mean(q,na.rm=TRUE),diff=abs(s/default_k-Eq))|>
  dplyr::ungroup()|>
  dplyr::group_by(s,r_y)|>
  dplyr::filter(diff==min(diff))|>
  dplyr::ungroup()|>
  ggplot(mapping = aes(x=q,group=interaction(s,r_y)))+
  geom_histogram(aes(y=..density..),color="black",alpha=.5)+
  geom_density()+
  geom_vline(mapping = aes(xintercept=s/default_k),color="red")+
  geom_vline(mapping = aes(xintercept=Mq),color="green")+
  geom_vline(mapping = aes(xintercept=Eq),color="blue")+
  facet_grid(s~r_y)
  

qs|>dplyr::group_by(s,r_y,i)|>
  dplyr::mutate(t=dplyr::row_number())|>
  dplyr::filter(dplyr::row_number()>burning)|>
  dplyr::ungroup()|>
  dplyr::group_by(s,r_y,i)|>
  dplyr::mutate(Mq=median(q,na.rm=TRUE),Eq=mean(q,na.rm=TRUE),diff=abs(s/default_k-Eq))|>
  dplyr::ungroup()|>
  dplyr::group_by(s,r_y)|>
  dplyr::filter(diff==max(diff))|>
  dplyr::ungroup()|>
  ggplot(mapping = aes(x=q,group=interaction(s,r_y)))+
  geom_histogram(aes(y=..density..),color="black",alpha=.5)+
  geom_density()+
  geom_vline(mapping = aes(xintercept=s/default_k),color="red")+
  geom_vline(mapping = aes(xintercept=Mq),color="green")+
  geom_vline(mapping = aes(xintercept=Eq),color="blue")+
  facet_grid(s~r_y)


qs|>dplyr::group_by(s,r_y,i)|>
  dplyr::mutate(t=dplyr::row_number())|>
  dplyr::filter(dplyr::row_number()>burning)|>
  dplyr::ungroup()|>
  dplyr::group_by(s,r_y,i)|>
  dplyr::mutate(Mq=median(q,na.rm=TRUE),Eq=mean(q,na.rm=TRUE),diff=abs(s/default_k-Eq))|>
  dplyr::ungroup()|>
  dplyr::group_by(s,r_y)|>
  dplyr::filter(i==1)|>
  dplyr::ungroup()|>
  ggplot(mapping = aes(x=q,group=interaction(s,r_y)))+
  geom_histogram(aes(y=..density..),color="black",alpha=.5)+
  geom_density()+
  geom_vline(mapping = aes(xintercept=s/default_k),color="red")+
  geom_vline(mapping = aes(xintercept=Mq),color="green")+
  geom_vline(mapping = aes(xintercept=Eq),color="blue")+
  facet_grid(s~r_y)


aa|>dplyr::group_by(s,r_y)|>
  dplyr::slice(which.min(abs(s/default_k-Eq)))|>
  dplyr::ungroup()->bb

bb|>plyr::dlply(~s+r_y,function(d){
  plot_q_2_f(qs|>dplyr::filter(s==d$s,r_y==d$r_y,i==d$i),s=d$s,r_y=d$r_y,i=d$i,burning=1000)})->bestplots
  
burning=1000
qs|>dplyr::group_by(s,r_y,i)|>
  dplyr::mutate(t=dplyr::row_number(),group=t%/%1000)|>
  dplyr::filter(dplyr::row_number()>burning)|>
  dplyr::ungroup()|>
  dplyr::group_by(s,r_y,i,group)|>
  dplyr::summarise(Eq=mean(q,na.rm=TRUE))|>
  dplyr::ungroup()|>
  dplyr::group_by(s,r_y)|>
  dplyr::mutate(meanEq=mean(Eq,na.rm=TRUE))|>
  dplyr::ungroup()|>
  ggplot(mapping = aes(x=Eq,fill=group,colour=group,group=group))+
  geom_histogram()+
  #    geom_histogram(aes(y=..density..),alpha=.5)+
  #    geom_density()+
  facet_grid(s~r_y)+
  geom_vline(mapping = aes(xintercept=s/k),color="red")+
  geom_vline(mapping = aes(xintercept=meanEq),color="blue")
qs|>  plot_q_2_f(burning=5000)

qs|>dplyr::filter(s==100,r_y==.5,i==1)|>
  dplyr::mutate(ii=dplyr::row_number())|>
  ggplot(aes(x=ii,y=q))+geom_line()
