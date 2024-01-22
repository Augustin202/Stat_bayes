#'@examples
#'tar_load(assignment_3_data)
#'datasetnames=names(assignment_3_data)


final_project_get_dictionnary<-function(datasetnames){
  if(!file.exists("extdata/dict.txt")){
    dir.create("extdata",showWarnings = FALSE)
    x=tempfile(fileext = ".txt")
    "https://www2.nber.org/pub/barro.lee/readme.txt"|>
      download.file(destfile=
                      "extdata/dict.txt")}
  "extdata/dict.txt"|>
    readLines()|>
    (function(x){x[x!=""]})()|>
    (`[`)(77:314)->dict0
  pos1<-plyr::aaply(dict0,1,function(text,pattern){gregexpr(pattern,text)[[1]][1]},pattern = ":")
  starts<-which(pos1==18)
  grep(":",dict0)
  ends=c(starts[-1]-1,length(dict0))
  dictvar0=sapply(1:length(starts),
                  function(i){
                    paste(trimws(dict0,which="right")[starts[i]:ends[i]],collapse=" ")})|>
    gsub(pattern = "Source:",replacement = "Source :")
  hassource=grepl("source ",tolower(dictvar0))  
  isy<-grepl("xx",dictvar0)
  is5y<-(!isy)&grepl("x",dictvar0)
  dictvar=data.frame(
    nom=substr(dictvar0,1,17)|>trimws()|>
      gsub(pattern="xx",replacement="65")|>
      gsub(pattern="x",replacement="1")|>
      tolower(),
    description=substr(dictvar0,19,nchar(dictvar0))|>
      trimws()|>
      sub(pattern = "Source :.*",replacement = "")|>
      sub(pattern = "\\(xx.*",replacement = ""),
    source=ifelse(hassource,
                  substr(dictvar0,19,nchar(dictvar0))|>trimws()|>sub(pattern = ".*Source :",replacement = ""),"")|>
      sub(pattern = "\\(xx.*",replacement = ""),
    `time frame`=ifelse(isy,"Year 1965",ifelse(is5y, "60--64","")))
  w=grep(",",dictvar$nom)
  dictvar2<-dictvar[w,]
  dictvar$nom[w]<-sub("\\,.*", "", dictvar$nom[w])|>trimws()
  dictvar2$nom<-sub(".*\\,", "", dictvar2$nom)|>trimws()
  dictvar=rbind(dictvar,
                dictvar2,
                data.frame(nom=c("outcome","intercept"),
                           description=c("national growth rates in GDP per capita for the periods 1965-1975","Constant variable equal to 1"),
                           source="",
                           `time frame`=""))|>
    dplyr::mutate(nom=sub(" .*", "", nom),
                  description=stringr::str_squish(description))|>
    dplyr::filter(is.element(tolower(nom),tolower(datasetnames)))|>
    dplyr::arrange(2*(nom!="outcome")+(nom!="intercept"),nom)|>
    (`rownames<-`)(NULL)
  }


final_project_get_data<-function(){
  data(GrowthData,package="hdm")|>get()
 }







#'@examples
#'tar_load(assignment_3_data)
#'tar_load(ré_q_grid)
#'nrep=10
#'burning=4
#'dataset=assignment_3_data
#'nchain=2
#'final_project_run_mcmc(dataset,r2_q_grid,nrep,burning,2)->samples




final_project_run_mcmc<-function(dataset,
                                 r2_q_grid,
                                 nrep,
                                 burning,
                                 nchain,
                                 a=1,
                                 b=1,aa=1,bb=1){

  xx<-dataset|>dplyr::select(-Outcome,-intercept)|>
    as.matrix()|>
    plyr::aaply(2,function(x){(x-mean(x))/sd(x)})|>
    t()

  y=dataset[["Outcome"]]|>as.matrix()
  u=dataset|>as.matrix()|>(`[`)(,"intercept",drop=FALSE)
  tt=nrow(xx)
  k=ncol(xx)
  
  barvx=barvx_f(dataset)
  #'@examples
  #'attach(list(x = xx,y = y,u = matrix(0,length(y),1),barvx = barvx,
  #'  tt = tt,k = k,phi = 0,r2_q_grid = r2_q_grid,a=1,b=1,aa=1,bb=1,
  #'  nrep=nrep,burning=burning,posterior="beta"))
  
  library(parallel)
  numberOfCores <- detectCores()
  
  mclapply(as.list(1:nchain),function(i){set.seed(i)
  Gibbs_q(x = xx,y = y,u = u,barvx = barvx,
             tt = tt,k = k,phi = 0,r2_q_grid = r2_q_grid,a=a,b=b,aa=aa,
             bb=bb,nrep=nrep,burning=burning,posterior="beta")}, 
  mc.cores = numberOfCores)
}


final_project_analyse_mcmc_sample<-function(){}

#'@examples
#'tar_load(merged_samples)
final_project_analyse_mcmc_sample<-function(merged_samples){
  merged_samples[1000+(100*(0:1000)),,1:2]->XX

  XX|>reshape2::melt()->YY
  (XX!=0)|>plyr::aaply(2,mean)|>sort()

  XX|>plyr::aaply(2,mean)
  
    
  YY|>
    dplyr::filter(Var2=="s_z")|>
    ggplot(aes(x=Var1,xend=Var1,y=0,yend=value,group=Var3))+
    facet_grid(~Var3)+
    geom_segment()
    
XX|>
  aperm(c(1,3:2))|>
  (function(x){matrix(x,prod(dim(x)[1:2]),dim(x)[3])})()->x
x|>
  sign()|>
    var()|>
    (function(x){x/sqrt(t(t(c(diag(x))))%*%t(c(diag(x))))})()->
    corr
  corr|>
    reshape2::melt()|>
    ggplot(aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()
  
  
  YY|>
    dplyr::filter(Var2=="intercept",Var3==1, value!=0)|>
    ggplot(aes(x=value))+
    geom_density()
#corr|>heatmap()  

  YY|>dplyr::filter(is.element(Var3,c()))
    
}

#'@examples
#'tar_load(merged_samples)
assign3_plot_acf_f<-function(merged_samples,output="output/assignment_3_fig2.pdf"){
    merged_samples$onesample[1:50000,c("q","s_z","r2","sigma_epsilon","phi1")]|>
    plyr::aaply(2,function(x){acf(x,lag.max=100,plot=FALSE)$acf})|>
    (function(x){dimnames(x)[2]=list(lag=0:100);x})()|>
    reshape2::melt()|>
    dplyr::mutate(            X1=dplyr::recode(X1,"q"="$q$",
                                     "s_z"="$\\sum_{i=1}^k z_i$",
                                     "sigma_epsilon"="$\\sigma_\\varepsilon$",
                                     "r2"="$R^2$",
                                     "phi1"="$\\phi_1$"))|>
    dplyr::filter(Var2>1)|>
    ggplot(aes(x=Var2,xend=Var2,y=0,yend=value,group=X1))+
    facet_grid(rows = vars(X1),scales="free_y")+
    geom_segment()+xlab('')+ylab('')+
    geom_hline(yintercept = qnorm(c(.025,.975))/sqrt(50000),color="red")->the_plot

    if(is.element("SweaveLst",installed.packages())){
      SweaveLst::graph2pdffile(texte="print(the_plot)",output =  output,
                               widthe=8,heighte=4)
    }else{
      ggsave(plot = the_plot,
             filename = output,
             height=9,
             width=15,
             units="cm")
    }
  output
}

#'@examples
#'tar_load(merged_samples)
assign3_plot_trace_f<-function(merged_samples,output="output/assignment_3_fig3.pdf"){
  merged_samples$allsamples[,c("q","s_z","r2","sigma_epsilon","phi1"),1:3,drop=FALSE]|>
    reshape2::melt()|>
    dplyr::mutate(          Var2=dplyr::recode(Var2,"q"="$q$",
                                               "s_z"="$\\sum_{i=1}^k z_i$",
                                               "sigma_epsilon"="$\\sigma_\\varepsilon$",
                                               "r2"="$R^2$",
                                               "phi1"="$\\phi_1$"))|>
    ggplot(aes(x=Var1,y=value,group=Var3))+
    facet_grid(Var2~Var3,scale="free_y")+
    geom_line()+xlab('')+ylab('')->the_plot
  
  if(is.element("SweaveLst",installed.packages())){
    SweaveLst::graph2pdffile(texte="print(the_plot)",output =  output,
                             widthe=8,heighte=4)
  }else{
    ggsave(plot = the_plot,
           filename = output,
           height=9,
           width=15,
           units="cm")
  }
  output
}




plot_1_data_f<-function(assignment_3_data,r2_q_grid){
  set.seed(1)
  sample1=final_project_run_mcmc(dataset=assignment_3_data,
                                 r2_q_grid=r2_q_grid,
                                 nrep=1000,
                                 burning=0,
                                 nchain=1,
                                 a=1,
                                 b=1,
                                 aa=1,
                                 bb=1)  
  sample2=final_project_run_mcmc(dataset=assignment_3_data,
                                 r2_q_grid=r2_q_grid,
                                 nrep=1000,
                                 burning=0,
                                 nchain=1,
                                 a=1,
                                 b=10,
                                 aa=1,
                                 bb=1)
  plot_1_data=abind::abind(sample1=sample1[[1]],sample2=sample2[[1]],along=3)
}

#'@examples
#'tar_load(assignment_3_data)
#'tar_load(,r2_q_grid)
#'plot_1_data=plot1_data_f(assignment_3_data,r2_q_grid)
#'assign3_plot_1_f(plot_1_data,"output/assignment_3_fig1.pdf")

assign3_plot1_f<-function(plot_1_data,output= "output/assignment_3_fig1.pdf"){
  assign3_plot_1<-
    plot_1_data[,c("q","s_z","r2","sigma_epsilon","phi1"),]|>
    reshape2::melt()|>
    dplyr::mutate(Var3=dplyr::recode(Var3,"sample1"="$a=1, b=1$","sample2"="$a=10, b=100$"),
                  Var2=dplyr::recode(Var2,"q"="$q$",
                                     "s_z"="$\\sum_{i=1}^k z_i$",
                                     "sigma_epsilon"="$\\sigma_\\varepsilon$",
                                     "r2"="$R^2$",
                                     "phi1"="$\\phi_1$"))|>
    ggplot(aes(x=Var1,xend=Var1,y=0,yend=value,group=Var2))+
    facet_grid(rows = Var2~Var3,scales="free_y")+
    geom_line(aes(y=value))+xlab('')+ylab('')
  
  if(is.element("SweaveLst",installed.packages())){
    SweaveLst::graph2pdffile(texte="print(assign3_plot_1)",output = output,
                             widthe=8,heighte=4)
  }else{
    ggsave(plot = assign3_plot_1,
           filename =output,
           height=9,
           width=15,
           units="cm")  }
  output
  
}
#'@examples
#'tar_load(merged_samples)
#'tar_load(r2_q_grid)
#'output="output/assignment_3_figCI.pdf"
#'.vars=c("q","s_z","r2","sigma_epsilon","phi1")

assign3_plot_ci_f<-function(merged_samples,output="output/assignment_3_figCI.pdf",r2_q_grid,.vars=c("q","s_z","r2","sigma_epsilon","phi1")){
  merged_samples$allsamples[,.vars,]|>
    reshape2::melt()|>
    dplyr::select(-Var1,-Var3)|>
    dplyr::mutate(
      Var2=dplyr::recode(Var2,"q"="$q$",
                                     "s_z"="$\\sum_{i=1}^k z_i$",
                                     "sigma_epsilon"="$\\sigma_\\varepsilon$",
                                     "r2"="$R^2$",
                                     "phi1"="$\\phi_1$"))->X
    
    center<-function(x){
      y=unique(x)
      c(3/2*y[1]-y[2]/2,
        (y[-1]+y[-length(y)])/2,
        3/2*y[length(y)]-1/2*y[length(y)-1])|>sort()
      }
    geom.text.size = 6
    theme.size =  geom.text.size  
    
X|>ggplot(aes(x=value,group=Var2)) + 
  mapply(FUN = function(x){
    z=unique(x$Var2)
    if(is.element(z,c("$q$","$R^2$","$\\phi_1$"))){
      if(z=="$q$"){
        pp=geom_line(data=data.frame(x=seq(min(x$value),max(x$value),length.out=1000),Var2="$q$")|>
                       dplyr::mutate(y=dbeta(x,10,100)),aes(x=x,y=y),color="blue",size=.75,linetype = "dashed")
      }
      if(z=="$R^2$"){
        pp=geom_line(data=data.frame(x=seq(min(x$value),max(x$value),length.out=1000),Var2="$R^2$")|>
                       dplyr::mutate(y=1),aes(x=x,y=y),color="blue",size=.75,linetype = "dashed")
      }
      if(z=="$\\phi_1$"){
        pp=geom_line(data=data.frame(x=seq(min(x$value),max(x$value),length.out=1000),Var2="$\\phi_1$")|>
                       dplyr::mutate(y=1),aes(x=x,y=y),color="blue",size=.75,linetype = "dashed")
      }
    }else{pp=geom_blank()}
    pp
  }, 
  plyr::dlply(X, ~Var2))+ 
  mapply(FUN = function(x){
    z=unique(x$Var2)
    if(!is.element(z,c("$q$","$\\sum_{i=1}^k z_i$","r2"="$R^2$"))){
      breaks = quantile(x$value,probs=seq(0,1,length.out=41))
      geom_histogram(data = x,aes(y=..density..), color=NA,fill=rep(c("black","red","black"),c(1,38,1)),alpha=.5, position="identity",
                     breaks = breaks)}else{
      if(z=="$\\sum_{i=1}^k z_i$"){
        breaks=unique(c(x$value-1/2,x$value+1/2))}
    if(z=="$R^2$"){
          breaks=center(r2_q_grid$r2)
          breaks[breaks>1]<-1
          breaks[breaks<0]<-0
          breaks<-c(breaks[1],breaks[-1][breaks[-length(breaks)]<=max(x$value)])
          breaks<-c(breaks[-length(breaks)][breaks[-1]>=min(x$value)],breaks[length(breaks)])
          breaks=unique(breaks)|>sort()}
    if(z=="$q$"){
      breaks=center(r2_q_grid$q)
      breaks[breaks>1]<-1
      breaks[breaks<0]<-0
      breaks<-c(breaks[1],breaks[-1][breaks[-length(breaks)]<=max(x$value)])
      breaks<-c(breaks[-length(breaks)][breaks[-1]>=min(x$value)],breaks[length(breaks)])
      breaks=unique(breaks)|>sort()}
    breaks2=breaks[-length(breaks)]|>sort()
      colors=c("black","red")[1+(breaks2>=quantile(x$value,probs=0.025)
                                         &breaks2<=quantile(x$value,probs=0.975))]
      geom_histogram(data = x,aes(y=..density..), fill=colors,alpha=.5, color=NA,position="identity",
                     breaks = breaks)    
    }
  }, 
         plyr::dlply(X, ~Var2))+ 
  mapply(FUN = function(x){
    vv=unique(x$Var2)
    z=quantile(x$value,c(.025,.975))
    conf=signif(mean((x$value>=z[1])&(x$value<=z[2])),4)
    confint=paste0("$",conf,"$ C.I.: $[",signif(z[1],4),",",signif(z[2],4),"]$")
    geom_text(data=data.frame(Var2=vv), 
              aes(label = confint, x = -Inf, y = Inf),
              hjust = 0, vjust = 3,  size=geom.text.size/.pt)
    },
  plyr::dlply(X, ~Var2))+theme_bw()+
  facet_wrap(facets = vars(Var2),nrow=1,scales="free")+xlab("")+ylab("")+ 
  theme(axis.text = element_text(size = theme.size, colour="black")) ->the_plot

if(is.element("SweaveLst",installed.packages())){
  SweaveLst::graph2pdffile(texte="print(the_plot)",output = output,
                           widthe=8,heighte=2)
}else{
  ggsave(plot = the_plot,
         filename =output,
         height=9,
         width=15,
         units="cm")  }
output


}


assign_3_convergence_diagnostics_f<-function(merged_samples,output){
  mcmc_object<-list(
    parameters.to.save=dimnames(merged_samples$allsamples)[[2]],
    n.iter=dim(merged_samples$allsamples)[1],
      BUGSoutput=list(n.chains=dim(merged_samples$allsamples)[3],
                      n.thin=100,
                      n.burnin=1000,
                      n.iter=101000,
                      sims.array=merged_samples$allsamples))
    
    
    (merged_samples$allsamples[,c(1:5,66),])|>plyr::alply(3,as.mcmc)|>as.mcmc.list()->X
    X|>gelman.diag()
    X|>gelman.plot()
    X[[1]]|>geweke.diag()
    X[[1]]|>geweke.plot()
}


#'@examples
#'tar_load(assignment_3_data)
#'tar_load(data_dictionnary)
assign_3_rank_f<-function(merged_samples,output,assignment_3_data,data_dictionnary){
  merged_samples$allsamples[,grep("beta",dimnames(merged_samples$allsamples)[[2]]),]|>
    plyr::aaply(c(1,3),rank)->rr
  rr|>plyr::aaply(3,mean)|>order(decreasing = TRUE)->oo
  names(assignment_3_data)[2+oo[c(1:2,length(oo))]]->.vars
  variables<-data_dictionnary|>
    (`rownames<-`)(data_dictionnary$nom)|>
    (`[`)(.vars,)
  set.seed(1)
  priorodds<-rbeta(10000000,10,100)|>(function(x){(1-(1-x)^2)/(1+(1-x)^2)})()|>mean()
  (1-.81)/(1+.81)
  post0=mean(abs(rr[,,oo[1]])>abs(rr[,,oo[2]]))
  BF01_1=((post0)/(1-post0))/priorodds
  post0=mean(abs(rr[,,oo[1]])>abs(rr[,,oo[length(oo)]]))
  BF01_2=((post0)/(1-post0))/priorodds
  
list(variables=variables,BF01_1=BF01_1,BF01_2=BF01_2)  
}
