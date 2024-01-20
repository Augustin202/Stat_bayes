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
#'tar_load_everything()
#'dataset=assignment_3_data
final_project_run_mcmc<-function(dataset,
                                 r2_q_grid,
                                 nrep,
                                 burning,
                                 nchain){

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
             tt = tt,k = k,phi = 0,r2_q_grid = r2_q_grid,a=1,b=1,aa=1,
             bb=1,nrep=nrep,burning=burning,posterior="beta")}, 
  mc.cores = numberOfCores)
}


final_project_analyse_mcmc_sample<-function(){}

#'@examples
#'
#'
#'tar_load_everything()
#'mcmc_sample=assignment_3_mcmc_sample
final_project_analyse_mcmc_sample<-function(mcmc_sample){
  assignment_3_mcmc_sample|>
    c(list(along=3))|>
    do.call(what=abind::abind)->XX

  XX|>reshape2::melt()->YY
    
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
