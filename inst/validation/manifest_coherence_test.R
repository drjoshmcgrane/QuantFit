# Manifest MH axis vs the model-based eps axis, identical data.
# MH HOLDS on: DM, MON, LCR, RM   |  MH VIOLATED on: IIO, UN
suppressMessages(library(QuantFit)); env<-asNamespace("QuantFit")
mh<-get(".manifest_mh_holds",env); epsh<-get(".manifest_mon_holds",env)
gen<-function(m,J,N,seed,sep=NULL){
  if(!is.null(sep)){ set.seed(seed); C<-3L; b<-sort(runif(J,-1.5,1.5)); th<-(seq_len(C)-2)*sep
    P<-plogis(outer(th,b,"-")); cl<-sample.int(C,N,replace=TRUE); d<-matrix(rbinom(N*J,1,P[cl,]),N,J)
  } else if(m=="RM"){ set.seed(seed); b<-sort(runif(J,-2,2)); th<-rnorm(N)
    d<-matrix(rbinom(N*J,1,plogis(outer(th,b,"-"))),N,J)
  } else { d<-simulate_responses(m,n_persons=N,n_items=J,n_classes=3,seed=seed); d<-if(is.list(d))d$data else d }
  storage.mode(d)<-"integer"; d }
cond<-list(list(l="DM HOLDS",m="DM",J=8L,s=NULL), list(l="MON HOLDS",m="MON",J=8L,s=NULL),
 list(l="LCR HOLDS",m="LCR",J=8L,s=NULL), list(l="RM HOLDS J6",m="RM",J=6L,s=NULL),
 list(l="DM sep0.6 HOLDS",m="DM",J=12L,s=0.6),
 list(l="IIO VIOL J8",m="IIO",J=8L,s=NULL), list(l="IIO VIOL J12",m="IIO",J=12L,s=NULL),
 list(l="UN VIOL",m="UN",J=8L,s=NULL))
R<-as.integer(Sys.getenv("MH_R","10")); od<-Sys.getenv("MH_OUT","mh_out"); dir.create(od,showWarnings=FALSE)
tk<-expand.grid(ci=seq_along(cond),r=seq_len(R))
run<-function(k){ ci<-tk$ci[k]; r<-tk$r[k]; cc<-cond[[ci]]
  f<-file.path(od,sprintf("m_%02d_%02d.csv",ci,r)); if(file.exists(f))return()
  d<-gen(cc$m,cc$J,1500L,7000+100*ci+r,cc$s)
  a<-tryCatch(mh(d,3L,49L,5L,TRUE,0.05,r),error=function(e)NULL)
  b<-tryCatch(epsh(d,3L,24L,5L,TRUE,0.01,r+7L),error=function(e)NULL)
  write.csv(data.frame(cond=cc$l,rep=r,
    mh=if(is.null(a))NA else a$holds, mh_p=if(is.null(a))NA else round(a$p,3),
    eps=if(is.null(b))NA else b$holds),f,row.names=FALSE) }
invisible(parallel::mclapply(seq_len(nrow(tk)),function(k) tryCatch(run(k),error=function(e)NULL),
  mc.cores=as.integer(Sys.getenv("MH_CORES","8"))))
cat("MH TEST DONE\n")
