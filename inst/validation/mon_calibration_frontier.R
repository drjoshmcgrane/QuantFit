# Compare TWO monotone-cone projections for the MON null, on IDENTICAL data:
#   pava : weighted isotonic -> POOLS violating classes to equal values, which
#          destroys class separation and inflates the null (suspected cause of
#          the IIO power loss at J=12).
#   sort : sort each item's class-probabilities along the ordered classes.
#          Enforces monotonicity while preserving each item's spread EXACTLY,
#          so class separation in the null matches the observed data.
# ...vs the incumbent eps rule.
suppressMessages(library(QuantFit)); env<-asNamespace("QuantFit")
rmt<-get("refit_model_type",env); co<-get(".class_orderings",env); monh<-get(".manifest_mon_holds",env)
down_of<-function(P,o) sum(pmax(0,-t(apply(P[,o,drop=FALSE],1,diff))))
stat_of<-function(P){o<-co(ncol(P)); min(apply(o,1,function(x) down_of(P,x)))/nrow(P)}
best_ord<-function(P){o<-co(ncol(P)); o[which.min(apply(o,1,function(x) down_of(P,x))),]}
.pava<-function(y,w){n<-length(y); if(n<=1L)return(y); val<-y; wt<-w; idx<-as.list(seq_len(n)); i<-1L
  while(i<length(val)){ if(val[i]>val[i+1L]+1e-12){ nw<-wt[i]+wt[i+1L]; nv<-(val[i]*wt[i]+val[i+1L]*wt[i+1L])/nw
    val[i]<-nv; wt[i]<-nw; idx[[i]]<-c(idx[[i]],idx[[i+1L]]); val<-val[-(i+1L)]; wt<-wt[-(i+1L)]; idx<-idx[-(i+1L)]
    if(i>1L)i<-i-1L } else i<-i+1L }
  out<-numeric(n); for(b in seq_along(val)) out[idx[[b]]]<-val[b]; out }
mon_null<-function(data,how,C=3L,B=49L,ns=5L,alpha=0.05,seed=NULL){
  data<-as.matrix(data); n<-nrow(data); J<-ncol(data)
  un<-rmt("UN",data,C,ns,TRUE); P<-un$item_probs; pi_c<-un$class_probs
  obs<-stat_of(P); o<-best_ord(P); P0<-P
  P0[,o]<-if(how=="pava") t(apply(P[,o,drop=FALSE],1,function(x) .pava(x,pi_c[o])))
          else            t(apply(P[,o,drop=FALSE],1,sort))
  P0<-pmin(pmax(P0,1e-6),1-1e-6); tP0<-t(P0)
  if(!is.null(seed)) set.seed(seed)
  nulls<-vapply(seq_len(B),function(b){ cls<-sample.int(C,n,replace=TRUE,prob=pi_c)
    sim<-matrix(stats::rbinom(n*J,1,tP0[cls,,drop=FALSE]),n,J); storage.mode(sim)<-"integer"
    f<-tryCatch(rmt("UN",sim,C,ns,TRUE),error=function(e)NULL)
    if(is.null(f)) NA_real_ else stat_of(f$item_probs) },numeric(1))
  nulls<-nulls[!is.na(nulls)]; p<-(1+sum(nulls>=obs))/(length(nulls)+1)
  list(p=p,holds=p>alpha) }
gen<-function(m,J,N,seed,sep=NULL){
  if(!is.null(sep)){ set.seed(seed); C<-3L; beta<-sort(runif(J,-1.5,1.5)); th<-(seq_len(C)-2)*sep
    P<-plogis(outer(th,beta,"-")); cl<-sample.int(C,N,replace=TRUE); d<-matrix(rbinom(N*J,1,P[cl,]),N,J)
  } else if(m=="RM"){ set.seed(seed); b<-sort(runif(J,-2,2)); th<-rnorm(N)
    d<-matrix(rbinom(N*J,1,plogis(outer(th,b,"-"))),N,J)
  } else { d<-simulate_responses(m,n_persons=N,n_items=J,n_classes=3,seed=seed); d<-if(is.list(d))d$data else d }
  storage.mode(d)<-"integer"; d }
cond<-list(list(lab="DM sep0.6 HOLDS",m="DM",J=12L,sep=0.6), list(lab="DM HOLDS",m="DM",J=8L,sep=NULL),
  list(lab="MON HOLDS",m="MON",J=8L,sep=NULL), list(lab="RM HOLDS J6",m="RM",J=6L,sep=NULL),
  list(lab="LCR HOLDS",m="LCR",J=8L,sep=NULL),
  list(lab="IIO VIOL J8",m="IIO",J=8L,sep=NULL), list(lab="IIO VIOL J12",m="IIO",J=12L,sep=NULL),
  list(lab="UN VIOL",m="UN",J=8L,sep=NULL))
R<-as.integer(Sys.getenv("P2_R","10")); od<-Sys.getenv("P2_OUT","p2_out"); dir.create(od,showWarnings=FALSE)
tk<-expand.grid(ci=seq_along(cond),r=seq_len(R))
run<-function(k){ ci<-tk$ci[k]; r<-tk$r[k]; cc<-cond[[ci]]
  f<-file.path(od,sprintf("q_%02d_%02d.csv",ci,r)); if(file.exists(f))return()
  d<-gen(cc$m,cc$J,1500L,6000+100*ci+r,cc$sep)
  a<-tryCatch(mon_null(d,"pava",seed=r),error=function(e)NULL)
  b<-tryCatch(mon_null(d,"sort",seed=r),error=function(e)NULL)
  e<-tryCatch(monh(d,3L,24L,5L,TRUE,0.01,r+7L),error=function(e)NULL)
  write.csv(data.frame(cond=cc$lab,rep=r,
    pava=if(is.null(a))NA else a$holds, sort=if(is.null(b))NA else b$holds,
    eps=if(is.null(e))NA else e$holds),f,row.names=FALSE) }
invisible(parallel::mclapply(seq_len(nrow(tk)),function(k) tryCatch(run(k),error=function(e)NULL),
  mc.cores=as.integer(Sys.getenv("P2_CORES","6"))))
cat("P2 DONE\n")
