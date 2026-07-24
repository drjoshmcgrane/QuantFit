# Hybrid selector on the TI&D DEVELOPMENT data (tid_data/ archive: real
# TI&D-simulated datasets, N=5000, J in {6,12,24,48}, dichotomous, known
# generating models 0-5 = UN/MON/IIO/DM/LCR/RM). Balanced sample over
# model x nI; J=48 excluded (LCR bridge grain ceiling(49/2)=25 classes,
# infeasible). Each dataset SUBSAMPLED to N=1500 rows (seeded) because the LR
# edge over full N=5000 exceeds this machine's ~20-min background-job ceiling.
# Resumable: one CSV per (id). Config matches the validated hybrid.
suppressMessages(library(QuantFit))
DD <- "tid_data"
gc_ <- read.csv(file.path(DD,"generatingConditions.csv"))
models <- c("UN","MON","IIO","DM","LCR","RM")
scale_of <- c(UN="nominal",MON="ordinal",IIO="ordinal",DM="ordinal",LCR="quant",RM="quant")
K  <- as.integer(Sys.getenv("TID_K","5"))
NS <- as.integer(Sys.getenv("TID_N","1500"))
nIs <- c(6L,12L)   # 24 excluded: LR edge >20min ceiling
# balanced pick: first K ids per (model,nI)
pick <- do.call(rbind, lapply(0:5, function(m) do.call(rbind, lapply(nIs, function(n){
  ids <- gc_$id[gc_$model==m & gc_$nI==n]; head(data.frame(id=ids, model=m, nI=n), K) }))))
out <- Sys.getenv("TIDH_OUT","tid_hybrid_out"); dir.create(out,showWarnings=FALSE)
load_d <- function(id){ e<-new.env(); load(file.path(DD,sprintf("TA%d.Rdata",id)),envir=e)
  d<-get(ls(e)[1],e)$obsData; storage.mode(d)<-"integer"; d }
run <- function(k){ cs<-pick[k,]; f<-file.path(out,sprintf("t_%d.csv",cs$id)); if(file.exists(f))return()
  d<-load_d(cs$id)
  set.seed(cs$id); if(nrow(d)>NS) d<-d[sample.int(nrow(d),NS),,drop=FALSE]
  tru<-models[cs$model+1]
  s<-tryCatch(select_model_manifest(d,n_classes=3L,B=49L,dm_quant="lr",mc.cores=1L,seed=1,verbose=FALSE)$selected,
              error=function(e)NA_character_)
  write.csv(data.frame(id=cs$id,nI=cs$nI,truth=tru,truth_scale=scale_of[tru],
    selected=s,selected_scale=if(is.na(s))NA else scale_of[s]),f,row.names=FALSE) }
cat("TID hybrid:",nrow(pick),"datasets (N subsampled to",NS,")\n")
invisible(parallel::mclapply(seq_len(nrow(pick)),function(k) tryCatch(run(k),error=function(e)NULL),
  mc.cores=as.integer(Sys.getenv("TID_CORES","8"))))
cat("TID HYBRID DONE\n")
