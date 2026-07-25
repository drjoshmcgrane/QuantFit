# Hybrid selector on the TI&D DEVELOPMENT data (tid_data/ archive: real
# TI&D-simulated datasets, N=5000, J in {6,12,24,48}, dichotomous, known
# generating models 0-5 = UN/MON/IIO/DM/LCR/RM). Balanced sample over
# model x nI; J=48 excluded (LCR bridge grain ceiling(49/2)=25 classes).
# Resumable: one CSV per (id). Config matches the validated hybrid.
#
# TID_N  rows per dataset. DEFAULT 5000 = the FULL data. (An earlier run used
#        1500; that subsampling is what produced the only blemish in the TI&D
#        result - RM/nI=6 -> IIO - and it was adopted on a since-disproven
#        belief about job kill limits. Set TID_N=1500 to reproduce it.)
# TID_NI comma-separated item counts, default "6,12,24".
suppressMessages(library(QuantFit))
DD <- "tid_data"
gc_ <- read.csv(file.path(DD,"generatingConditions.csv"))
models <- c("UN","MON","IIO","DM","LCR","RM")
scale_of <- c(UN="nominal",MON="ordinal",IIO="ordinal",DM="ordinal",LCR="quant",RM="quant")
K  <- as.integer(Sys.getenv("TID_K","5"))
NS <- as.integer(Sys.getenv("TID_N","5000"))
# nI set is configurable. NOTE: an earlier run excluded nI=24 and subsampled to
# N=1500 on the belief that longer jobs were killed by a machine/harness ceiling.
# That belief was TESTED AND DISPROVEN (a probe ran 14+ min and an 8-worker job
# 16+ min uninterrupted; the earlier kills coincided with several concurrent jobs
# oversubscribing the cores). Run ONE job at a time and full N=5000 is viable --
# which matters because the N=1500 subsampling is what produced the only blemish
# in the TI&D result (the RM/nI=6 -> IIO artifact).
nIs <- as.integer(strsplit(Sys.getenv("TID_NI","6,12,24"), ",")[[1]])
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
  s<-tryCatch(select_model_hybrid(d,n_classes=3L,B=49L,mc.cores=1L,seed=1,verbose=FALSE)$selected,
              error=function(e)NA_character_)
  write.csv(data.frame(id=cs$id,nI=cs$nI,truth=tru,truth_scale=scale_of[tru],
    selected=s,selected_scale=if(is.na(s))NA else scale_of[s]),f,row.names=FALSE) }
cat("TID hybrid:",nrow(pick),"datasets (N subsampled to",NS,")\n")
invisible(parallel::mclapply(seq_len(nrow(pick)),function(k) tryCatch(run(k),error=function(e)NULL),
  mc.cores=as.integer(Sys.getenv("TID_CORES","8"))))
cat("TID HYBRID DONE\n")
