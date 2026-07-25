# Validate the poset rung end-to-end on the real TI&D datasets that BOTH
# selectors call nominal. Expect: true-UN -> antichain; misclassified (IIO/MON)
# -> partial.
suppressMessages(library(QuantFit))
h<-do.call(rbind,lapply(list.files("QuantFit/inst/validation/tid_realdata_results","h_",full.names=TRUE),read.csv))
both<-h[h$lattice=="UN" & h$hybrid=="UN", c("id","truth","nI")]
od<-Sys.getenv("PR_OUT"); dir.create(od,showWarnings=FALSE)
run<-function(k){ cs<-both[k,]; f<-file.path(od,sprintf("r_%d.csv",cs$id)); if(file.exists(f))return()
  e<-new.env(); load(sprintf("tid_data/TA%d.Rdata",cs$id),envir=e)
  d<-get(ls(e)[1],e)$obsData; storage.mode(d)<-"integer"
  r<-tryCatch(select_model_hybrid(d,n_classes=3L,B=49L,mc.cores=1L,seed=1),error=function(e)NULL)
  if(is.null(r)||is.null(r$poset)) return()
  write.csv(data.frame(id=cs$id,truth=cs$truth,nI=cs$nI,selected=r$selected,
    shape=r$poset$shape,comparable=r$poset$comparable,lo=r$poset$lo),f,row.names=FALSE) }
invisible(parallel::mclapply(seq_len(nrow(both)),function(k) tryCatch(run(k),error=function(e)NULL),mc.cores=8))
cat("PR DONE\n")
