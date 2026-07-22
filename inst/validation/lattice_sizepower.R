# Focused size/power confirmation for the lattice selector (LC only, fast).
# 6 models x J{6,12,24} x {sigma 0.5,1,2 for RM; 1 else} x 8 reps.
# Question: recovery holds across seeds (not just the 8 hand-picked cases)?
suppressMessages(library(QuantFit))
cells <- rbind(
  expand.grid(model="RM", J=c(6L,12L,24L), sigma=c(0.5,1,2), rep=1:8, stringsAsFactors=FALSE),
  expand.grid(model=c("UN","MON","IIO","DM","LCR"), J=c(6L,12L,24L), sigma=1, rep=1:8, stringsAsFactors=FALSE))
cells$id <- seq_len(nrow(cells)); set.seed(918273); cells$seed <- sample.int(.Machine$integer.max, nrow(cells))
out <- "lattice_sp_tid032"; dir.create(out, showWarnings=FALSE)
gen <- function(g){ set.seed(g$seed)
  if(g$model=="RM"){b<-runif(g$J,-2,2);th<-rnorm(1500,0,g$sigma)
    d<-matrix(rbinom(1500*g$J,1,plogis(outer(th,b,"-"))),1500,g$J)}
  else {d<-simulate_responses(g$model,n_persons=1500,n_items=g$J,n_classes=3,seed=g$seed)
    d<-if(is.list(d))d$data else d}; storage.mode(d)<-"integer"; d }
run1 <- function(i){ g<-cells[cells$id==i,]; f<-file.path(out,sprintf("S%03d.csv",i))
  if(file.exists(f)) return(invisible())
  d<-gen(g)
  lc<-tryCatch(suppressWarnings(select_model_ll(d,n_classes=1:5,B=99,n_starts=5,
      boot_n_starts=5,method="lattice",severity=FALSE,
      seed=g$id*1000003L,mc.cores=1)),error=function(e)NULL)
  row<-data.frame(id=i,model=g$model,J=g$J,sigma=g$sigma,rep=g$rep,
    selected=if(is.null(lc))NA else lc$selected,err=if(is.null(lc))"lc" else "")
  write.csv(row,f,row.names=FALSE) }
cat("lattice size/power:",nrow(cells),"datasets\n")
invisible(parallel::mclapply(cells$id,function(i)
  tryCatch(run1(i),error=function(e)NULL),mc.cores=6))
cat("SIZEPOWER DONE\n")
