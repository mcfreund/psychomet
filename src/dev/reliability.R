## about ----
## 
## sets up environment for qc_group.rmd and renders report.
## 
## mike freund, 2020-01-20


## setup ----

library(here)
library(magrittr)
library(gifti)
library(cifti)
library(abind)
library(mikeutils)

## functions

are.equal.dims <- function(l) {
  m <- vapply(l, dim, numeric(3))
  all.equal(m[, 1], m[, 2])
}

are.equal.dimnames <- function(l) {
  l2 <- lapply(l, function(.) unlist(dimnames(.)))
  all.equal(l2[[1]], l2[[2]])
}

whitening <- function(E, shrink = TRUE) {
  
  S <- cov(E)  ## sample cov
  # corrplot::corrplot(cov2cor(S), method = "color")
  H <- diag(nrow(S))  ## target to shrink towards (Ledoit-Wolf's 'F')
  
  if (shrink)  {
    k <- tawny::shrinkage.intensity(E, H, S)
    lambda <- max(c(0, min(k / nrow(E), 1)))  ## shrinkage factor
  } else lambda <- 0
  
  S_hat <- lambda * H + (1 - lambda) * S  ## shrunken matrix
  
  expm::sqrtm(solve(S_hat))  ## mahalanobis whitening matrix
  
}

collate_surface_params <- function (name, space = "hcp", xlabels = NULL, pattern = NULL, warn = TRUE)
{
  if (space == "hcp") {
    n.vertices <- 32492
  }
  else if (space == "fsave") {
    n.vertices <- 10242
  }
  else {
    stop("!space %in% c('hcp', 'fsave')")
  }
  xlabels.actual <- afni("3dinfo", paste0("-label ", name))
  xlabels.actual <- unlist(strsplit(xlabels.actual, "\\|"))
  if (!is.null(xlabels) && !is.null(pattern)) {
    stop("xlabels and pattern cannot both be specified")
  }
  else if (is.null(xlabels) && is.null(pattern)) {
    xinds <- seq_len(length(xlabels.actual))
  }
  else if (!is.null(pattern)) {
    xinds <- grep(pattern, xlabels.actual)
    if (length(xinds) < 1) {
      if (warn) return(NA) else stop("pattern has no matches in sub-brick labels")
    }
  }
  else {
    xinds <- which(xlabels %in% xlabels.actual)
    if (length(xinds) < 1)
      if (warn) return(NA) else stop("no labels match")
  }
  gii <- gifti::read_gifti(name)
  d <- gii$data[xinds]
  m <- matrix(unlist(d, use.names = FALSE), nrow = length(d),
              byrow = TRUE)
  n.params <- length(xinds)
  if (!all(dim(m) == c(n.params, n.vertices))) {
   if (warn) return(NA) else stop("dims not expected!")
  }
  dimnames(m) <- list(param = xlabels.actual[xinds], vertex = NULL)
  m
}

## variables

parcellation <- read_atlas("schaefer400")

dir.ccp.hcp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
dir.subsubj <- "/data/nil-external/ccp/witzermanm/AFNI_ANALYSIS_SUBSUBJECT"
dir.results <- file.path(dir.subsubj, "RESULTS_RUNWISE")

subjs <- list.dirs(dir.results, recursive = FALSE, full.names = FALSE)
# subjs <- subjs[1]
tasks <- c("Axcpt", "Cuedts", "Stern", "Stroop")
# tasks <- tasks[1]
sessi <- c("baseline", "proactive", "reactive")
sessi.short <- c("Bas", "Pro", "Rea")
sessi <- sessi[1]
sessi.short <- sessi.short[1]

# label.params.main <- list(
#   Axcpt  = c("AY", "BX", "BY", "AX", "Bng", "Ang"),
#   Cuedts = c("InConInc", "InconNoInc", "ConInc", "ConNoInc"),
#   Stern  = c("LL5RN", "not5RN", "LL5NN", "not5NN"),
#   Stroop = c("biasInCon", "PC50InCon", "biasCon", "pc50Con")
# )
# label.params.alt <- list(
#   Axcpt  = c("AY", "BX", "BY", "AX", "Bng", "Ang"),
#   Cuedts = c("SwitchNoInc", "SwitchInc", "RepeatNoInc", "RepeatInc"),
#   Stern  = c("LL5RN", "not5RN", "LL5NN", "not5NN"),
#   Stroop = c("PC50InCon", "pc50Con")
# )
label.contrast <- list(
  Axcpt = c("HI_LO_conf", "Nogo_Go"),
  Cuedts = c("InCon_Con", "Switch_Repeat"),
  Stern = c("RN_NN_all", "RN_NN_LL5"),
  Stroop = c("InCon_Con_PC50bias", "InCon_Con_PC50")
)
label.glm <- list(
  Axcpt  = c("Cues", "Cues"),
  Cuedts = c("CongruencyIncentive", "SwitchIncentive"),
  Stern  = c("ListLength", "ListLength"),
  Stroop = c("Congruency", "Congruency")
)


# n.knots <- 0:7
# n.tr <- 610
# n.vertices <- 32492

r.stats <- c("b", "t", "prewhitened")
d.stats <- c("mean", "r", "d", "rprew", "mahal")

r <- array(
  NA,
  dim = c(
    subj  = length(subjs), 
    roi   = length(parcellation$key), 
    task  = length(tasks), 
    contr = 2,  ## set manually; two for each (main, alt)
    stat  = length(r.stats)
    ),
  dimnames = list(
    subj  = subjs, 
    roi   = parcellation$key, 
    task  = tasks, 
    contr = c("main", "alt"),
    stat  = r.stats
    )
)  ## for reliabilities

d <- array(
  NA,
  dim = c(
    subj  = length(subjs), 
    roi   = length(parcellation$key), 
    task  = length(tasks), 
    contr = 2,
    stat  = length(d.stats),
    run   = 2
    ),
  dimnames = list(
    subj  = subjs, 
    roi   = parcellation$key, 
    task  = tasks, 
    contr = c("main", "alt"),
    stat  = d.stats,
    run   = c("run1", "run2")
    )
)  ## for effects


for (subj.i in seq_along(subjs)) {
  # subj.i = 1
  
  name.subj.i <- subjs[subj.i]
  
  for (task.j in seq_along(tasks)) {
    # task.j = 1
    
    name.task.j <- tasks[task.j]
    
    ## build paths to GLM
    
    dirs <- combopaste(
      file.path(dir.results, name.subj.i, name.task.j),
      paste0("/", sessi, "/", sessi, "_", label.glm[[name.task.j]], "_EVENTS_censored")
    )
    
    dirs.run1 <- paste0(dirs, "_run1")
    dirs.run2 <- paste0(dirs, "_run2")
    
    
    ## read ----
    
    stats <- setNames(vector("list", length(label.contrast[[name.task.j]])), label.contrast[[name.task.j]])
    param <- stats
    resid <- stats
    
    for (glm.i in seq_along(stats)) {
      # glm.i = 1
      
      ## read contrasts
      
      pat <- label.contrast[[task.j]][glm.i] %>% combopaste(., c(".._Coef", ".._Tstat")) %>% paste0(collapse = "|")
      
      stats[[glm.i]] <- abind(
        run1 = cbind(
          collate_surface_params(
            file.path(dirs.run1[glm.i], paste0("stats_", subjs[subj.i], "_L.func.gii")),
            pattern = pat
          ),
          collate_surface_params(
            file.path(dirs.run1[glm.i], paste0("stats_", subjs[subj.i], "_R.func.gii")),
            pattern = pat
          )
        ),
        run2 = cbind(
          collate_surface_params(
            file.path(dirs.run2[glm.i], paste0("stats_", subjs[subj.i], "_L.func.gii")),
            pattern = pat
          ),
          collate_surface_params(
            file.path(dirs.run2[glm.i], paste0("stats_", subjs[subj.i], "_R.func.gii")),
            pattern = pat
          )
        ),
        rev.along = 0
      )
      names(dimnames(stats[[glm.i]])) <- c("stat", "vertex", "run")
      
      ## read estimates
      ##TODO (embed in above code)
      ## - pull all necessary params; then rearrange once in memory
      
      ## read residuals
      
      resid[[glm.i]] <- abind(
        run1 = cbind(
          read_gifti2matrix(file.path(dirs.run1[glm.i], paste0("wherr_", subjs[subj.i], "_L.func.gii"))),
          read_gifti2matrix(file.path(dirs.run1[glm.i], paste0("wherr_", subjs[subj.i], "_R.func.gii")))
        ),
        run2 = cbind(
          read_gifti2matrix(file.path(dirs.run2[glm.i], paste0("wherr_", subjs[subj.i], "_L.func.gii"))),
          read_gifti2matrix(file.path(dirs.run2[glm.i], paste0("wherr_", subjs[subj.i], "_R.func.gii")))
        ),
        along = 3
      )
      names(dimnames(resid[[glm.i]])) <- c("tr", "vertex", "run")
      
    }
    
    ## temporarily dump in list for processing
    
    l <- list(stats = stats, resid = resid)
    

    
    ## ingore / skip errors
    
    failed.to.pull.data <- is.na(unlist(lapply(c(l$stats, l$resid), is.na)))
    if (any(failed.to.pull.data)) next  ## skip to next task!
      
    ## check if all good ----      
    # 
    # is.obj.ok <- vapply(l, function(.) are.equal.dims(.), logical(1))
    # 
    # if (any(!is.obj.ok)) stop("obj not ok")
    
    
    ## format ----
    
    ## reshape to 4D array

    la <- lapply(  # la: list of arrays
      l,
      function(.) {
        array(
          unlist(., use.names = FALSE),
          dim = c(dim(.[[1]]), 2),
          dimnames = c(dimnames(.[[1]]), list(contr = c("main", "alt")))
          
        )
      }
    )
    
    dimnames(la$stats)$stat <- c("b", "t")
    
    ## re-assign each element of list to original object name
    
    for (ii in seq_along(la)) assign(names(la)[ii], la[[ii]])
    
    ## clean up env
    
    rm(l, la, ii)
    gc()
    
    
    ## loop for ROI statisics ----
    
    for (roi.i in seq_along(parcellation$key)) {
      # roi.i = 1
      
      ## define mask
      
      mask.i <- parcellation$atlas == roi.i
      
      for (contr.i in c("main", "alt")) {
        # contr.i = "main"
        
        ## get error timecourses
        
        E <- resid[, mask.i, , contr.i]
        
        ## remove bad timepoints and bad vertices:

        n.tr <- dim(E)[1]
        n.vertex <- dim(E)[2]
        E.bin <- E == 0
        is.allzero <- rowSums(apply(E.bin, c("vertex", "run"), sum) == n.tr) > 0
        is.censored <- rowSums(apply(E.bin, c("tr", "run"), sum) == n.vertex) > 0
        E <- E[!is.censored, !is.allzero, ]
        
        ## estimate prewhitening matrix
        
        W <- array(NA, dim = c(dim(E)[2], dim(E)[2], 2), dimnames = list("n" = NULL, "n" = NULL, "run" = NULL))
        for (run.i in 1:2) W[, , run.i] <- whitening(E[, , run.i], shrink = FALSE)
        
        ## pattern reliabilities ----
        
        betas <- stats["b", mask.i, , contr.i]
        tstat <- stats["t", mask.i, , contr.i]
        
        ## also remove from betas and tstats!
        betas <- betas[!is.allzero, ]
        tstat <- tstat[!is.allzero, ]
        
        r.b <- cor(betas)[1, 2]
        r.t <- cor(tstat)[1, 2]
        r.w <- cor(W[, , 1] %*% betas[, 1], W[, , 2] %*% betas[, 2])
        
        r[subj.i, roi.i, task.j, "main", ] <- c(r.b, r.t, r.w)

        
        ## get effects ----
        
        
        ## univariate
        
        d[subj.i, roi.i, task.j, "main", "mean", ] <- colMeans(betas)
        
        ## pearson
        
        ## pearson w/ univ nn
        
        ## pearson w/ multiv nn
        
        ## euclidean
        
        ## mahanlanobis (prewhitened euclidean's)
      
      
      }
      
    }  ## end ROI loop
    
    print(paste0("task ", tasks[task.j]))
    
  }  ## end task loop
  
  print(paste0("****subj ", subjs[subj.i], "****"))
  
}  ## end subject loop


## save ----

saveRDS(r, here("out", "runwise", "qc_group", "reliability_hilo_baseline_schaefer400.rds"))
saveRDS(d, here("out", "runwise", "qc_group", "estimates-win-run_hilo_baseline_schaefer400.rds"))


## render ----

# rmarkdown::render(
#   here("analyses", "runwise", "qc_subj", "qc_subj.rmd")
#   )
# # output_file
# output_dir
# output_format = "html_document"
