## about ----
## 
## sets up environment for qc_group.rmd and renders report.
## 
## mike freund, 2020-01-20

## TODO:
## CHECK ALL PATHS FIRST!

## setup ----

library(here)
library(magrittr)
library(gifti)
library(cifti)
library(abind)
library(data.table)
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

whitening <- function(E, shrinkage = "ledoitwolf") {
  
  S <- cov(E)  ## sample cov
  # corrplot::corrplot(cov2cor(S), method = "color")
  H <- diag(nrow(S))  ## target to shrink towards (Ledoit-Wolf's 'F')
  
  if (shrink %in% c("lw", "ledoitwolf", "LW"))  {
    k <- tawny::shrinkage.intensity(E, H, S)
    lambda <- max(c(0, min(k / nrow(E), 1)))  ## shrinkage factor
  }
  
  S_hat <- lambda * H + (1 - lambda) * S  ## shrunken matrix
  
  solve(S_hat)  ## mahalanobis whitening matrix^2
  
}

# collate_surface_params <- function(name, space = "hcp", xlabels = NULL, pattern = NULL, warn = TRUE)
# {
#   if (space == "hcp") {
#     n.vertices <- 32492
#   }
#   else if (space == "fsave") {
#     n.vertices <- 10242
#   }
#   else {
#     stop("!space %in% c('hcp', 'fsave')")
#   }
#   xlabels.actual <- afni("3dinfo", paste0("-label ", name))
#   xlabels.actual <- unlist(strsplit(xlabels.actual, "\\|"))
#   if (!is.null(xlabels) && !is.null(pattern)) {
#     stop("xlabels and pattern cannot both be specified")
#   }
#   else if (is.null(xlabels) && is.null(pattern)) {
#     xinds <- seq_len(length(xlabels.actual))
#   }
#   else if (!is.null(pattern)) {
#     xinds <- grep(pattern, xlabels.actual)
#     if (length(xinds) < 1) {
#       if (warn) return(NA) else stop("pattern has no matches in sub-brick labels")
#     }
#   }
#   else {
#     xinds <- which(xlabels %in% xlabels.actual)
#     if (length(xinds) < 1)
#       if (warn) return(NA) else stop("no labels match")
#   }
#   gii <- gifti::read_gifti(name)
#   d <- gii$data[xinds]
#   m <- matrix(unlist(d, use.names = FALSE), nrow = length(d),
#               byrow = TRUE)
#   n.params <- length(xinds)
#   if (!all(dim(m) == c(n.params, n.vertices))) {
#     if (warn) return(NA) else stop("dims not expected!")
#   }
#   dimnames(m) <- list(param = xlabels.actual[xinds], vertex = NULL)
#   m
# }


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
# sessi <- sessi[1]
# sessi.short <- sessi.short[1]
n.knots <- 8

## create "label table" (lt) for loop indices

lt <- rbind(
  data.table(
    task     = "Axcpt", 
    variable = c("AY", "BX", "BY", "AX", "Bng", "Ang", "error", "button1", "button2"),
    glm.name = c(rep("Cues_EVENTS_censored", 7), rep("Buttons_censored", 2))
    # contrast = c("Acue_Bcue", "HI_LO_conf", "Nogo_Go", "error_correct", "B1_B2"),
  ),
  data.table(
    task     = "Cuedts", 
    variable = c(
      "InConInc", "InconNoInc", "ConInc", "ConNoInc", "error",
      "SwitchInc", "SwitchNoInc", "RepeatInc", "RepeatNoInc",
      "button1", "button2"
      ),
    glm.name = c(
      rep("CongruencyIncentive_EVENTS_censored", 5), 
      rep("SwitchIncentive_EVENTS_censored", 4),
      rep("Buttons_censored", 2)
      )
    # contrast = c("Inc_NoInc", "InCon_Con", "Switch_Repeat", "error_correct", "B1_B2"),
  ),
  data.table(
    task     = "Stern", 
    variable = c("LL5NP", "LL5NN", "LL5RN", "not5NP", "not5NN", "not5RN", "error", "button1", "button2"),
    glm.name = c(rep( "ListLength_EVENTS_censored", 7), rep("Buttons_censored", 2))
    # contrast = c("RN_NN_all", "RN_NN_LL5", "RN_NN_not5", "not5_LL5", "error_correct", "B1_B2"),
  ),
  data.table(
    task     = "Stroop", 
    variable = c("PC50Con", "PC50InCon", "biasCon", "biasInCon", "error", "button1", "button2"),
    glm.name = c(rep( "ListLength_EVENTS_censored", 5), rep("Buttons_censored", 2))
    # contrast = c("InCon_Con_bias", "InCon_Con_PC50", "InCon_Con_PC50bias", "error_correct"),
  )
)

contrast <- list(
  Axcpt  = c("Acue_Bcue", "HI_LO_conf", "Nogo_Go", "error_correct", "B1_B2"),
  Cuedts = c("Inc_NoInc", "InCon_Con", "Switch_Repeat", "error_correct", "B1_B2"),
  Stern  = c("RN_NN_all", "RN_NN_LL5", "RN_NN_not5", "not5_LL5", "error_correct", "B1_B2"),
  Stroop = c("InCon_Con_bias", "InCon_Con_PC50", "InCon_Con_PC50bias", "error_correct")
)

statistics.wrun <- c(combo.paste(c("cor", "euc"), c("raw", "unn", "mnn")), "mean")
statistics.crun <- c(combo.paste(c("cor", "euc", "cvcor", "cveuc"), c("raw", "unn", "mnn")))

## lists of arrays

d.within <- vector("list", length(tasks)) %>% setNames(tasks)  ## d within
d.across <- d.within  ## d between

for (task.i in seq_along(tasks)) {

  d.within[[task.i]] <- array(
    NA,
    dim = c(
      subj  = length(subjs), 
      sess  = length(sessi),
      roi   = length(parcellation$key), 
      cont  = length(contrast[[task.i]]),
      stat  = length(statistics.wrun),
      knot  = n.knots,
      run   = 2
    ),
    dimnames = list(
      subj  = subjs, 
      sess  = sessi,
      roi   = parcellation$key, 
      cont  = contrast[[task.i]],
      stat  = statistics.wrun,
      knot  = paste0("knot", seq_len(n.knots)),
      run   = c("run1", "run2")
    )
  )
  
  d.across[[task.i]] <- array(
    NA,
    dim = c(
      subj  = length(subjs), 
      sess  = length(sessi),
      roi   = length(parcellation$key), 
      cont  = length(contrast[[task.i]]),
      stat  = length(statistics.crun),
      knot  = n.knots,
      run   = 2
    ),
    dimnames = list(
      subj  = subjs, 
      sess  = sessi,
      roi   = parcellation$key, 
      cont  = contrast[[task.i]],
      stat  = statistics.crun,
      knot  = paste0("knot", seq_len(n.knots)),
      run   = c("run1", "run2")
    )
  )
  
}

r <- d.within  ## lis r for split-half (across runs) reliabilities


## loop ----


for (task.i in seq_along(tasks)) {
  # task.i = 1
  
  name.task.i <- tasks[task.i]
  
  for (sess.i in seq_along(sessi)) {
    # sess.i = 1
    
    name.sess.i <- sessi[sess.i]
    
    for (subj.i in seq_along(subjs)) {
      ## subj.i = 1
      
      name.subj.i <- subjs[subj.i]
      
      lt.i <- lt[task == name.task.i, c("variable", "glm.name")]
      
      for (name.glm.i in unique(lt.i$glm.name)) {
        # name.glm.i = unique(lt.i$glm.name)[1]
      
        vars.i <- lt.i[glm.name == name.glm.i]$variable
          
        dirs <- combopaste(
          file.path(dir.results, name.subj.i, name.task.i),
          paste0("/", name.sess.i, "/", name.sess.i, "_", name.glm.i)
        )

        betas <- vector("list", 4) %>% setNames(combo_paste(c("run1", "run2"), c("L", "R")))
        resid <- betas
        
        for (name.run.i in c("run1", "run2")) {
          # name.run.i = "run1"
          for (name.hemi.i in c("L", "R")) {
            # name.hemi.i = "L"
            
            f.betas <- file.path(paste0(dirs, "_", name.run.i), paste0("stats_", subjs[subj.i], "_", name.hemi.i, ".func.gii"))
            f.resid <- file.path(paste0(dirs, "_", name.run.i), paste0("wherr_", subjs[subj.i], "_", name.hemi.i, ".func.gii"))
            
            file.is.missing <- !file.exists(f.betas) | !file.exists(f.resid)
            if (file.is.missing) next  ## go to next task!
            
            name.run.hemi.i <- paste0(name.run.i, "_", name.hemi.i)
            
            betas[[name.run.hemi.i]] <- collate_surface_params(
              f.betas, 
              pattern = combopaste(vars.i, c("#[0-9]_Coef")) %>% paste0(collapse = "|")
              )
            
            resid[[name.run.hemi.i]] <- read_gifti2matrix(f.resid)
            
          }  ## end hemi loop
        }  ## end run loop
        
        betas <- abind(
          run1 = cbind(betas[["run1_L"]], betas[["run1_R"]]),
          run2 = cbind(betas[["run2_L"]], betas[["run2_L"]]),
          rev.along = 0
        )
        resid <- abind(
          run1 = cbind(resid[["run1_L"]], resid[["run1_R"]]),
          run2 = cbind(resid[["run2_L"]], resid[["run2_L"]]),
          rev.along = 0
        )
        
        
        ## estimation ----
        
        for (roi.i in seq_along(parcellation$key)) {
          # roi.i = 1
          name.roi.i <- parcellation$key[roi.i]
          betas.i <- betas[, parcellation$atlas == roi.i, ]
          
          for (knot.i in seq_len(n.knots)) {
            # knot.i = 1
            
            ## mask (for knot and region)
            
            rows.knot.i <- grep(paste0("#", knot.i - 1, "_Coef"), rownames(betas.i))  ## minus one b/c 0-based ind.
            betas.ii <- betas.i[rows.knot.i, , ]
            
            for (contrast.i in contrast[[task.i]]) {
              
              ## group effects
              
              ## check order of rownames
              
              ## get contrasts
              contrast
              m <- rbind(
                ##    ax  ay  ang  bx  by  bng  error
                cbind(1/2, 1/2, 0, -1/2, -1/2, 0, 0),
                cbind(-1/2, 1/2, 0, 1/2, -1/2, 0, 0),
                cbind(-1/4, -1/4, 1/2, -1/4, -1/4, 1/2, 0),
                cbind(-1/4, -1/4, 0, -1/4, -1/4, 0, 1/2)
              )
              dimnames(m) <- list(contrast = contrast[[task.i]][-5], param = c("AX", "AY", "Ang", "BX", "BY", "Bng", "error"))
              
              contrasts.i <- abind(
                run1 = m %*% betas.ii[, , "run1"],
                run2 = m %*% betas.ii[, , "run2"],
                rev.along = 0
              )
              
              
              
              ## within-run
              
              diag(cor(t(contrasts.i[, , 1]), t(contrasts.i[, , 2])))
              dist(contrasts.i[, , 1], contrasts.i[, , 2])
              
              ## across-run
              
              cor(t(betas.ii[, , 1]), t(betas.ii[, , 2]))  ## cross-run (correlation, raw)
              
              cor(t(betas.ii[, , 1]), t(betas.ii[, , 2]))  ## cross-run (correlation, cross-validated)
              
              cor(t(betas.ii[, , 1]), t(betas.ii[, , 2]))  ## cross-run (euclidean, raw)
              
              
              
              betas.ii
              
              
              
              betas.i[]
              
              names(dimnames(d.within$Axcpt))
              
              
              d.within[[name.task.i]][subj.i, sess.i, roi.i, , , knot.i, ]  ## subj, sess, roi, cont, stat, knot, run
              
              
              
              dimnames(d.within$Axcpt)
              dimnames(d.across$Axcpt)
              

            }
            
            
            
            
          }
          

          
        }
        
        
        ## reliabilities
        
        ## TODO: get reliabilities and effect sizes
        ## TODO: store in d.wrun, d.crun, and r
        
      }  ## glm loop end
      
      
      
    }  ## subj loop end
    
    
  }  ## end session loop
  
  
  ## save task array
}







for (subj.i in seq_along(subjs)) {
  # subj.i = 1
  
  name.subj.i <- subjs[subj.i]
  
  for (task.j in seq_along(tasks)) {
    # task.j = 1
    
    name.task.j <- tasks[task.j]
    
    ## build paths to GLM
    
    dirs <- combopaste(
      file.path(dir.results, name.subj.i, name.task.j),
      paste0("/", sessi, "/", sessi, "_", label.glm[name.task.j], "_EVENTS_censored")
    )
    
    f1l <- file.path(paste0(dirs, "_run1"), paste0("stats_", subjs[subj.i], "_L.func.gii"))
    f1r <- file.path(paste0(dirs, "_run1"), paste0("stats_", subjs[subj.i], "_R.func.gii"))
    f2l <- file.path(paste0(dirs, "_run2"), paste0("stats_", subjs[subj.i], "_L.func.gii"))
    f2r <- file.path(paste0(dirs, "_run2"), paste0("stats_", subjs[subj.i], "_R.func.gii"))
    
    f1l.w <- file.path(paste0(dirs, "_run1"), paste0("wherr_", subjs[subj.i], "_L.func.gii"))
    f1r.w <- file.path(paste0(dirs, "_run1"), paste0("wherr_", subjs[subj.i], "_R.func.gii"))
    f2l.w <- file.path(paste0(dirs, "_run2"), paste0("wherr_", subjs[subj.i], "_L.func.gii"))
    f2r.w <- file.path(paste0(dirs, "_run2"), paste0("wherr_", subjs[subj.i], "_R.func.gii"))
    
    ## read ----
    
    file.is.missing <- any(
      !vapply(
        list(f1l, f1r, f2l, f2r, f1l.w, f1r.w, f2l.w, f2r.w),
        file.exists,
        logical(1)
      )
    )
    if (file.is.missing) next  ## go to next task!
    
    ## read contrasts
      
    pat <- label.contrast[task.j] %>% combopaste(., c(".._Coef", ".._Tstat")) %>% paste0(collapse = "|")

    stats <- abind(
      run1 = cbind(collate_surface_params(f1l, pattern = pat), collate_surface_params(f1r, pattern = pat)),
      run2 = cbind(collate_surface_params(f2l, pattern = pat), collate_surface_params(f2r, pattern = pat)),
      rev.along = 0
    )
    dimnames(stats) <- list(stat = c("b", "t"), vertex = NULL, run = c("run1", "run2"))
    
    ## read estimates
    ##TODO (embed in above code)
    ## - pull all necessary params; then rearrange once in memory
    
    ## read residuals
    
    resid <- abind(
      run1 = cbind(read_gifti2matrix(f1l.w), read_gifti2matrix(f1r.w)),
      run2 = cbind(read_gifti2matrix(f2l.w), read_gifti2matrix(f2r.w)),
      rev.along = 0
    )
    names(dimnames(resid)) <- c("tr", "vertex", "run")

    ## ingore / skip errors
    
    failed.to.pull.data <- class(stats) != "array" | class(resid) != "array"
    if (failed.to.pull.data) next  ## skip to next task!
    

    ## loop for ROI statisics ----
    
    ## TODO: parellelize the prewhitening...
    for (roi.i in seq_along(parcellation$key)) {
      # roi.i = 1
      
      ## define mask
      
      mask.i <- parcellation$atlas == roi.i
      
      ## get error timecourses
      
      E <- resid[, mask.i, ]
      
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
      
      betas <- stats["b", mask.i, ]
      tstat <- stats["t", mask.i, ]
      
      ## also remove from betas and tstats!
      betas <- betas[!is.allzero, ]
      tstat <- tstat[!is.allzero, ]
      
      r.b <- cor(betas)[1, 2]
      r.t <- cor(tstat)[1, 2]
      r.w <- cor(W[, , 1] %*% betas[, 1], W[, , 2] %*% betas[, 2])
      
      r[subj.i, roi.i, task.j, ] <- c(r.b, r.t, r.w)  ## order corresponds!
      
      
      ## get effects ----
      
      
      ## univariate
      
      d[subj.i, roi.i, task.j, "mean", ] <- colMeans(betas)
      
      ## pearson
      
      ## pearson w/ univ nn
      
      ## pearson w/ multiv nn
      
      ## euclidean
      
      ## mahanlanobis (prewhitened euclidean's)
      
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
