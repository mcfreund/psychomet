
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
suppressWarnings(source(here("code", "_atlases.R")))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))


subjs <- subjs[!subjs %in% "432332"]


glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name.glm = c(
    "baseline_Cues_EVENTS_censored_shifted",
    "baseline_CongruencySwitch_EVENTS_censored_shifted",
    "baseline_ListLength_EVENTS_censored_shifted",
    "baseline_Congruency_EVENTS_censored_shifted"
  ),
  stringsAsFactors = FALSE
)
glminfo <- as.data.table(glminfo)
lev1 = list(
  Axcpt = "BX",
  Cuedts = c("InConInc", "InConNoInc"),
  Stern = "LL5RN",
  Stroop = c("biasInCon", "PC50InCon")
)
lev2 = list(
  Axcpt = "BY",
  Cuedts = c("ConInc", "ConNoInc"),
  Stern = "LL5NN",
  Stroop = c("biasCon", "PC50Con")
)


## loop ----

cl <- makeCluster(nrow(glminfo))
registerDoParallel(cl)
time.start <- Sys.time()
z <- foreach(
  glm.i = seq_len(nrow(glminfo)), .inorder = FALSE, .verbose = TRUE,
  .combine = c,
  .packages = c("mikeutils", "here", "data.table", "ggplot2", "grid", "gridExtra", "dplyr", "abind")
  ) %dopar% {
  # glm.i = 1
  
  res <- enlist(parcellation$key)
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  
  out.dir <- here("out", "fprint", "unimulti_hilo_target_schaefer400-07")
  if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
  
  fig.dir <-file.path(out.dir, "dmat_figs", paste0(name.task.i, "_", name.glm.i))
  if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)
  
  ## read betas:
  
  betas.i <- readRDS(
    here::here("out", "glms", paste0("betas_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS"))
  )
  
  betas.i <- betas.i[, , , !dimnames(betas.i)$subj %in% "432332", ]  ## remove subj with missing data
  
  ## average across target TRs:
  betas.i <- abind(
    apply(betas.i[, lev1[[glm.i]], target.trs[[glm.i]], , ], c("vertex", "subj", "run"), mean),
    apply(betas.i[, lev2[[glm.i]], target.trs[[glm.i]], , ], c("vertex", "subj", "run"), mean),
    along = 0
  )   ## condition, vertex, subj, run
  names(dimnames(betas.i)) <- c("condition", "vertex", "subj", "run")
  dimnames(betas.i)$condition <- c("lev1", "lev2")
  
  for (do.scale in c(FALSE)) {
    
    D <- 
      array(
        NA, 
        dim = c(length(subjs), length(subjs), length(parcellation$key), length(subjs), 2), 
        dimnames = list(subj_run1 = subjs, subj_run2 = subjs, parcellation$key, ndim = NULL, type = c("noncv", "cv"))
      )
    
    for (parcel.i in seq_along(parcellation$key)) {
      # parcel.i = 1
      
      is.parcel <- schaefer10k == parcel.i
      B <- betas.i[, is.parcel, , ]
      
      ## remove vertices that are unresponsive in any subject:
      
      has.signal.allsubjs <- apply(abs(B) <= .Machine$double.eps, "vertex", sum) < 1
      B <- B[, has.signal.allsubjs, , ]
      
      
      if (do.scale) {
        
        nvert <- dim(B)[2]
        B <- sweep(
          B,
          c(1, 3, 4),
          sqrt(apply(B, c("condition", "subj", "run"), function(x) sum(x^2)) / nvert),
          "/"
          )  ## divisive normalization by root mean square
        
      }
      
      B_contrast <- B["lev1", , , ] - B["lev2", , , ]  ## get contrast
      
      B1 <- t(B_contrast[, , 1])  ## separate by runs
      B2 <- t(B_contrast[, , 2])
      
      A1 <- prcomp(B1)$x  ## projections of data onto PCs ('scores')
      A2 <- prcomp(B2)$x
      
      ndims <- min(dim(B)[2], length(subjs))
      
      for (ndim.i in seq_len(ndims)) {
        # ndim.i = 6
        
        A1_n <- cbind(A1[, seq_len(ndim.i)], 1)  ## get reduced-dimension subspace (and append intercept)
        A2_n <- cbind(A2[, seq_len(ndim.i)], 1)
        
        B1_hat <- A1_n %*% coef(.lm.fit(x = A1_n, y = B1))  ## least-squares projection
        B2_hat <- A2_n %*% coef(.lm.fit(x = A2_n, y = B2))
        
        B1_hat_cv <- A2_n %*% coef(.lm.fit(x = A2_n, y = B1))  ## cross-validate: proejct into subspace of other run
        B2_hat_cv <- A1_n %*% coef(.lm.fit(x = A1_n, y = B2))
        
        
        D[, , parcel.i, ndim.i, "noncv"] <- pdist2(B1_hat, B2_hat) / ncol(B1)  ## divide by number of features (verts)
        D[, , parcel.i, ndim.i, "cv"] <- pdist2(B1_hat_cv, B2_hat_cv) / ncol(B1)
        
      }
      
    }
    
    ## save
    
    suffix <- switch(do.scale + 1, "", "scaled_")  ## false, true
    
    saveRDS(D, file.path(out.dir, paste0("euclidean_svd_",  suffix, name.task.i, "_", name.glm.i, ".RDS")))
    
    
  }
  
  
  NULL
  
  
}
stopCluster(cl)
(time.run <- Sys.time() - time.start)




