
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
  
  for (do.scale in c(TRUE, FALSE)) {
    
    suffix <- switch(do.scale + 1, "", "scaled_")  ## false, true
    
    D <- 
      array(
        NA, 
        dim = c(length(subjs), length(subjs), length(parcellation$key), 2), 
        dimnames = list(subj_run1 = subjs, subj_run2 = subjs, parcellation$key, type = c("univariate", "multivariate"))
      )
    
    for (parcel.i in seq_along(parcellation$key)) {
      # parcel.i = 1
      
      is.parcel <- schaefer10k == parcel.i
      B <- betas.i[, is.parcel, , ]
      
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
      
      ## multivariate:
      
      D[, , parcel.i, "multivariate"] <- pdist2(B1, B2) / ncol(B1)  ## divide by number of features (vertices)
      
      
      ##  univariate:
      
      B1_bar <- cbind(rowMeans(B1))  ## get means
      B2_bar <- cbind(rowMeans(B2))
      
      D[, , parcel.i, "univariate"] <- pdist2(B1_bar, B2_bar)  / ncol(B1_bar) ## b/c univariate, denominator == 1
      
      
      ## save 
      
      
      p <- arrangeGrob(
        matplot(D[, , parcel.i, "multivariate"]) + labs(title = "multivariate"),
        matplot(D[, , parcel.i, "univariate"]) + labs(title = "univariate"),
        top = parcellation$key[parcel.i],
        nrow = 1
      )
      
      ggsave(
        file.path(fig.dir, paste0("euclidean_", suffix, parcellation$key[parcel.i], ".pdf")), 
        p,
        width = 14, height = 8, units = "cm",
        device = "pdf"
        )
      
    }
    
    
    saveRDS(D, file.path(out.dir, paste0("euclidean_",  suffix, name.task.i, "_", name.glm.i, ".RDS")))
    
    
  }
  
  
  NULL
  
  
}
stopCluster(cl)
(time.run <- Sys.time() - time.start)




