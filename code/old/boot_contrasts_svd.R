
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
suppressWarnings(source(here("code", "_atlases.R")))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))


subjs <- subjs[!subjs %in% "432332"]
nresamp <- 1E4
out.dir <- here("out", "fprint", "unimulti_hilo_target_schaefer400-07")
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

lt <- lower.tri(diag(length(subjs)))
ut <- upper.tri(diag(length(subjs)))


idi <- function(d, ii, lt, ut) {
  # ii <- sample.int(54, replace = TRUE); d <- D_i[, , 1]; lt = lower.tri(diag(length(ii))); 
  # ut = upper.tri(diag(length(ii)))
  
  di <- d[ii, ii]
  mean(c(di[lt], di[ut])) - mean(diag(di))
  
}

idi_contrast <- function(d, ii, lt, ut) {
  # ii <- sample.int(54, replace = TRUE); d <- D_i; 
  # lt = lower.tri(diag(length(ii))); ut = upper.tri(diag(length(ii)))
  
  di <- d[ii, ii, ]
  u <- di[, , "univariate"]
  m <- di[, , "multivariate"]
  
  contrast.m <- mean(c(m[lt], m[ut])) - mean(diag(m))
  contrast.u <- mean(c(u[lt], u[ut])) - mean(diag(u))
  
  contrast.m - contrast.u
  
}

freqp <- function(x, alternative = "two.sided") {
  
  greater <- sum(x < 0) / length(x)
  lower <- 1 - greater
  two.sided <- 2*min(greater, lower)
  
  p <- c(greater = greater, lower = lower, two.sided = two.sided)
  
  p[[alternative]]
  
}



## loop ----

time.start <- Sys.time()
l <- enlist(tasks)
for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 1
  
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  
  l.scale <- enlist(c("unscaled"))
  for (scale.i in seq_along(l.scale)) {
    # scale.i = 1; ## 1 == unscaled
    
    suffix <- switch(scale.i, "", "scaled_")  ## SAME ORDER AS names(l.scale)!!!!
    
    D <- readRDS(file.path(out.dir, paste0("euclidean_svd_",  suffix, name.task.i, "_", name.glm.i, ".RDS")))
    
    
    # time.start <- Sys.time()
    cl <- makeCluster(n.cores/2)
    registerDoParallel(cl)
    l.parcel <- foreach(
      parcel.i = seq_along(parcellation$key), 
      .packages = c("boot", "data.table"),
      # .combine = function(x) data.table::rbindlist(x, use.names = TRUE, idcol = "parcel"),
      .final = function(x) setNames(x, parcellation$key)
      ) %dopar% {
    
    # for (parcel.i in seq_along(parcellation$key)) {
      # parcel.i = 1
      
      D_i <- D[, , parcel.i, , ]
      ndims <- which(apply(D_i, 3, function(x) any(!is.na(x))))
      
      istats <- vector("list", length(ndims))
      for (ndim.i in ndims) {
        # ndim.i = 1
        
        res.noncv <- boot(D_i[, , ndim.i, "noncv"], idi, R = nresamp, parallel = "multicore")
        res.cv    <- boot(D_i[, , ndim.i, "cv"], idi, R = nresamp, parallel = "multicore")

        bca.cv    <- boot.ci(res.cv, type = "bca")$bca
        bca.noncv <- boot.ci(res.noncv, type = "bca")$bca
        
        istats[[ndim.i]] <- 
          
          list(
            
            estimate = res.noncv$t0, se = sd(res.noncv$t), 
            p = freqp(res.noncv$t, "greater"), lb = bca.noncv[4], ub = bca.noncv[5],
            
            estimate_cv = res.cv$t0, se_cv = sd(res.cv$t), 
            p_cv = freqp(res.cv$t, "greater"), lb_cv = bca.cv[4], ub_cv = bca.cv[5]
            
          )
        
      }
      
      
      rbindlist(istats, use.names = TRUE, idcol = "ndim")  ## return
      
      
    }
    stopCluster(cl)
    # (time.end <- Sys.time() - time.start)
    
    l.scale[[scale.i]] <- rbindlist(l.parcel, use.names = TRUE, idcol = "parcel")
    
  }
  
  
  l[[glm.i]] <- rbindlist(l.scale, use.names = TRUE, idcol = "transform")
  cat(paste0("----\n", glm.i, " of ", nrow(glminfo), " glms done\n"))
  print(Sys.time() - time.start)
  cat("(since script start)\n")
  
}
(time.end <- Sys.time() - time.start)

## combine and save:

results <- rbindlist(l, use.names = TRUE, idcol = "task")

fwrite(results, file.path(out.dir, paste0("idi_group_svd.csv")))


(time.run <- Sys.time() - time.start)

