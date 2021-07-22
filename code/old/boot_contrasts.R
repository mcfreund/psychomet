
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
  
  
  l.scale <- enlist(c("unscaled", "scaled"))
  for (scale.i in 1:2) {
    # scale.i = 1; ## 1 == unscaled
    
    suffix <- switch(scale.i, "", "scaled_")  ## SAME ORDER AS names(l.scale)!!!!
    
    D <- readRDS(file.path(out.dir, paste0("euclidean_",  suffix, name.task.i, "_", name.glm.i, ".RDS")))
    
    
    cl <- makeCluster(n.cores/2)
    registerDoParallel(cl)
    time.start <- Sys.time()
    l.parcel <- foreach(
      parcel.i = seq_along(parcellation$key), 
      .packages = "boot",
      .final = function(x) setNames(x, parcellation$key)
      ) %dopar% {
    
    # for (parcel.i in seq_along(parcellation$key)) {
      # parcel.i = 1
      
      D_i <- D[, , parcel.i, ]
      
      lt <- lower.tri(diag(length(subjs)))
      ut <- upper.tri(diag(length(subjs)))
      
      ## get means and CIs:
      
      res.u <- boot(D_i[, , "univariate"], idi, R = nresamp, parallel = "multicore")
      res.m <- boot(D_i[, , "multivariate"], idi, R = nresamp, parallel = "multicore")
      res.contrast <- boot(D_i, idi_contrast, R = nresamp, parallel = "multicore")
      
      istats <- as.data.frame(
        rbind(
          boot.ci(res.u, type = "bca")$bca[4:5],
          boot.ci(res.m, type = "bca")$bca[4:5],
          boot.ci(res.contrast, type = "bca")$bca[4:5]
        )
      )
      names(istats) <- c("lb", "ub")
      istats$variable <- c("univariate", "multivariate", "contrast")
      istats$estimate <- c(res.u$t0, res.m$t0, res.contrast$t0)
      
      ## get p-values:
      
      istats$p <- 
        c(
          freqp(res.u$t, "greater"), 
          freqp(res.m$t, "greater"),
          freqp(res.contrast$t, "two.sided")
        )
      
      istats
      
    }
    stopCluster(cl)
    
    l.scale[[scale.i]] <- bind_rows(l.parcel, .id = "parcel")
    
  }
  
  
  l[[glm.i]] <- bind_rows(l.scale, .id = "transform")
  print(paste0(glm.i, " of ", nrow(glminfo), " glms done\n"))
  
  
}
(time.end <- Sys.time() - time.start)

## combine and save:

results <- bind_rows(l, .id = "task")

fwrite(results, file.path(out.dir, paste0("idi_group.csv")))


(time.run <- Sys.time() - time.start)




