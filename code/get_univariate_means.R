library(here)
library(here)
library(dplyr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mikeutils)
library(progress)

source(here("code", "_constants.R"))
source(here("code", "_atlases.R"))
source(here("code", "_funs.R"))

## input vars ----

atlas <- "schaefer"
do_network <- TRUE

rois <- split(key_schaefer$parcel, key_schaefer$network)

subjs <- subjs_ub55[!subjs_ub55 %in% "432332"]

glminfo <- data.table(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name_glm = c(
    "baseline_Cues_EVENTS_censored_shifted",
    "baseline_CongruencySwitch_EVENTS_censored_shifted",
    "baseline_ListLength_EVENTS_censored_shifted",
    "baseline_Congruency_EVENTS_censored_shifted"
  )
)

hi = list(
  Axcpt = "BX",
  Cuedts = c("InConInc", "InConNoInc"),
  Stern = "LL5RN",
  Stroop = "biasInCon"
)
lo = list(
  Axcpt = "BY",
  Cuedts = c("ConInc", "ConNoInc"),
  Stern = "LL5NN",
  Stroop = "biasCon"
)




## wrangle beta contrasts into array B ----

B <- 
  array(
    NA,
    dim = c(n_vert, length(subjs), length(tasks), 2),
    dimnames = list(vertex = NULL, subj = subjs, tasks = tasks, run = c("run1", "run2"))
  )
for (glm_i in seq_len(nrow(glminfo))) {
  # glm_i = 2
  
  name_glm_i <- glminfo[glm_i]$name_glm
  name_task_i <- glminfo[glm_i]$task
  
  ## read betas:
  
  betas_i <- 
    readRDS(
      here::here(
        "..", "ub55", "out", "glms", 
        paste0("betas_", name_task_i, "_", name_glm_i,  ".RDS")
        )
  )
  
  ## remove subj with missing data, get target TRs:
  
  betas_i <- betas_i[, , target_trs[[glm_i]], dimnames(betas_i)$subj %in% subjs, ]
  betas_i <- aperm(betas_i, c("vertex", "subj", "reg", "run", "tr"))
  
  ## contrast and save:
  
  betas_i_hi <- rowMeans(betas_i[, , hi[[glm_i]], , ], dims = 3) ## ave across TR, relevant regs
  betas_i_lo <- rowMeans(betas_i[, , lo[[glm_i]], , ], dims = 3) ## ave across TR, relevant regs
  
  B[, , name_task_i, ] <- betas_i_hi - betas_i_lo
  
}




## estimate projections ----

m <- 
  array(
    NA,
    dim = c(length(tasks), length(rois), length(subjs), 2),
    dimnames = list(
      task = tasks, roi = names(rois), 
      subj = subjs, run = c("run1", "run2")
      )
  )


for (subj_i in seq_along(subjs))  {
  # subj_i = which(subjs == "DMCC6418065")
  
  res <- enlist(names(rois))
  
  name_subj_i <- subjs[subj_i]
  B_subj_i <- B[, subj_i, , ]
  
  
  ## read noise models (whitening matrices):
  
  noise <- enlist(tasks)
  for (task_i in seq_along(tasks)) {
    noise[[task_i]] <- readRDS(
      here(
        "..", "ub55", "out", "invcov", name_subj_i,
        paste0(
          "invcov_", tasks[task_i], "_baseline_aggressive1_EVENTS_censored_shifted", 
          "_est-concat", 
          "_parc-", switch(do_network + 1, "parcels400", "network7"), 
          ".RDS"
        )
      )
    )
  }
  
  noise <- bind_rows(noise, .id = "task")
  
  
  ## estimate projections:
  
  for (roi_i in seq_along(rois)) {
    # roi_i = 1
    
    name_roi_i <- names(rois)[roi_i]
    
    
    ## mask:
    
    which_parcels <- match(rois[[roi_i]], key_schaefer$parcel)  ## works with both network and parcel level
    is_roi <- schaefer10k %in% which_parcels
    B_roi_i <- B_subj_i[is_roi, , ]
    
    ## get good verts:
    
    W <- setNames(noise$invcov[noise$roi == name_roi_i], noise$task[noise$roi == name_roi_i])
    good_verts <- Reduce(intersect, lapply(W, function(x) attr(x, "which.vert")))
    if (is.null(good_verts)) next
    
    U <- B_roi_i[good_verts, , ]  ## pattern matrix    
    
	## univariate (means)
	
	m[, roi_i, subj_i, ] <- colMeans(B_roi_i)
	
	
  }
    
}



## save ----

saveRDS(m, here("out", "icc", paste0("means", "_parc-schaefer07", ".RDS")))
