if (interactive()) {
  
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
  
}


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



## setup ----



## functions

get_good_verts <- function(m, vert_names_keep) {
  ## works on inverse covariance matrices.
  ## for a given ROI*subj, the invcov matrix for each task may be estimated from a different set of vertices.
  ## (e.g., due to having no signal in some vertices for some tasks).
  ## this function makes it possible to subset each invcov matrix to the rows/cols that correspond to the 
  ## intersection of these sets of vertices.
  ##
  ## NB: vert_names_keep should be a CHARACTER vector.
  
  vert_names <- as.character(attr(m, "which.vert"))
  dimnames(m) <- list(row = vert_names, col = vert_names)
  m[vert_names_keep, vert_names_keep]
  
}


## weight matrices for averaging across folds

## fold: task*run; hold-one-task-out cross-validation
## avg is post-multiplied to matrix of voxel-by-taskrun matrix of beta coefficients
## avg creates task axes in "hold-one-task-out" manner, by averaging across task*runs of held-out task.

avg <- matrix(1, length(taskruns), length(tasks))  ## averages
dimnames(avg) <- list(train = taskruns, test = tasks)
for (name_task_i in tasks) avg[grep(name_task_i, rownames(avg)), name_task_i] <- 0
avg <- sweep(avg, 2, colSums(avg), "/")



## progress bar

n_iter <- length(subjs)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n_iter, clear = FALSE, width = 120
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
  
  B[, , name_task_i, ] <- betas_i_hi - betas_i_lo  ## TODO: move to subsequent loop so can standardize.
  
}




## estimate projections ----


p <- 
  array(
    NA,
    dim = c(length(tasks), length(rois), length(subjs), 2),
    dimnames = list(
      task = tasks, parcel = names(rois), 
      subj = subjs, run = c("run1", "run2")
      )
  )


for (subj_i in seq_along(subjs))  {
  # subj_i = 1
  
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
    W <- lapply(W, get_good_verts, vert_names_keep = as.character(good_verts))
    
    U <- B_roi_i[good_verts, , ]  ## pattern matrix
    U <- apply(U, 2:3, scale)  ## z-score standardize patterns per task*run
    # U <- U / sqrt(nrow(U))  ## scale by number of vertices (for comparability across ROIs)
    
    
    ## prewhiten each axes per task
    
    U_w <- U
    for (task_i in seq_along(tasks)) {
      # task_i = 1
      U_w[, task_i, ] <- crossprod(W[[task_i]], U[, task_i, ])
    }
    
    # par(mfrow = c(1, 2))
    # image(cor(U))
    # image(cor(U_w))
    
    
    ## reshape to 2D array:
    U_w <- U_w
    dim(U_w) <- c(nrow(U_w), length(tasks)*2)  ## concatenate runs columnwise
    colnames(U_w) <- combo_paste(tasks, c("run1", "run2"))
    U_w <- U_w[, taskruns]
    
    ## train (get axes):
    
    a <- U_w %*% avg  ## hi-lo axes
    colnames(a) <- tasks
    a <- scale2unit(a)  ## scale each axis by it's length (root sum of squares)
    
    ## project onto axes:
    
    for (task_i in seq_along(tasks)) {
      
      a_i <- a[, task_i]
      u_i <- U[, task_i, ]
      
      ## divide by number of vertices for comparison across ROIs:
      p[task_i, roi_i, subj_i, ] <- crossprod(a_i, u_i) / nrow(U)
      
    }

  }
  
  pb$tick()  ## progress bar
  
}



## save ----

if (!dir.exists(here("out", "icc"))) dir.create(here("out", "icc"))
saveRDS(
  p, 
  here(
    "out", "icc", 
    paste0(
      "hiloaxis-projections", 
      "_parc-schaefer07",
      "_prew-concat-runs",
      "_stand-axes",
      ".RDS"
      )
    )
  )
