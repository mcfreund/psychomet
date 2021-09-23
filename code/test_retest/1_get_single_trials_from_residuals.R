## setup ----

library(colorout)
library(here)
library(dplyr)
library(tibble)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mikeutils)
library(ggplot2)
library(purrr)
library(gifti)

source(here("code", "_constants.R"))
source(here("code", "_atlases.R"))
source(here("code", "_funs.R"))
source(here("code", "_read_behav.R"))


subjs <- subjs_wave12
do_waves <- c(1, 2)
resid_type <- "errts"  ## "wherr" or "errts"; probably don't want to use wherr for this purpose.

## iterators for dev:
# wave_i <- 1
# task_i <- 4
# session_i <- 1
# subj_i <- 1


## run ----

## 1. build averaging matrix A, of dimension trial by TR. ----
## when applied to BOLD time series E (TR by vertex), A aggregates across target TRs within each trial.
## that is, Z = A %*% E  is a trial-by-vertex matrix.
## To do this, rows (TRs) of A must be >0 when the TR belongs to the target window of a given trial, 0 otherwise.
## Further, rows must sum to one.

waves <- waves[do_waves]

A_list <- enlist(combo_paste(waves, tasks, sessions, subjs))

for (wave_i in seq_along(waves)) {
  
  name_wave_i <- waves[wave_i]
  
  dir_image <- file.path("/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS", wavedir_image[wave_i], "fMRIPrep_AFNI_ANALYSIS")
  dir_evts <- file.path("/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/EVTS", wavedir_evts[wave_i])
  
  
  for (task_i in seq_along(tasks)) {
    
    name_task_i <- tasks[task_i]
    
    for (session_i in seq_along(sessions)) {
      
      name_session_i <- sessions[session_i]
      
      n_tr <- n_trs[paste0(name_task_i, "_", name_session_i)]
      n_trial <- n_trialspr[paste0(name_task_i, "_", name_session_i)]*2
      
      
      ## read in stim times, censor list (all subjs):
      
      fname.stimtimes <-
        file.path(
          dir_evts, subjs, "evts", 
          paste0(subjs, "_", name_task_i, "_", name_session_i, "_allTrials.txt")
          )
      stimtimes <- lapply(fname.stimtimes, function(x) as.matrix(fread(x, header = FALSE, sep = " ")))
      # stimtimes <- lapply(fname.stimtimes, function(x) as.matrix(read.table(x, header = FALSE, sep = " ")))
      
      fname.censor <-
        file.path(dir_image, subjs,  "INPUT_DATA", name_task_i, name_session_i, "movregs_FD_mask.txt")
      censor_list <- lapply(fname.censor, fread, header = FALSE, colClasses = "integer")
      
      
      ## build averaging matrix A:
      
      for (subj_i in seq_along(subjs)) {
        
        onsets <- c(
          stimtimes[[subj_i]][1, ], 
          stimtimes[[subj_i]][2, ] + 1.2*n_tr/2    ## shift run 2 onsets by duration (in seconds) of run 1
          )
        
        if (sum(duplicated(onsets)) > 0) {
          print(
            paste0(
              "identical onsets: ", 
              subjs[subj_i], " ", name_task_i, " ", name_session_i, ", ", sum(duplicated(onsets))
              )
            )
        }
        overlapping_trials <- which(duplicated(onsets) | duplicated(onsets, fromLast = TRUE))
        intervals <- floor(onsets / 1.2) + 1  ## gives image (interval) to which onset belongs
  
        
        ## build indicator matrix (dummy):
        
        dummy <- matrix(0, nrow = n_tr, ncol = n_trial, dimnames = list(tr = NULL, trial = NULL))
        for (target_tr_i in seq_along(target_trs[[name_task_i]])) {
          
          ## assign unique value to each onset (i.e., each trial)
          
          target_tr_val <- target_trs[[name_task_i]][target_tr_i] + 1 ## plus 1 because of TENTzero/afni 0 based index
          trial_nums <- seq_along(intervals)

          tr_grid <- numeric(n_tr)
          tr_grid[intervals + target_tr_val] <- trial_nums
          
          ## set previous, so built in model.mat; for censoring overlapping trials
          tr_grid[which(tr_grid == overlapping_trials[2]) - 1] <- overlapping_trials[1]
  
          dummy <- dummy + model.matrix(~ as.factor(tr_grid))[, -1]  ## remove intercept
          
        }
        
        
        ## censor:
        
        dummy[, overlapping_trials] <- 0  ## censor overlapping trials
        dummy[censor_list[[subj_i]] < 1, ] <- 0  ## censor TRs
        
        ## scale:
        
        divisor <- colSums(dummy)
        dummy[, divisor < 1] <- NA ## when this is zero, all TRs of the trial were censored due to motion
        if (any(divisor < 1)) print("got one!")  ## just to look
        A <- sweep(dummy, 2, divisor, "/")  ## normalize elements of each column (trial) to sum to one
        
        ## save:
        
        nm <- paste0(waves[wave_i], "_", name_task_i, "_", name_session_i, "_", subjs[subj_i])
        A_list[[nm]] <- A
        # image(A)
        
      }
      
      
    }
  }
}



## <<<<<< BEGIN CHECKING DUMMY MATRIX >>>>>>

# read in design matrices (all subjs):
# fname.xmat <-
#   file.path(
#     dir_image, subjs,
#     "1TRpK_SURFACE_RESULTS", name_task_i, paste0(name_session_i, "_ON_MIXED_censored"), "X.xmat_L.1D"
#     )
# 
# ## NEED TO FORCE PYTHON2 USE
# # https://stackoverflow.com/questions/7237415/python-2-instead-of-python-3-as-the-temporary-default-python
# # mkdir ~/bin
# # PATH=~/bin:$PATH
# # ln -s /usr/bin/python2 ~/bin/python
# # rm ~/bin/python
# X <- read_xmat(fname.xmat[subj_i])
# X <- X > 0.9  ## threshold (b/c some non-'active' TRs will nevertheless have >0 value, though quite small)
# class(X) <- "numeric"
# is_censored <- is_equal(rowSums(X), 0)
# sum(is_censored)
# xmat_intervals <- unlist(apply(X[, c("ON_TRIALS#8", "ON_TRIALS#9", "ON_TRIALS#10")] > 0, 2, which), use.names = FALSE)
# onset_intervals <- unlist(apply(dummy > 0, 2, which))
# sum(!is_equal(sort(xmat_intervals), sort(onset_intervals)))
# sort(setdiff(xmat_intervals, onset_intervals))
# sort(setdiff(onset_intervals, xmat_intervals))
# which(censor_list[[subj_i]] < 1)
# length(xmat_intervals)
# length(onset_intervals)
# which(X[, "ON_TRIALS#0"] > 0) == (intervals + 1)
# (X[, "ON_TRIALS#10"] > 0) == (tr_grid > 0)
# m <- cbind(
#   # xmat = which(X[, "ON_TRIALS#0"] > 0),
#   # intervals + 1  ## plus 1 because of TENTzero
#   xmat = X[, "ON_TRIALS#1"] > 0,
#   own = tr_grid > 0
#   )
# image(m)
# plot(m)
# abline(0, 1)


## <<<<<< END CHECKING DUMMY MATRIX >>>>>>




## <<<<<< BEGIN CHECKING WITH PREVIOUS IMPLEMENTATION >>>>>>

A_list_ub55 <- readRDS(here("..", "ub55", "A_list.RDS"))
subjs_common <- intersect(A_list_ub55$subj, subjs)
A_list_ub55 <- A_list_ub55 %>% filter(subj %in% subjs_common)

subj_val = "448347"
task_val = "Stroop"

named_vector <- function(type, nms) setNames(vector(type, length(nms)), nms)

r <- named_vector("logical", combo_paste(tasks, subjs_common))
for (task_val in tasks) {

  for (subj_val in subjs_common) {

    a1 <- A_list_ub55 %>% filter(subj == subj_val, task == task_val) %>% select(mat)
    a1_run1 <- a1[1, "mat"][[1]][[1]]
    a1_run2 <- a1[2, "mat"][[1]][[1]]

    a1 <- rbind(
      cbind(a1_run1, array(0, dim = dim(a1_run2))),
      cbind(array(0, dim = dim(a1_run1)), a1_run2)
    )

    nm <- paste0("wave1_", task_val, "_baseline_", subj_val)
    a2 <- t(A_list[[nm]])

    r[[paste0(task_val, "_", subj_val)]] <- isTRUE(all.equal(c(a1), c(a2)))
    # rowSums(is.na(a1)) > 0
    # a2[rowSums(is.na(a1)) > 0, ]


  }

}

## only diffs seem to be final trials at end of run, which sometimes are NA in original implementation.
## checking the censor lists, these trials/trs are not censored, so they should be in there.



## <<<<<< END CHECKING WITH PREVIOUS IMPLEMENTATION >>>>>>




# name_subj_i <- subjs[subj_i]
# name_wave_i <- waves[wave_i]
# name_task_i <- tasks[task_i]

cl <- makeCluster(n_core / 2)
registerDoParallel(cl)
res <- foreach(
  subj_i = seq_along(subjs),
  .final = function(x) setNames(x, subjs),
  .verbose = TRUE,
  .packages = c("data.table", "here", "mikeutils")
) %dopar% {

# for (subj_i in seq_along(subjs)) {
  
  name_subj_i <- subjs[subj_i]
  
  for (wave_i in seq_along(waves)) {
    
    name_wave_i <- waves[wave_i]
  
    for (task_i in seq_along(tasks)) {
      
      name_task_i <- tasks[task_i]
      
      for (session_i in seq_along(sessions)) {
        
        
        
        ## get constants, averaging matrix:
        
        name_session_i <- sessions[session_i]
        name_glm <- paste0(name_session_i, "_", "null_2rpm_", waves[wave_i])
        
        dir_glm <- here("out", "glms", name_subj_i, "RESULTS", name_task_i, name_glm)
        
        n_tr <- n_trs[paste0(name_task_i, "_", name_session_i)]
        n_trial <- n_trialspr[paste0(name_task_i, "_", name_session_i)]*2
        
        nm <- paste0(waves[wave_i], "_", name_task_i, "_", name_session_i, "_", subjs[subj_i])
        A <- A_list[[nm]]
        
        
        ## load data:
        
        eps_name <- here(dir_glm, paste0(resid_type, "_", name_subj_i, "_", c("L", "R"), "_REML.func.gii"))  ## L, R
        if (any(!file.exists(eps_name))) stop("no file!")
        
        E <- cbind(read_gifti2matrix(eps_name[1]), read_gifti2matrix(eps_name[2]))  ## L, R
        
        dims_bad <- any(dim(E) != c(n_tr, n_vert))
        if (dims_bad) stop ("bad dims: error time-series")
        
        ## average and save:
        
        B <- crossprod(A, E)  ## trial by vertex matrix B
        # image(B)
        saveRDS(B, here(dir_glm, paste0(resid_type, "_trials_target_epoch.RDS")))
        
        
        
      }
      
    }
    
  }
  
}
stopCluster(cl)





## initialize list for holding A matrices:

# A_list <- expand.grid(
#   task = tasks,
#   subj = subjs,
#   run = 1:2,
#   mat = vector("list", 1)
# )
# A_list <- as_tibble(A_list)
# 
# 
# for (task.i in seq_along(tasks)) {
#   # task.i = 1
#   
#   ## get event regressors:
#   
#   X <- X_list[[task.i]]
#   X <- X[, -grep("Pol|block|movregs", colnames(X)), , ]
#   X <- X > 0.9  ## threshold (b/c some non-'active' TRs will nevertheless have >0 value, though quite small)
#   
#   
#   ## get matrix that dummy codes for target TR windows:
#   
#   targets <- X[, grep(paste0("#", target_trs[[task.i]], "$", collapse = "|"), colnames(X)), , ]
#   
#   
#   ## create averaging matrix A per subj*run:
#   
#   for (subj.i in seq_along(subjs)) {
#     # subj.i <- 1; run.i <- 1
#     
#     name.subj.i <- subjs[subj.i]
#     
#     if (name.subj.i %in% c("432332", "DMCC5820265")) next
#     
#     for (run.i in 1:2) {
#       
#       
#       A <- matrix(
#         0, 
#         nrow = n_trialspr[task.i], ncol = n.trs[task.i]/2, 
#         dimnames = list(trial = paste0("trial_", 1:n_trialspr[task.i]), tr = NULL)
#       )
#       
#       targets_i <- targets[, , name.subj.i, run.i]
# 
#       
#       ## loop over rows (TRs) of targets_i, filling in relevant columns of A:
#       
#       trial.i <- 1  ## start trial counter
#       two_trials_counter <- 0  ## start multiple simultaneous trials couter (see below)
#       # for (tr.i in seq_len(nrow(targets_i))) {
#       tr.i <- 0
#       while (tr.i < nrow(targets_i)) {  ## while loop for control over iteration (can repeat if necessary)
#         tr.i <- tr.i + 1  ## iterate tr
#         # print(tr.i)
#         
#         is_target_tr <- any(targets_i[tr.i, ])
#         if (!is_target_tr) next
#         
#         condition_tr_i <- gsub("#[0-9]*$", "", colnames(targets_i)[targets_i[tr.i, ]])  ## get name of cond, rm knot #
#         knot_tr_i <- as.numeric(gsub(".*#", "", colnames(targets_i)[targets_i[tr.i, ]]))  ## get knot #, remove name of cond
#         
#         ## handling cases when multiple trials have overlapping target TRs (in Stroop):
#         
#         has_multiple_targets <- length(condition_tr_i) > 1
#         # if (has_multiple_targets) stop("bad TR: multiple regressors are in 'target' window!")
#         if (has_multiple_targets) {
#           
#           if (length(condition_tr_i) == 2) {
#             ## when a double count/overlapping target trs is *first* encountered, we want to take the earliest knot.
#             ## then we want to iterate the trial count, but hold the tr.i back one iteration, so the later knot at tr.i
#             ## would be picked up on the next iteration.
#             ## on the subsequent iteration, we don't want to hold the trial count back.
#             ## so, we need to track whether this is the earlier or later knot (via is_first_of_two).
#             ## seems the only subj with this issue is 672756...
# 
#             two_trials_counter <- two_trials_counter + 1
#             is_first_of_two <- two_trials_counter == 1  ## behavior only needs to change if two
#             
#             if (is_first_of_two) {
#               ind <- which.max(knot_tr_i)  ## take max b/c first knot of overlap will be a knot towards the end of the trial
#             } else {
#               ind <- which.min(knot_tr_i)
#             }
#             
#             condition_tr_i <- condition_tr_i[ind]
#             knot_tr_i <- knot_tr_i[ind]
#             
#             print(paste0("2 targets: ", tasks[task.i], " ", name.subj.i, ", tr ", tr.i))
# 
#           } else stop("cant handle this right now")
#           
#         } else {
#           
#           ## for resetting counters (following )
#           two_trials_counter <- 0
#           is_first_of_two <- two_trials_counter == 1  ## behavior only needs to change if true
#           
#         }
#         
#         ## mark trial and save info:
#         
#         A[trial.i, tr.i] <- 1  ## mark trial
#         ## label row with condition and trial number info:
#         # rownames(A)[trial.i] <- paste0(condition_tr_i, "_", formatC(trial.i, width = 3, format = "d", flag = "0"))
#         rownames(A)[trial.i] <- condition_tr_i
#         
#         if (is_first_of_two) tr.i <- tr.i - 1  ## re-do this iteration to catch the second trial with target at this TR
#         
#         is_final_knot <- knot_tr_i == max(target_trs[[task.i]])  ## iterate trial counter if final knot
#         if (is_final_knot) trial.i <- trial.i + 1
#         
#       }
#       
#       ## normalize rows to sum to one:
#       
#       divisor <- rowSums(A)
#       divisor[divisor == 0] <- NA ## when this is zero, means TRs of the trial were fully censored due to motion
#       A <- sweep(A, 1, divisor, "/")
#       
#       is_unexpeted_res <- any(!unique(rowSums(A)) %in% c(1, NA))  ## a final check
#       if (is_unexpeted_res) stop("unexpected result: rows do not sum to one")
#       
#       
#       ## save:
#       
#       which.el <- which(A_list$task == tasks[task.i] & A_list$subj == subjs[subj.i] & A_list$run == run.i)
#       A_list$mat[[which.el]] <- A
#       
#     }
#   }
#   
#   
#   
# }




## 2. apply A to residuals and save resulting trial-wise BOLD activity estimates. ----

# 
# cl <- makeCluster(n.core / 2)
# registerDoParallel(cl)
# res <- foreach(
#   subj.i = seq_along(subjs),
#   .final = function(x) setNames(x, subjs),
#   .verbose = TRUE,
#   .packages = c("data.table", "here", "mikeutils")
# ) %dopar% {
#   # subj.i = 1
#   
#   name.subj.i <- subjs[subj.i]
#   
#   if (name.subj.i == "432332") return()
#       
#   for (task.i in seq_along(tasks)) {
#     # task.i = 1
#     
#     name.task.i <- glminfo_null[task.i]$task
#     name.glm.i <- glminfo_null[task.i]$name.glm
#     
#     
# 
#     for (run.i in 1:2) {
#       # run.i = 1
#       
#       B <- array(
#         NA,
#         c(n_trialspr[task.i], n.vert, 2),
#         dimnames = list(trial = NULL, vertex = NULL)
#         )
#       
#       eps.name <- here(
#         "out", "glms", name.subj.i, "RESULTS", name.task.i, paste0(name.glm.i, "_", run.i),
#         paste0(resid_type, "_", name.subj.i, "_", run.i, "_", c("L", "R"), "_REML.func.gii")
#       )  ## LEFT then RIGHT
#       
#       if (any(!file.exists(eps.name))) return(NA)
#       
#       E <- cbind(
#         read_gifti2matrix(eps.name[1]),
#         read_gifti2matrix(eps.name[2])
#       )  ## LEFT then RIGHT
#       
#       dims.bad <- any(dim(E) != c(n.trs[name.task.i]/2, n.vert))
#       if (dims.bad) stop ("bad dims: error time-series")
#       
#       ii <- which(A_list$subj == subjs[subj.i] & A_list$task == tasks[task.i] & A_list$run == run.i)
#       A <- A_list$mat[[ii]]
#       
#       B <- A %*% E
#       ## store condition names as attribute so can be wrangled into array later:
#       attr(B, paste0("order", run.i)) <- rownames(A)
#       dimnames(B)$trial <- NULL
#       
#       saveRDS(
#         B,
#         here(
#           "out", "glms", name.subj.i, "RESULTS", name.task.i, paste0("baseline_null_", run.i), 
#           paste0(resid_type, "_trials_target_epoch.RDS")
#           )
#         )
#       
#     }
#       
#     
#   }
# 
#   
# } ## end subj loop
# stopCluster(cl)




# dir_image <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/DMCC_Phase3/fMRIPrep_AFNI_ANALYSIS"
# dir_evts <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/EVTS/DMCC3/"
# glminfo_null <- data.table(
#   expand.grid(
#     task = tasks,
#     session = sessions
#   )
# )
# glminfo_null$name <- paste0(glminfo_null$session, "_null_2rpm")
# 
# 
# glminfo_events <- data.table(
#   expand.grid(
#     task = tasks,
#     session = sessions
#   )
# )
# glminfo_events$name <- paste0(glminfo_events$session, "_ON_MIXED_censored")
