## execute this script from _master.R file


## TODO
## ONs
## prewhitening
## change singletrial averaging code:
##    - fix code to handle no-signal verts / check for no-signal verts!
##    - dplyr -> data.table




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





## estimate axes and wrangle into array B ----


betas <- 
  array(
    NA,
    dim = c(n_vert, length(subjs), length(tasks), 2, 2),
    dimnames = list(vertex = NULL, subj = subjs, tasks = tasks, trialtype = c("lev1", "lev2"), run = c("run1", "run2"))
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
  
  betas_i <- betas_i[, , target_trs[[glm_i]], dimnames(betas_i)$subj %in% subjs, , drop = FALSE]
  betas_i <- aperm(betas_i, c("vertex", "subj", "run", "tr", "reg"))  ## permute for averaging across reg and tr
  
  ## average across TR, relevant regs:
  
  betas_i_lev1 <- rowMeans(betas_i[, , , , lev1[[glm_i]], drop = FALSE], dims = 3)  ## keeps: vertex, subj, run
  betas_i_lev2 <- rowMeans(betas_i[, , , , lev2[[glm_i]], drop = FALSE], dims = 3)
  
  betas[, , name_task_i, "lev1", ] <- betas_i_lev1
  betas[, , name_task_i, "lev2", ] <- betas_i_lev2
  
}
rm(name_glm_i, name_task_i)




## estimate projections ----




cl <- makeCluster(n_cores / 2)
registerDoParallel(cl)
res <- foreach(
  subj_i = seq_along(subjs),
  .final = function(x) setNames(x, subjs),
  .verbose = TRUE,
  .packages = c("data.table", "here", "mikeutils", "dplyr", "abind")
  ) %dopar% {
# for (subj_i in seq_along(subjs))  {
  # subj_i = 44
  
  
  
  name_subj_i <- subjs[subj_i]
  betas_subj_i <- betas[, subj_i, , , ]
  
  
  ## read ----
  
  
  ## read noise models (whitening matrices):
  
  fname_noise <- here(
    "..", "ub55", "out", "invcov", name_subj_i,
    paste0(
      "invcov_", tasks, "_baseline_aggressive1_EVENTS_censored_shifted", 
      "_est-concat", 
      "_parc-", switch(do_network + 1, "parcels400", "network7"), 
      ".RDS"
    )
  )
  noise <- bind_rows(lapply(fname_noise, readRDS), .id = "task")

  
  
  ## read trials and trial order (residuals from null model)
  
  resids <- enlist(tasks)
  trialorders <- enlist(tasks)
  # resporders <- enlist(tasks)
  for (name_task_i in tasks) {
    # name_task_i = "Cuedts"
    
    fname_resids <- here(
      "..", "ub55", 
      "out", "glms", name_subj_i, "RESULTS", name_task_i, paste0("baseline_null_", 1:2), 
      "errts_trials_target_epoch.RDS"
    )  ## both runs
    resids_i <- lapply(fname_resids, readRDS)
    trialorders[[name_task_i]] <- do.call(cbind, lapply(resids_i, attr, "order"))  ## save trial orders
    resids[[name_task_i]] <- abind(resids_i, rev.along = 0)  ## bind into array; vert by trial by run
    
  }
  
  
  

  ## project ----
  
  
  if (do_wthntsk) proj_wthn_subj_i <- enlist(names(rois))
  if (do_btwntsk) proj_btwn_subj_i <- enlist(names(rois))
  
  for (roi_i in seq_along(rois)) {
    # roi_i = 1
    
    name_roi_i <- names(rois)[roi_i]
    
    
    
    ## mask:
    
    which_parcels <- match(rois[[roi_i]], key_schaefer$parcel)  ## works with both network and parcel level
    is_roi <- schaefer10k %in% which_parcels
    
    betas_roi_i <- betas_subj_i[is_roi, , , ]
    
    resids_roi_i <- lapply(
      resids, 
      function(x) {
        x <- x[, is_roi, ]
        names(dimnames(x)) <- c("trial", "vertex", "run") 
        x
        }
      )

    
    
    ## get good verts:
    
    # W <- setNames(noise$invcov[noise$roi == name_roi_i], noise$task[noise$roi == name_roi_i])
    # good_verts <- Reduce(intersect, lapply(W, function(x) attr(x, "which.vert")))
    # 
    # if (is.null(good_verts)) next
    # 
    # W <- lapply(W, get_good_verts, vert_names_keep = as.character(good_verts))
    # betas_roi_i <- betas_roi_i[good_verts, , , ]  ## pattern matrix
    # resids_roi_i <- lapply(resids_roi_i, function(x) x[, good_verts, ])
    
    
    
    ## center resids at mean of condition means ('mean pattern'):

    if (do_center) {
      for (name_task_i in tasks) {
        for (run_i in 1:2) {
          # name_task_i = "Cuedts"; run_i = 1
  
          is_missing_trial <- rowSums(is.na(resids_roi_i[[name_task_i]][, , run_i])) == ncol(resids_roi_i[[name_task_i]])
  
          # resids_roi_i[!is_missing_trial, , run_i]
  
          avgmat <- model.matrix(~ 0 + trialorders[[name_task_i]][!is_missing_trial, run_i])
          colnames(avgmat) <- gsub("^trialorders.*i]", "", colnames(avgmat))
          
          avgmat_lev1 <- avgmat[, lev1[[name_task_i]], drop = FALSE]
          avgmat_lev2 <- avgmat[, lev2[[name_task_i]], drop = FALSE]
  
          avgmat_lev1 <- rowSums(sweep(avgmat_lev1, 2, colSums(avgmat_lev1), "/")) / ncol(avgmat_lev1)
          avgmat_lev2 <- rowSums(sweep(avgmat_lev2, 2, colSums(avgmat_lev2), "/")) / ncol(avgmat_lev2)
  
          mean_pattern_mat <- (avgmat_lev1 + avgmat_lev2)/2
          mean_pattern <- c(crossprod(resids_roi_i[[name_task_i]][!is_missing_trial, , run_i], mean_pattern_mat))
  
          resids_roi_i[[name_task_i]][!is_missing_trial, , run_i] <-
            sweep(resids_roi_i[[name_task_i]][!is_missing_trial, , run_i], 2, mean_pattern)
  
        }
      }
    }
    
    
    ## estimate axes:
    
    if (contrast == "hilo") {
      
      axis_roi_i <- betas_roi_i[, , "lev1", ] - betas_roi_i[, , "lev2", ]  ## get axes:
      
    } else if (contrast == "ons") {
      
      axis_roi_i <- betas_roi_i[, , "lev1", ] + betas_roi_i[, , "lev2", ]  ## get axes:
      
    }
    # axis_roi_i <- apply(betas_roi_i[, , "lev1", ], 2:3, scale) - apply(betas_roi_i[, , "lev2", ], 2:3, scale)
    
    axis_roi_i <- apply(axis_roi_i, 2:3, scale2unit)
    
    
    ## within-task projections (btwn run):
    
    if (do_wthntsk) {

      projs_wthn <- enlist(tasks)
      for (task_i in seq_along(tasks)) {

        U <- resids_roi_i_center[[task_i]]
        trialorders_i <- trialorders[[task_i]]
        
        ## (scale by nvert to aid comparison across ROIs):
        proj1 <- (U[, , 1] %*% axis_roi_i[, task_i, 2]) / nrow(axis_roi_i)  ## test: run 1
        proj2 <- (U[, , 2] %*% axis_roi_i[, task_i, 1]) / nrow(axis_roi_i)  ## test: run 2
        
        projs1 <- data.table(value = c(proj1), trialtype = trialorders_i[, 1])
        projs2 <- data.table(value = c(proj2), trialtype = trialorders_i[, 2])
        
        projs <- bind_rows(projs1, projs2, .id = "run")
        projs <- projs %>% group_by(run) %>% mutate(trial.num = 1:n())
        
        projs_wthn[[task_i]] <- projs

      }
      
      projs_wthn <- bind_rows(projs_wthn, .id = "task")
      
    }
    
    
    
    ## btwn-task projections:
    
    if (do_btwntsk) {
      
      
      ## reshape axes for averaging across tasks
      
      dim(axis_roi_i) <- c(dim(axis_roi_i)[1], prod(dim(axis_roi_i)[2:3]))  ## concatenate runs columnwise
      colnames(axis_roi_i) <- combo_paste(tasks, c("run1", "run2"))
      axis_roi_i <- axis_roi_i[, rownames(avg)]  ## same order as avg matrix
      
      ## average across tasks (training folds):
      
      axis_roi_i <- axis_roi_i %*% avg  ## lev1-lev2 axes
      axis_roi_i <- scale2unit(axis_roi_i)
      
      ## test:
      
      projs_btwn <- enlist(tasks)
      for (task_i in seq_along(tasks)) {
        
        U <- resids_roi_i_center[[task_i]]
        trialorders_i <- trialorders[[task_i]]
        
        U <- aperm(U, c(2, 1, 3))
        dim(U) <- c(dim(U)[1], prod(dim(U)[2:3]))  ## concatenate runs columnwise

        projs <- crossprod(U, axis_roi_i[, task_i]) / nrow(axis_roi_i)
        
        projs <- data.table(
          value = c(projs), 
          trialtype = c(trialorders_i), 
          run = as.character(rep(1:2, each = nrow(projs) / 2))
          )
        
        projs <- projs %>% group_by(run) %>% mutate(trialnum = 1:n())
        
        projs_btwn[[task_i]] <- projs
        
      }
      
      projs_btwn <- bind_rows(projs_btwn, .id = "task")
      
      
    }
    
    
    
    ## save (end ROI loop):
    
    if (do_wthntsk) proj_wthn_subj_i[[name_roi_i]] <- projs_wthn
    if (do_btwntsk) proj_btwn_subj_i[[name_roi_i]] <- projs_btwn
    
  }
  
  
  
  ## bind and return (end subj loop):
  
  res_i <- enlist(c("wthntsk", "btwntsk"))
  if (do_wthntsk) res_i$wthntsk <- bind_rows(proj_wthn_subj_i, .id = "roi")
  if (do_btwntsk) res_i$btwntsk <- bind_rows(proj_btwn_subj_i, .id = "roi")
  
  
  res_i
   
  
}
stopCluster(cl)



d <- bind_rows(lapply(res, bind_rows, .id = "comparison"), .id = "subj")
d <- as.data.table(d)


if (do_btwntsk) {
  
  saveRDS(
    d[comparison == "btwntsk"], 
    here(
      "out", "icc", 
      paste0(
        "projs-trials", 
        "_axis-", contrast,
        "_cval-btwntsk", 
        "_", atlas, "_", roi_set,
        "_prew-", switch(do_prew + 1, "none", paste0("concat-runs-", prewtype)),
        "_stand-", switch(do_stand + 1, "none", "axes"),
        "_center-", switch(do_center + 1, "none", "pattern"),
        "_resid-", resid_type,
        ".RDS"
      )
    )
  )
  
}


if (do_wthntsk) {

  saveRDS(
    d[comparison == "wthntsk"], 
    here(
      "out", "icc", 
      paste0(
        "projs-trials", 
        "_axis-", contrast,
        "_cval-wthntsk", 
        "_", atlas, "_", roi_set,
        "_prew-", switch(do_prew + 1, "none", paste0("concat-runs-", prewtype)),
        "_stand-", switch(do_stand + 1, "none", "axes"),
        "_center-", switch(do_center + 1, "none", "pattern"),
        "_resid-", resid_type,
        ".RDS"
      )
    )
  )
  
}





# d %>%
#   
#   filter(comparison == "wthntsk", trialtype %in% c(unlist(lev1), unlist(lev2))) %>%
#   # filter(comparison == "wthntsk", trialtype %in% c("button_1", "button_2")) %>%
#   
#   group_by(subj, roi, task, trialtype) %>%
#   summarize(value = mean(value, na.rm = TRUE)) %>%
#   
#   mutate(is_lev1 = trialtype %in% unlist(lev1)) %>%
#   
#   ggplot(aes(trialtype, value)) +
#   
#   # stat_summary(fun.data = "mean_cl_boot", aes(color = trialtype)) +
#   stat_summary(fun.data = "mean_cl_boot", aes(color = ifelse(is_lev1, "firebrick2", "black"))) +
#   scale_color_identity() +
#   
#   facet_grid(vars(roi), vars(task), scales = "free_x")
# 
# 
# d %>%
#   
#   # filter(comparison == "wthntsk", trialtype %in% c(unlist(lev1), unlist(lev2))) %>%
#   filter(comparison == "btwntsk", trialtype %in% c("button_1", "button_2")) %>%
#   
#   group_by(subj, roi, task, trialtype) %>%
#   summarize(value = mean(value, na.rm = TRUE)) %>%
#   
#   # mutate(is_lev1 = trialtype %in% unlist(lev1)) %>%
#   
#   ggplot(aes(trialtype, value)) +
#   
#   stat_summary(fun.data = "mean_cl_boot", aes(color = trialtype)) +
#   # stat_summary(fun.data = "mean_cl_boot", aes(color = ifelse(is_lev1, "firebrick2", "black"))) +
#   # scale_color_identity() +
#   
#   facet_grid(vars(roi), vars(task), scales = "free_x")
# 
# 
# 
# # d %>%
# #   
# #   filter(comparison == "btwntsk", trialtype %in% c(unlist(lev1), unlist(lev2))) %>%
# #   
# #   group_by(subj, roi, task, trialtype) %>%
# #   summarize(value = mean(value)) %>%
# #   
# #   mutate(is_lev1 = trialtype %in% unlist(lev1)) %>%
# #   
# #   ggplot(aes(trialtype, value)) +
# #   
# #   stat_summary(fun.data = "mean_cl_boot", aes(color = ifelse(is_lev1, "firebrick2", "black"))) +
# #   scale_color_identity() +
# #   
# #   facet_grid(vars(roi), vars(task), scales = "free_x")
# 
#   
# d %>%
#   
#   # filter(comparison == "btwntsk", trialtype %in% c(unlist(lev1), unlist(lev2))) %>%
#   filter(comparison == "wthntsk", trialtype %in% c("button_1", "button_2")) %>%
#   
#   
#   group_by(subj, roi, task, trialtype) %>%
#   summarize(value = mean(value, na.rm = TRUE)) %>%
#   
#   group_by(roi, task, trialtype) %>% 
#   mutate(value = value/sd(value)) %>%
#   # mutate(is_lev1 = trialtype %in% unlist(lev1), value = value/sd(value)) %>%
#   
#   ggplot(aes(trialtype, value)) +
#   
#   stat_summary(fun = "mean", geom = "col", aes(color = trialtype)) +
#   
#   # stat_summary(fun = "mean", geom = "col", aes(color = ifelse(is_lev1, "firebrick2", "black"))) +
#   # stat_summary(fun.data = "mean_cl_boot", aes(color = ifelse(is_hi, "firebrick2", "black"))) +
#   # scale_color_identity() +
#   
#   facet_grid(vars(roi), vars(task), scales = "free")
# 
# 
# 
# 
# 
# 
# ## save res ----
# 
# 
# 
# 
# 
# 
# 
# 
# stats_wthn <- d %>%
#   
#   filter(comparison == "wthntsk", trialtype %in% c(unlist(lev1), unlist(lev2))) %>%
#   
#   group_by(subj, roi, task, trialtype) %>%
#   summarize(value = mean(value, is.na = FALSE)) %>%
#   mutate(is_lev1 = trialtype %in% unlist(lev1)) %>%
#   
#   group_by(subj, roi, task, is_lev1) %>%
#   summarize(value = mean(value)) %>%
#   
#   tidyr::pivot_wider(values_from = "value", names_from = "is_lev1", names_prefix = "is_lev1_") %>%
#   group_by(roi, task) %>%
#   mutate(value = is_lev1_TRUE - is_lev1_FALSE) %>%
#   
#   summarize(
#     stat = t.test(value)$statistic,
#     p = t.test(value)$p.value
#   ) %>%
#   group_by(task) %>%
#   mutate(p.fdr = p.adjust(p, "fdr"))
# 
# 
# stats_wthn %>% filter(p.fdr < 0.05)
# 
# 
# 
# stats_btwn <- d %>%
#   
#   filter(comparison == "btwntsk", trialtype %in% c(unlist(lev1), unlist(lev2))) %>%
#   
#   group_by(subj, roi, task, trialtype) %>%
#   summarize(value = mean(value, is.na = FALSE)) %>%
#   mutate(is_lev1 = trialtype %in% unlist(lev1)) %>%
#   
#   group_by(subj, roi, task, is_lev1) %>%
#   summarize(value = mean(value)) %>%
#   
#   tidyr::pivot_wider(values_from = "value", names_from = "is_lev1", names_prefix = "is_lev1_") %>%
#   group_by(roi, task) %>%
#   mutate(value = is_lev1_TRUE - is_lev1_FALSE) %>%
#   
#   summarize(
#     stat = t.test(value)$statistic,
#     p = t.test(value)$p.value
#   ) %>%
#   group_by(task) %>%
#   mutate(p.fdr = p.adjust(p, "fdr"))
# 
# 
# stats_btwn %>% filter(p.fdr < 0.05)
# 
# 
# 
# 
# 
# stat_wthn <- d %>%
#   
#   filter(comparison == "wthntsk", trialtype %in% c("button_1", "button_2")) %>%
#   
#   group_by(subj, roi, task, trialtype) %>%
#   summarize(value = mean(value, is.na = FALSE)) %>%
#   group_by(subj, roi, task, trialtype) %>%
#   summarize(value = mean(value)) %>%
#   
#   tidyr::pivot_wider(values_from = "value", names_from = "trialtype") %>%
#   group_by(roi, task) %>%
#   mutate(
#     value = button_1 - button_2,
#     # value = value / sd(value)
#   ) %>%
#   summarize(
#     stat = t.test(value)$statistic,
#     p = t.test(value)$p.value
#     ) %>%
#   group_by(task) %>%
#   mutate(p.fdr = p.adjust(p, "fdr"))
# 
# 
# stat_wthn %>% filter(p.fdr < 0.05)
# 
# 
# 
# 
# stat_btwn <- d %>%
#   
#   filter(comparison == "btwntsk", trialtype %in% c("button_1", "button_2")) %>%
#   
#   group_by(subj, roi, task, trialtype) %>%
#   summarize(value = mean(value, is.na = FALSE)) %>%
#   group_by(subj, roi, task, trialtype) %>%
#   summarize(value = mean(value)) %>%
#   
#   tidyr::pivot_wider(values_from = "value", names_from = "trialtype") %>%
#   group_by(roi, task) %>%
#   mutate(
#     value = button_1 - button_2,
#     # value = value / sd(value)
#   ) %>%
#   summarize(
#     stat = t.test(value)$statistic,
#     p = t.test(value)$p.value
#   ) %>%
#   group_by(task) %>%
#   mutate(p.fdr = p.adjust(p, "fdr"))
# 
# 
# stat_btwn %>% filter(p.fdr < 0.1)






## scratch ----








# rowSums(is_equal(betas_roi_i, 0))

# if (do_stand) U <- apply(U, 2:3, scale)  ## z-score standardize patterns per task*run
# U <- U / sqrt(nrow(U))  ## scale by number of vertices (for comparability across ROIs)



## prewhiten each axes per task

# U_w <- U
# if (do_prew) {
#   if (prewtype %in% c("multi", "univ")) {
#     if (prewtype == "univ") {  ## set off-diagonal to zero
#       for (task_i in seq_along(tasks)) W[[task_i]][row(W[[task_i]]) != col(W[[task_i]])] <- 0
#     }
#     for (task_i in seq_along(tasks)) U_w[, task_i, ] <- crossprod(W[[task_i]], U[, task_i, ])
#   }
# }

#par(mfrow = c(1, 2))
#image(cor(U))
#image(cor(U_w))
#cor(U_w[, , 1], U[, , 1])

## reshape to 2D array:

# U_w <- U_w
# dim(U_w) <- c(nrow(U_w), length(tasks)*2)  ## concatenate runs columnwise
# colnames(U_w) <- combo_paste(tasks, c("run1", "run2"))
# U_w <- U_w[, taskruns]


# if (do_optwt) {
#   ## TODO:
#   ## create convex combination for each task / hold-out fold, to use as weighting for "training".
#   
#   r <- 
#     cbind(
#       cor(U_w[, "Axcpt_run1"], U_w[, "Axcpt_run2"]),
#       cor(U_w[, "Cuedts_run1"], U_w[, "Cuedts_run2"]),
#       cor(U_w[, "Stern_run1"], U_w[, "Stern_run2"]),
#       cor(U_w[, "Stroop_run1"], U_w[, "Stroop_run2"])
#     )
#   r[r<0] <- 0
#   
#   avg[] <- 0
#   for (task_i in seq_along(tasks)) {  ## hold-out task (test)
#     lambda <- r[-task_i] / sum(r[-task_i])
#     names(lambda) <- tasks[-task_i]
#     for (task_j in tasks[-task_i]) {  ## training task
#       avg[grep(task_j, rownames(avg)), task_i] <- lambda[task_j] / 2
#     }
#   }
#   
#   if (any(is.nan(avg))) print(paste0(subjs[subj_i], " ", names(rois)[roi_i], "\n"))
#   #stop("bad combo.")
#   
# }








# projs %>%
#   
#   ggplot(aes(trialnum, value, fill = trialtype)) +
#   
#   geom_hline(yintercept = 0) +
#   geom_line(aes(group = run)) +
#   geom_point(size = 4, shape = 21, color = "white") +
#   
#   facet_grid(cols = vars(run)) +
#   scale_fill_brewer(type = "qual", palette = 2)
# 
# projs %>%
#   
#   ggplot(aes(trialtype, value, color = trialtype)) +
#   stat_summary(fun.data = "mean_cl_boot") +
#   scale_color_brewer(type = "qual", palette = 2)
# 
# projs %>%
#   filter(trialtype %in% c("BY", "BX")) %>%
#   mutate(is_hi = value > 0) %>%
#   ungroup %>%
#   select(trialtype, is_hi) %>%
#   table






# resp <- behav[[name_task_i]][subj == subjs[subj_i], c("run", "trial.num", "resp")] %>% arrange(run, trial.num)
# resporders[[name_task_i]] <- cbind(
#   paste0("button_", resp[run == 1]$resp),
#   paste0("button_", resp[run == 2]$resp)
# )

## check trial orders against behavioral data:

# resp <- behav[[name_task_i]][subj == subjs[subj_i], c("run", "trial.num", "resp", "cue", "number", "letter")]
# test$cue <- ifelse(test$cue == "n", "num", ifelse(test$cue == "l", "let", NA))
# test$number <- ifelse(test$number == "e", "eve", ifelse(test$number == "o", "odd", NA))
# test$letter <- ifelse(test$letter == "v", "vow", ifelse(test$letter == "c", "con", NA))
# test$condition <- paste0(test$cue, "_", test$number, "_", test$letter)
# is_mismatched_condition <- 
#   any(c(as.matrix(cbind(test[run == 1, c("condition")], test[run == 2, c("condition")])) != trialorders$Cuedts))
