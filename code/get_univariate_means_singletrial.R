
means <- enlist(subjs)
for (subj_i in seq_along(subjs)) {
  # subj_i = 44
  
  name_subj_i <- subjs[subj_i]

  
  ## read ----
  
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
  
  means_i <- enlist(names(rois))
  for (roi_i in seq_along(rois)) {
    # roi_i = 1
    
    name_roi_i <- names(rois)[roi_i]
    
    
    ## mask:
    
    which_parcels <- match(rois[[roi_i]], key_schaefer$parcel)  ## works with both network and parcel level
    is_roi <- schaefer10k %in% which_parcels
    
    resids_roi_i <- lapply(
      resids, 
      function(x) {
        x <- x[, is_roi, ]
        names(dimnames(x)) <- c("trial", "vertex", "run") 
        x <- aperm(x, c("trial", "run", "vertex"))
        rowMeans(x, dims = 2)
      }
    )
    
    resids_roi_i <- lapply(resids_roi_i, function(x) as.data.table(reshape2::melt(x))) %>% rbindlist(idcol = "task")
    
    means_i[[roi_i]] <- resids_roi_i
  
  }
  
  
  ## bind and return (end subj loop):
  
  means[[subj_i]] <- rbindlist(means_i, idcol = "roi")
  
  
}

means <- rbindlist(means, idcol = "subj")


## save ----

saveRDS(
  means, 
  here(
    "out", "icc", 
    paste0(
      "means-trials", 
      "_", atlas, "-", roi_set,
      "_resid-", resid_type,
      ".RDS"
    )
  )
)
  
