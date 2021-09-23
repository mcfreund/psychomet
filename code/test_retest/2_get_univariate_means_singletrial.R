#source(here::here("code", "test_retest", "get_univariate_means_singletrial_master.R"))


library(colorout)
library(here)
library(here)
library(dplyr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mikeutils)
library(ggplot2)
library(purrr)

source(here("code", "_constants.R"))
source(here("code", "_atlases.R"))
source(here("code", "_funs.R"))
source(here("code", "_read_behav.R"))


## input vars ----

atlas <- "schaefer"
subjs <- subjs_wave12
resid_type <- "errts"  ## "wherr" or "errts"
do_waves <- c(1, 2)


## execute ----

waves <- waves[do_waves]
# roi_set <- "parcel400"
# subj_i <- 1
# roi_i = 1
# wave_i = 1
# task_i = 1
# session_i = 1
# name_task_i <- tasks[task_i]
# name_wave_i <- waves[wave_i]
# name_session_i <- sessions[session_i]


for (roi_set in c("network07", "parcel400")) {

  if (roi_set == "network07") {
    
    rois <- split(key_schaefer$parcel, key_schaefer$network)
    
  } else if (roi_set == "parcel400") {
    
    rois <- split(key_schaefer$parcel, key_schaefer$parcel)
    
  }
  
  
  for (wave_i in seq_along(waves)) {
    
    name_wave_i <- waves[wave_i]
    
    for (task_i in seq_along(tasks)) {
      
      name_task_i <- tasks[task_i]
      
      for (session_i in seq_along(sessions)) {
        
        name_session_i <- sessions[session_i]
        
        for (subj_i in seq_along(subjs)) {
          
          
          name_subj_i <- subjs[subj_i]
          
          dir_glm <- 
            here("out", "glms", name_subj_i, "RESULTS", name_task_i, paste0(name_session_i, "_null_2rpm_", name_wave_i))
          
          fname_resids <- here(dir_glm, paste0(resid_type, "_trials_target_epoch.RDS"))
          resids <- readRDS(fname_resids)
          
          means <- vapply(
            
            seq_along(rois),
            
            function(roi_i) {
              
              which_parcels <- match(rois[[roi_i]], key_schaefer$parcel)  ## works with both network and parcel level
              is_roi <- schaefer10k %in% which_parcels
              rowMeans(resids[, is_roi])
              
            },
            
            numeric(nrow(resids))
            
          )
          
          dimnames(means) <- list(trial = paste0("trial", 1:nrow(means)), roi = names(rois))
          means <- as.data.table(means)
          
          
          filename <- here(dir_glm, paste0("means-trials", "_", atlas, "-", roi_set, "_resid-", resid_type, ".csv"))
          fwrite(means, filename)
          
          
        }
        
      }
      
    }
    
  }
  
  print("done")
  
	
}
