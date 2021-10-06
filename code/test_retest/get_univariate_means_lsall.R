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
library(gifti)

source(here("code", "_constants.R"))
source(here("code", "_atlases.R"))
source(here("code", "_funs.R"))
source(here("code", "_read_behav.R"))


## input vars ----

atlas <- "schaefer"
subjs <- subjs_wave12
do_waves <- c(1, 2)
tasks <- "Stroop"
glmname = "lsall_2rpm"
rois <- key_schaefer$parcel



## execute ----

waves <- waves[do_waves]
roi_set <- "parcel400"


## read betas


res <- read_results(
  
  waves = waves, tasks = tasks, sessions = sessions, subjs = subjs, glmname = "lsall_2rpm",
  
  filename_fun = function(name_subj_i, ...) paste0("STATS_", name_subj_i, "_", c("L", "R"), "_REML.func.gii"),
  
  read_fun = function(fname) abind(lapply(fname, read_gifti2matrix), along = 2)
  
  )

# Reduce(union, lapply(res, dim))


## get regressor labels

fnames <- 
  c(
    here("out", "glms", "130518", "RESULTS", "Stroop", "baseline_lsall_2rpm_wave1", "STATS_130518_R_REML.func.gii"),
    here("out", "glms", "130518", "RESULTS", "Stroop", "proactive_lsall_2rpm_wave1", "STATS_130518_R_REML.func.gii"),
    here("out", "glms", "130518", "RESULTS", "Stroop", "reactive_lsall_2rpm_wave1", "STATS_130518_R_REML.func.gii")
  )


subbrick_labels <- lapply(fnames, get_subbrick_labels)
stopifnot(identical(subbrick_labels[[1]], subbrick_labels[[2]]))  ## check if proactive == reactive

## extract coef vectors

is_coef_baspro <- grepl("_Coef", subbrick_labels[[1]])
res_baspro <- lapply(
  res[grepl("baseline|proactive", names(res))],
  function(x) x[is_coef_baspro, ]
)

is_coef_rea <- grepl("_Coef", subbrick_labels[[3]])
res_rea <- lapply(
  res[grepl("reactive", names(res))],
  function(x) x[is_coef_baspro, ]
)

res <- c(res_baspro, res_rea)

res <- res[sort(names(res))]




## take means

means <- enlist(names(res))
for (el_i in seq_along(res)) {
  
  res_i <- res[[el_i]]
  
  
  means[[el_i]] <- vapply(
    
    seq_along(rois),
    
    function(roi_i) {
      
      which_parcels <- match(rois[[roi_i]], key_schaefer$parcel)  ## works with both network and parcel level
      is_roi <- schaefer10k %in% which_parcels
      rowMeans(res_i[, is_roi])
      
    },
    
    numeric(nrow(res_i))
    
  )
  
  
  
}




## write to file


filenames <- get_filenames_results(
  
  waves = waves, tasks = tasks, sessions = sessions, subjs = subjs, glmname = "lsall_2rpm",
  
  filename_fun = function(...) paste0("means-trials", "_", atlas, "-", roi_set, "_lsall.csv")

)

lapply(names(filenames), function(x) fwrite(means[[x]], filenames[x]))

