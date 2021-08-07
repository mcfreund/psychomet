#source(here::here("code", "icc_trial_level", "get_univariate_means_singletrial_master.R"))



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

## exclude these guys for now?

subjs <- subjs_ub55[!subjs_ub55 %in% c("432332", "DMCC5820265", "DMCC9441378", "DMCC8260571")]



## execute ----



for (roi_set in c("network07", "parcel400")) {
  
  if (roi_set == "network07") {

    rois <- split(key_schaefer$parcel, key_schaefer$network)
    
  } else if (roi_set == "parcel400") {
    
    rois <- split(key_schaefer$parcel, key_schaefer$parcel)
    
  }
  
	for (resid_type in c("errts", "wherr")) {
	  source(here::here("code", "icc_trial_level", "get_univariate_means_singletrial_master.R"))
	}
  

}
