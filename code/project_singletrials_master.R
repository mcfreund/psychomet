#source(here::here("code", "project_singletrials_master.R"))


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
resid_type <- "wherr"  ## errts, wherr
# do_network <- TRUE  ## for noise mods
do_wthntsk <- TRUE
do_btwntsk <- TRUE
do_stand <- TRUE
do_prew <- TRUE
prewtype <- "multi-alltask"  ## multi-alltask, multi-trials, "univ...
do_optwt <- FALSE


## exclude these guys for now?

subjs <- subjs_ub55[!subjs_ub55 %in% c("432332", "DMCC5820265", "DMCC9441378", "DMCC8260571")]

## glms to get axes from:

glminfo <- data.table(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name_glm = c(
    "baseline_Cues_EVENTS_censored_shifted",
    "baseline_cueletnum_EVENTS_censored_shifted",
    # "baseline_CongruencySwitch_EVENTS_censored_shifted",
    "baseline_ListLength_EVENTS_censored_shifted",
    "baseline_Congruency_EVENTS_censored_shifted"
  )
)
lev1 = list(
  Axcpt = "BX",
  Cuedts = c("num_odd_vow", "num_eve_con", "let_odd_vow", "let_eve_con"),
  # Cuedts = c("InConSwitch", "InConRepeat"),
  Stern = "LL5RN",
  Stroop = "biasInCon"
)
lev2 = list(
  Axcpt = "BY",
  Cuedts = c("num_odd_con", "num_eve_vow", "let_odd_con", "let_eve_vow"),
  # Cuedts = c("ConSwitch", "ConRepeat"),
  Stern = "LL5NN",
  Stroop = "biasCon"
)





## execute ----


## roi_set: networkplusmd, parcel
## prewhiten
## residtype
## contrastcontrast <- "ons"  ## "hilo", "ons"


for (roi_set in c("network07", "parcel")) {
  
  if (roi_set == "network07") {

    rois <- split(key_schaefer$parcel, key_schaefer$network)
    
  } else if (roi_set == "parcel") {
    if (prewtype == "multi-alltask") next
    rois <- split(key_schaefer$parcel, key_schaefer$parcel)
    
  }
  
  
	for (contrast in c("hilo", "ons")) {

		if (contrast == "hilo") {

			do_center <- TRUE  ## centers patterns at mean-pattern
			source(here("code", "project_singletrials.R"))

		} else if (contrast == "ons") {

		  do_center <- FALSE
			source(here("code", "project_singletrials.R"))

		}
	  
	  print("done")
	  
	}
  

}





#do_prew <- TRUE  ## prewhiten activity patterns?
#prewtype <- "multi"  
#do_stand <- TRUE  ## standardize contrast patterns prior to computing projections/prewhitening
#do_optwt <- FALSE  ## compute weighting of tasks based on estimated reliability
#do_wthntsk <- TRUE  ## calculate within-task (cross-run) projections; only for prewhitened


# for (do_prew in c(TRUE, FALSE)) {
# 	## prewhiten activity patterns?
# 	
# 	for (do_stand in c(TRUE, FALSE)) {
# 	## standardize contrast patterns prior to computing projections/prewhitening
# 			
# 		for (do_optwt in c(TRUE, FALSE)) {
# 		## compute weighting of tasks based on estimated reliability
# 			
# 			if (do_prew) {
# 			
# 				do_wthntsk <- TRUE  ## calculate within-task (cross-run) projections
# 				for (prewtype in c("multi", "univ")) source(here("code", "est_hiloaxis.R"))
# 				
# 			} else {
# 			
# 				do_wthntsk <- FALSE
# 				source(here("code", "est_hiloaxis.R"))
# 				
# 			}
# 			
# 		}
# 		
# 	}
# 	
# }


## scratch ----

## for "action axis" (L-R)

# tasks <- c("Axcpt", "Cuedts", "Stern")
# taskruns <- sort(combo_paste(tasks, c("run1", "run2")))
# 
# glminfo <- data.table(
#   task = c("Axcpt", "Cuedts", "Stern"),
#   name_glm = c(
#     "baseline_Cues_EVENTS_censored_shifted",
#     "baseline_cueletnum_EVENTS_censored_shifted",
#     "baseline_ListLength_EVENTS_censored_shifted"
#   )
# )
# lev1 = list(  ## left button press ("button 1")
#   Axcpt = "BY",
#   Cuedts = c("num_eve_vow", "num_eve_con", "let_odd_vow", "let_eve_vow"),  ## left/index: num_eve, let_vow
#   Stern = c("LL5NP", "not5NP")
# )
# lev2 = list(  ## right button press ("button 2")
#   Axcpt = "AX",
#   Cuedts = c("num_odd_vow", "num_odd_con", "let_odd_con", "let_eve_con"),  ## right/middle: num_odd, let_con
#     Stern = c("LL5NN", "not5NN")
# )