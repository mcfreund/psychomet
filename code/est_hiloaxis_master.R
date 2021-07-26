#source(here::here("code", "est_hiloaxis_master.R"))

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
## each list element defines an ROI; parcels within elements define constituent parcels of ROI:
rois <- split(key_schaefer$parcel, key_schaefer$network)
do_network <- TRUE  ## network or parcel level?

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

subjs <- subjs_ub55[!subjs_ub55 %in% "432332"]



## execute ----


#do_prew <- TRUE  ## prewhiten activity patterns?
#prewtype <- "multi"  
#do_stand <- TRUE  ## standardize contrast patterns prior to computing projections/prewhitening
#do_optwt <- FALSE  ## compute weighting of tasks based on estimated reliability
#do_wthntsk <- TRUE  ## calculate within-task (cross-run) projections; only for prewhitened


for (do_prew in c(TRUE, FALSE)) {
	## prewhiten activity patterns?
	
	for (do_stand in c(TRUE, FALSE)) {
	## standardize contrast patterns prior to computing projections/prewhitening
			
		for (do_optwt in c(TRUE, FALSE)) {
		## compute weighting of tasks based on estimated reliability
			
			if (do_prew) {
			
				do_wthntsk <- TRUE  ## calculate within-task (cross-run) projections
				for (prewtype in c("multi", "univ")) source(here("code", "est_hiloaxis.R"))
				
			} else {
			
				do_wthntsk <- FALSE
				source(here("code", "est_hiloaxis.R"))
				
			}
			
		}
		
	}
	
}

