## system info

n_core <- parallel::detectCores()

## paths

dir_atlas <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"



# if (nodename == "ccplinux1") {
#   
#   # dir_atlas <- "/data/nil-external/ccp/freund/atlases"
#   # dir_schaefer <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
#   # dir_mmp <-      "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
#   dir_atlas <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
#   
#   
# } else if (nodename == "CCP-FREUND") {
#   ## mike freund's (i.e., ccp's) thinkpad
#   ## reliant on box drive
#   ## assumes box drive location at ./Users/mcf/Box
#   
#   # dir_atlas <- "C:/local/atlases"
#   # dir_schaefer <- dir_atlas
#   # dir_mmp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
#   
# } else if (nodename == "PUTER") {
#   
#   # dir_atlas <- "C:/Users/mcf/Documents/atlases"
#   # dir_schaefer <- file.path(dir_atlas, "ATLASES")
#   # dir_mmp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
#   
# }

## image, design, analysis info

n_vert <- 20484  ## surface hcp mesh
# n_trs <- c(
#   Axcpt   = 1220,
#   # Axcpt_proactive  = 1220,
#   # Axcpt_reactive   = 1220,
#   Cuedts  = 1300,
#   # Cuedts_proactive = 1300,
#   # Cuedts_reactive  = 1300,
#   Stern   = 1200,
#   # Stern_proactive  = 1200,
#   # Stern_reactive   = 1200,
#   Stroop  = 1080
#   # Stroop_proactive = 1080,
#   # Stroop_reactive  = 1180
# )  ## FOR BASELINE/PROACTIVE SESSION ONLY
n_trs <- c(
  Axcpt_baseline   = 1220,
  Axcpt_proactive  = 1220,
  Axcpt_reactive   = 1220,
  Cuedts_baseline  = 1300,
  Cuedts_proactive = 1300,
  Cuedts_reactive  = 1300,
  Stern_baseline   = 1200,
  Stern_proactive  = 1200,
  Stern_reactive   = 1200,
  Stroop_baseline  = 1080,
  Stroop_proactive = 1080,
  Stroop_reactive  = 1180
)  ## BOTH RUNS INCLUDED
n_trialspr <- c(
  Axcpt_baseline = 72, 
  Axcpt_proactive = 72, 
  Axcpt_reactive = 72, 
  Cuedts_baseline = 54, 
  Cuedts_proactive = 54, 
  Cuedts_reactive = 54, 
  Stern_baseline = 45, 
  Stern_proactive = 45, 
  Stern_reactive = 45, 
  Stroop_baseline = 108,
  Stroop_proactive = 108,
  Stroop_reactive = 120
  )  ## number of trials (events) per subj*run


tasks <- c("Axcpt", "Cuedts", "Stern", "Stroop")
taskruns <- sort(mikeutils::combo_paste(tasks, c("run1", "run2")))
sessions <- c("baseline", "proactive", "reactive")
waves <- c("wave1", "wave2", "wave3")
wavedir_image <- c("HCP_SUBJECTS_BACKUPS", "DMCC_Phase3")
wavedir_evts <- c("DMCC2", "DMCC3", "DMCC4")

subjs_ub55 <- data.table::fread(here::here("..", "ub55", "in", "ub55_subjects.txt"))[[1]]
subjs_wave12 <- fread(here("in", "subjects_wave12_complete_2021-09-01.txt"))$V1
subjs_retest_tang <- c(
  "178243",
  "130518",
  "393550",
  "178647",
  "877168",
  "DMCC5009144",
  "DMCC6627478",
  "DMCC9478705",
  "DMCC6671683",
  "DMCC2609759",
  "DMCC6418065",
  "DMCC8033964",
  "DMCC6721369",
  "DMCC9441378",
  "DMCC6705371"
)

## TRs of interest

## old:
# target_trs <- list(
#   Axcpt = 7:9,
#   Cuedts = 9:10,
#   Stern = 11:12,
#   Stroop = 2:4
# )
# cue_trs <- list(
#   Axcpt = 8:9,  ### ????
#   Cuedts = 9:10, ### ????
#   Stern = 7:8,
#   Stroop = 2:4
# )


## from jo:
# Axcpt: Cues, BX high, BY low. 17, 8:10
# Cuedts: CongruencyIncentive, InConNoInc high, ConNoInc low. 19, 9:11
# Stern: ListLength, LL5RN high, LL5NN low. 21, 12:14.
# Stroop: Congruency, biasInCon high, biasCon low. 13, 3:5

target_trs <- list(
  Axcpt = 8:10,
  Cuedts = 9:11,
  Stern = 12:14,
  Stroop = 3:5
)
