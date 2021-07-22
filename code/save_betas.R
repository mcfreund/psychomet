
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))


## input: subjs, task, glmname
## output: RDS files
# https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/


dirs <- expand.grid(subj = subjs, task = tasks, session = "baseline", stringsAsFactors = FALSE)
glminfo <- data.frame(
  task = "Cuedts",
  # task = c("Axcpt", "Cuedts", "Stern", "Stroop"),#, "Stroop"),
  name.glm = c(
    "baseline_cueletnum_EVENTS_censored_shifted"
    # "baseline_CongruencySwitch_EVENTS_censored_shifted",
    # "baseline_ListLength_EVENTS_censored_shifted",
    # "baseline_Congruency_EVENTS_censored_shifted"
    # "baseline_fix-item_EVENTS_censored"
  )
)
glminfo <- as.data.table(glminfo)



for (glm.i in seq_len(nrow(glminfo))) {
  
  betas.i <- read_betas(subjs, glminfo[glm.i]$task, glminfo[glm.i]$name.glm, dir.analysis)
  saveRDS(betas.i, here("out", "glms", paste0("betas_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")))
  
}


for (glm.i in seq_len(nrow(glminfo))) {
  
  sustained.i <- read_coefs(
    subjs, 
    glminfo[glm.i]$task, 
    glminfo[glm.i]$name.glm, 
    dir.analysis, 
    .getregs = function(x) grepl("block#0_Coef", x)
  )
  
  saveRDS(
    sustained.i, 
    here("out", "glms", paste0("sustained_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS"))
  )
  
}

