library(dplyr)
library(data.table)
library(here)
library(purrr)

source(here("code", "_constants.R"))

## read and combine:

behav_wave1 <- list(
  Axcpt  = fread(here("in", "behav", "dmcc2_behavior-and-events_axcpt_2021-09-13.csv")),
  Cuedts = fread(here("in", "behav", "dmcc2_behavior-and-events_cuedts_2021-09-13.csv")),
  Stern  = fread(here("in", "behav", "dmcc2_behavior-and-events_sternberg_2021-09-13.csv")),
  Stroop = fread(here("in", "behav", "dmcc2_behavior-and-events_stroop_2021-09-13.csv"))
)


behav_wave2 <- list(
  Axcpt  = fread(here("in", "behav", "dmcc3_behavior-and-events_axcpt_2021-09-23.csv")),
  Cuedts = fread(here("in", "behav", "dmcc3_behavior-and-events_cuedts_2021-09-23.csv")),
  Stern  = fread(here("in", "behav", "dmcc3_behavior-and-events_sternberg_2021-09-23.csv")),
  Stroop = fread(here("in", "behav", "dmcc3_behavior-and-events_stroop_2021-09-23.csv"))
)

behav <- map2(behav_wave1, behav_wave2, ~(rbind(.x, .y, idcol = "wave")))



## subset subjs

behav <- lapply(behav, function(x) x[subj %in% subjs_wave12])



## format task-unique vals/cols


## Axcpt:

behav$Axcpt <- behav$Axcpt %>% 
  select(-acc) %>% 
  rename(
    trialtype = trial.type,
    rt = target.rt, acc = target.acc, cresp = target.cresp, resp = target.resp
    )


## Cuedts:

behav$Cuedts[trial.type == "i"]$trial.type <- "InCon"
behav$Cuedts[trial.type == "c"]$trial.type <- "Con"

behav$Cuedts[target.color.orig == "green"]$incentive <- "Inc"
behav$Cuedts[target.color.orig == "black"]$incentive <- "NoInc"

behav$Cuedts$trialtype <- paste0(behav$Cuedts$trial.type, behav$Cuedts$incentive)


## Stern:

behav$Stern$resp <- ifelse(behav$Stern$resp == 1, 2, ifelse(behav$Stern$resp == 2, 1, NA))
behav$Stern$cresp <- ifelse(behav$Stern$cresp == 1, 2, ifelse(behav$Stern$cresp == 2, 1, NA))

behav$Stern$trialtype <- paste0(ifelse(behav$Stern$load == 5, "LL5", "not5"), behav$Stern$trial.type)


## Stroop:

behav$Stroop[trial.type == "i"]$trial.type <- "InCon"
behav$Stroop[trial.type == "c"]$trial.type <- "Con"

behav$Stroop[pc %in% c("mi", "mc")]$pc <- "bias"
behav$Stroop[pc == "pc50"]$pc <- "PC50"

behav$Stroop$trialtype <- paste0(behav$Stroop$pc, behav$Stroop$trial.type)

behav$Stroop[
  wave == 2 & session == "bas" & subj == "DMCC6705371" & run == 1 & trial.num == 59
  ]$response.final <- "pink"
  
behav$Stroop$resp <- 
  ifelse(
    behav$Stroop$acc.final %in% c("unintelligible", "no.response"), 
    behav$Stroop$acc.final, 
    behav$Stroop$response.final
    )

behav$Stroop <- behav$Stroop %>% rename(cresp = correctcolor)



## extract cols and bind


cols <- c("subj", "wave", "session", "run", "trial.num",  "trialtype", "cresp", "resp", "acc", "rt")
behav <- lapply(behav, function(x) x[, ..cols])

behav <- rbindlist(behav, idcol = "task")


## format vals:

behav[session == "bas"]$session <- "baseline"
behav[session == "pro"]$session <- "proactive"
behav[session == "rea"]$session <- "reactive"

behav <- behav[order(task, subj, wave, session, run, trial.num)]
behav[, trialnum := 1:.N, by = c("task", "subj", "wave", "session")]  ## not over run
behav[, c("run", "trial.num") := NULL]

behav[, wave := paste0("wave", wave)]


## make hilo column:

behav$hilo <- ifelse(behav$trialtype %in% hi, "hi", ifelse(behav$trialtype %in% lo, "lo", NA))

behav$hilo_all <- NA
behav$hilo_all[grepl("^BX$|^InCon|RN$|InCon$", behav$trialtype)] <- "hi"
behav$hilo_all[grepl("^BY$|^Con|NN$|biasCon|PC50Con$", behav$trialtype)] <- "lo"

# table(behav$hilo, behav$trialtype)
# table(behav$hilo_all, behav$trialtype)


## write:

# map2(behav, here("in", "behav", paste0("behavior-and-events_wave12_", tasks, ".csv")), fwrite)
fwrite(behav, here("in", "behav", paste0("behavior-and-events_wave12_alltasks.csv")))
# table(behav$hilo_all, behav$trialtype, behav$task)