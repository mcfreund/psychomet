library(dplyr)
library(data.table)
library(here)
library(purrr)

source(here("code", "_constants.R"))

behav_wave1 <- list(
  Axcpt  = fread(here("in", "behav", "dmcc2_behavior-and-events_axcpt_2021-09-13.csv")) %>% rename(resp = target.resp),
  Cuedts = fread(here("in", "behav", "dmcc2_behavior-and-events_cuedts_2021-09-13.csv")),
  Stern  = fread(here("in", "behav", "dmcc2_behavior-and-events_sternberg_2021-09-13.csv")),
  Stroop = fread(here("in", "behav", "dmcc2_behavior-and-events_stroop_2021-09-13.csv"))
)

behav_wave1$Stern$resp <- ifelse(behav_wave1$Stern$resp == 1, 2, ifelse(behav_wave1$Stern$resp == 2, 1, NA))
behav_wave1$Stern$cresp <- ifelse(behav_wave1$Stern$cresp == 1, 2, ifelse(behav_wave1$Stern$cresp == 2, 1, NA))

behav_wave1 <- lapply(behav_wave1, function(x) x[subj %in% subjs_wave12])


behav_wave2 <- list(
  Axcpt  = fread(here("in", "behav", "dmcc3_behavior-and-events_axcpt_2021-09-23.csv")) %>% rename(resp = target.resp),
  Cuedts = fread(here("in", "behav", "dmcc3_behavior-and-events_cuedts_2021-09-23.csv")),
  Stern  = fread(here("in", "behav", "dmcc3_behavior-and-events_sternberg_2021-09-23.csv")),
  Stroop = fread(here("in", "behav", "dmcc3_behavior-and-events_stroop_2021-09-23.csv"))
)

behav_wave2$Stern$resp <- ifelse(behav_wave2$Stern$resp == 1, 2, ifelse(behav_wave2$Stern$resp == 2, 1, NA))
behav_wave2$Stern$cresp <- ifelse(behav_wave2$Stern$cresp == 1, 2, ifelse(behav_wave2$Stern$cresp == 2, 1, NA))

behav_wave2 <- lapply(behav_wave2, function(x) x[subj %in% subjs_wave12])


behav <- map2(behav_wave1, behav_wave2, ~(rbind(.x, .y, idcol = "wave")))


## format cols:

behav$Cuedts[trial.type == "i"]$trial.type <- "InCon"
behav$Cuedts[trial.type == "c"]$trial.type <- "Con"

behav$Cuedts[target.color.orig == "green"]$incentive <- "Inc"
behav$Cuedts[target.color.orig == "black"]$incentive <- "NoInc"

behav$Cuedts$trialtype <- paste0(behav$Cuedts$trial.type, behav$Cuedts$incentive)


behav$Stroop[trial.type == "i"]$trial.type <- "InCon"
behav$Stroop[trial.type == "c"]$trial.type <- "Con"

behav$Stroop[pc %in% c("mi", "mc")]$pc <- "bias"
behav$Stroop[pc == "pc50"]$pc <- "PC50"

behav$Stroop$trialtype <- paste0(behav$Stroop$pc, behav$Stroop$trial.type)


behav$Stern$trialtype <- paste0(ifelse(behav$Stern$load == 5, "LL5", "not5"), behav$Stern$trial.type)

behav$Axcpt <- rename(behav$Axcpt, trialtype = trial.type)



## write:

map2(behav, here("in", "behav", paste0("behavior-and-events_wave12_", tasks, ".csv")), fwrite)
