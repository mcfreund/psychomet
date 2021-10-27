
library(colorout)
library(here)
library(reticulate)
library(tidyr)
library(dplyr)
library(magrittr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mikeutils)
library(progress)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(viridis)
library(purrr)

theme_set(theme_half_open())

source(here("code", "_constants.R"))
source(here("code", "_atlases.R"))
source(here("code", "_funs.R"))
source(here("code", "_read_behav.R"))

hi <- c("BX", "InCon", "RN")
lo <- c("BY", "Con", "NN")


subjs <- subjs_wave12


## read data

d <- read_results(
  c("wave1", "wave2"), tasks, "baseline", subjs,
  glmname = "null_2rpm", 
  filename_fun = function(...) "means-trials_schaefer-parcel400_resid-errts.csv",
  read_fun = fread
)
d <- rbindlist(d, idcol = "id")
d <- d[, trialnum := 1:.N, by = "id"]  ## make trial number column
d <- separate(d, id, c("wave", "task", "session", "subj"))



## check parcel means:

a <- readRDS(here("out", "icc", "means-trials_schaefer-parcel_resid-errts.RDS"))
subjs_common <- intersect(a$subj, subjs)
a <- a[subj %in% subjs_common]
w <- melt(d[wave == "wave1" & subj %in% subjs_common], id.vars = c("trialnum", "subj", "task", "wave", "session"))
w[, c("wave", "session") := NULL]
w <- rename(w, trial = trialnum, roi = variable)
a <- arrange(a, subj, task, roi, run, trial)
a[, trial := 1:.N, by = c("subj", "roi", "task")]
a[, c("run") := NULL]
m <- merge(a, w, by = c("subj", "roi", "task", "trial"))

m$values_are_equal <- is_equal(m$value.x, m$value.y)


m_sum <- m[, .(
  error = sqrt(mean((value.x - value.y)^2, na.rm = TRUE)),
  r = cor(value.x, value.y, use = "complete")
), by = c("subj", "roi", "task")]

m_sum %>%
  ggplot(aes(subj, r)) +
  geom_boxplot() +
  facet_grid(cols = vars(task)) +
  coord_flip()

## 448347 is quite weird: Axcpt, Cuedts, and Stroop
## 130518 is too: Stern, Stroop
## 203418: stroop


## check residuals:

subj_val <- "448347"
task_val <- "Stern"
hemi_val <- "L"


named_vector <- function(type, nms) setNames(vector(type, length(nms)), nms)

r <- named_vector("numeric", combo_paste(tasks, subjs_common, c("L", "R")))
for (task_val in tasks) {
  
  for (subj_val in subjs_common) {
    
    for (hemi_val in c("L", "R")) {
      
      gii_ub55_lh <- rbind(
        read_gifti2matrix(
          file.path(
            "/data/nil-external/ccp/freund/ub55/out/glms", subj_val, "RESULTS", task_val, "baseline_null_1",
            paste0("errts_", subj_val, "_1_", hemi_val, "_REML.func.gii")
          )
        ),
        read_gifti2matrix(
          file.path(
            "/data/nil-external/ccp/freund/ub55/out/glms", subj_val, "RESULTS", task_val, "baseline_null_2",
            paste0("errts_", subj_val, "_2_", hemi_val, "_REML.func.gii")
          )
        )
      )
      
      gii_psych_lh <- rbind(
        read_gifti2matrix(
          file.path(
            "/data/nil-external/ccp/freund/psychomet/out/glms", subj_val, "RESULTS", task_val, "baseline_null_2rpm_wave1",
            paste0("errts_", subj_val, "_", hemi_val, "_REML.func.gii")
          )
        )
      )
      
      giidat <- data.frame(ub = c(gii_ub55_lh), ps = c(gii_psych_lh))
      
      r[paste0(task_val, "_", subj_val, "_", hemi_val)] <- cor(giidat)[1, 2]
      
      
    }
    
  }
  
}

giidat %>%
  ggplot(aes(ub, ps)) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline()


## check vertex means:

subj_val <- "448347"
task_val <- "Stroop"

r_trial <- named_vector("numeric", combo_paste(tasks, subjs_common))
for (task_val in tasks) {
  
  for (subj_val in subjs_common) {
    
    ## psychomet:
    
    dir_glm <- here("out", "glms", subj_val, "RESULTS", task_val, "baseline_null_2rpm_wave1")
    fname_resids <- here(dir_glm, "errts_trials_target_epoch.RDS")
    resids <- readRDS(fname_resids)
    
    ## ub55:
    
    dir_glm_ub <-
      here("..", "ub55", "out", "glms", subj_val, "RESULTS", task_val, paste0("baseline_null_", 1:2))
    fname_resids_ub <- here(dir_glm_ub, "errts_trials_target_epoch.RDS")  ## both runs
    resids_ub <- lapply(fname_resids_ub, readRDS)
    resids_ub <- abind(resids_ub, along = 1)
    
    trialdat <- data.frame(ub = c(resids), ps = c(resids_ub))
    r_trial[paste0(task_val, "_", subj_val)] <- cor(trialdat, use = "complete")[1, 2]
    
    
  }
  
}

trial <- 206
data.frame(ub = resids[trial, ], ps = resids_ub[trial, ]) %>%
  ggplot(aes(ub, ps)) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline()



trialdat %>%
  ggplot(aes(ub, ps)) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline()

