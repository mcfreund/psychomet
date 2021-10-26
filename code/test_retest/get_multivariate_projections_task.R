#!/usr/bin/env Rscript


## --------------------------------------------------------------------------------------------------------------------
## Script: get_multivariate_projections_session.R 
## Author: Mike Freund (mfreundc@gmail.com)
## Input: 
##  - do_network (boolean)
##  - shrinkage_factor (boolean)
## Output:
##  - ./out/test_retest/*.RDS file of projections
##  ...
## Notes:
##  uses penalized linear discriminant analysis to estimate multivariate hi-lo contrast vector, and project trial-level 
##  data onto vector.
## --------------------------------------------------------------------------------------------------------------------



## setup ----



## packages and sourced variables


library(colorout)
library(here)
library(dplyr)
library(tidyr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mikeutils)
library(ggplot2)
library(purrr)
library(gifti)
# library(corpcor)
# library(klaR)
library(mda)

source(here("code", "_constants.R"))
source(here("code", "_atlases.R"))
source(here("code", "_funs.R"))


## command args

args <- commandArgs(trailingOnly = TRUE)
do_network <- as.integer(args[1]) == 1
shrinkage_factor <- as.integer(args[2])

if (interactive()) {
  do_network <- TRUE
  shrinkage_factor <- 1000
}

message(paste0(do_network, " ", shrinkage_factor))


## other constants

subjs <- subjs_wave12
do_waves <- 1:2
sessions <- "baseline"
tasks <- c("Axcpt", "Cuedts", "Stroop")
train_tasks <- c("Axcpt", "Cuedts")
test_task <- "Stroop"

n_resamples <- 100
n_cores <- 5
fun_resample <- min

vertex_cols <- paste0("V", 1:n_vert)

waves <- waves[do_waves]

if (do_network) {
  rois <- split(key_schaefer$parcel, get_network(key_schaefer$parcel))
  fname_rois <- "_schaefer-7network"
} else {
  rois <- split(key_schaefer$parcel, key_schaefer$parcel)
  fname_rois <- "_schaefer-400parcel"
}


## read trial data

behav <- fread(here::here("in", "behav", "behavior-and-events_wave12_alltasks.csv"))
behav <- behav[task %in% tasks, c("subj", "wave", "task", "session", "trialtype", "hilo_all", "trialnum")]
behav <- rename(behav, hilo = hilo_all)
behav[hilo == ""]$hilo <- NA

## read trial data

behav <- fread(here::here("in", "behav", "behavior-and-events_wave12_alltasks.csv"))
behav <- behav[task %in% tasks, c("subj", "wave", "task", "session", "trialtype", "hilo_all", "trialnum")]
behav <- rename(behav, hilo = hilo_all)
behav[hilo == ""]$hilo <- NA
behav <- behav[task %in% tasks & session %in% sessions & wave %in% waves]
behav$hilo <- factor(behav$hilo, levels = c("lo", "hi"))  ## prep for regression
contrasts(behav$hilo) <- matrix(c(-0.5, 0.5), nrow = 2, dimnames = list(c("lo", "hi"), "hi"))




## execute ----


l <- enlist(combo_paste(sessions, waves, subjs))
if (interactive()) {
  subj_i <- 3
  wave_i <- 1
  session_i <- 1
  task_i <- 1
  roi_i <- which(rois == "LH_SalVentAttn_FrOperIns_3")
}


for (subj_i in seq_along(subjs)) {

  for (wave_i in seq_along(waves)) {

    for (session_i in seq_along(sessions)) {




      subj_val <- subjs[subj_i]
      wave_val <- waves[wave_i]
      session_val <- sessions[session_i]

      ## read trial-wise ("glm-free") fmri coefs:

      res <- read_results(
        waves = wave_val, tasks = tasks, sessions = session_val, subjs = subj_val,
        glmname = "null_2rpm",
        filename_fun = function(...) "errts_trials_target_epoch.RDS",
        read_fun = readRDS
      )

      fmri <- lapply(res, function(x) as.data.table(x))
      names(fmri) <- tasks





      ## prepare for model:

      data_clean <- enlist(tasks)
      for (task_i in seq_along(fmri)) {


        task_val <- tasks[task_i]

        fmri_val <- fmri[[task_val]]
        fmri_val <- fmri_val[, trialnum := 1:.N]  ## make trial number column
        behav_val <- behav[subj == subj_val & wave == wave_val & session == session_val & task == task_val]


        d <- merge(fmri_val, behav_val, by = "trialnum")  ## bind to trial info

        ## remove missing/censored trials:

        is_ok_trial <- !is.na(rowSums(d[, ..vertex_cols]))
        # print(noquote(paste0("removing ", sum(!is_ok_trial), " trials due to missing values")))
        d <- d[is_ok_trial, ]

        d_hilo <- d[!is.na(d$hilo), ]  ## remove non hilo trials (no-go in axcpt)


        ## separate into two dts (corresponding rows)

        d_fmri <- d_hilo[, ..vertex_cols]
        d_trial <- d_hilo[, -..vertex_cols]

        d_trial$is_run1 <- d_trial$trialnum < (n_trialspr[paste0(task_val, "_", session_val)] + 1)


        ## regress nuisance variance:

        hilo <- model.matrix(~ hilo, d_trial)
        X <- cbind(hilo * d_trial$is_run1, hilo * !d_trial$is_run1)

        mu <- coef(.lm.fit(x = X, y = as.matrix(d_fmri)))
        # image(mu)
        mu_bar  <- X[, c(1, 3)] %*% mu[c(1, 3), ]  ## mean of hi/lo means per run
        d_fmri_c <- as.matrix(d_fmri) - mu_bar
        # a <- d_fmri[, mean(V1), by = interaction(is_run1, hilo[, 2] > 0)]$V1  ## for checking
        # all.equal(c(mean(a[1:2]), mean(a[3:4])), unique(mu_bar[, 1]))


        data_clean[[task_i]] <- cbind(d_fmri_c, d_trial)  ## FMRI DATA FIRST, so can easily index by integer later

      }





      ## extract ROI and fit model:

      projs <- matrix(NA, nrow = nrow(d_trial), ncol = length(rois), dimnames = list(trial = NULL, roi = names(rois)))
      projs <- as.data.table(projs)


      for (roi_i in seq_along(rois)) {

        roi_val <- names(rois)[roi_i]


        ## extract:

        which_parcels <- match(rois[[roi_i]], key_schaefer$parcel)  ## works with both network and parcel level
        is_roi <- schaefer10k %in% which_parcels

        d_roi <- lapply(data_clean, function(x) x[, ..is_roi])
        trialinfo <- lapply(data_clean, function(x) x[, c("hilo", "is_run1")])

        ## remove bad verts:

        good_vertices_task <- lapply(d_roi, function(x) which(!is_equal(x[, lapply(.SD, var)], 0)))
        good_vertices <- Reduce(intersect, good_vertices_task)
        d_roi_good <- lapply(d_roi, function(x) x[, ..good_vertices])



        ## fit:

        ## training obs:

        d_test <- d_roi[[test_task]]
        info_test <- trialinfo[[test_task]]

        d_train <- d_roi[train_tasks]
        info_train <- trialinfo[train_tasks]
        # info_train <- lapply(info_train, function(x) x[, ind := 1:.N])
        info_train <- rbindlist(info_train, idcol = "task")
        info_train[, ind := 1:.N]
        d_train <- rbindlist(d_train)

        ## resample

        groups_list <- split(info_train$ind, interaction(info_train[, c("task", "hilo", "is_run1")]))
        resample_to <- Reduce(fun_resample, lapply(groups_list, length))
        inds <- replicate(n_resamples, unlist(lapply(groups_list, sample, size = resample_to)), simplify = FALSE)

        x <- d_train
        y <- info_train$hilo
        x_test <- as.matrix(d_test)

        fits <- mclapply(
          inds,
          function(.inds) {
            fda(y[.inds] ~ as.matrix(x[.inds]), method = gen.ridge, lambda = shrinkage_factor)
          },
          mc.cores = n_cores
        )

        projs_i <- vapply(
          fits,
          function(.x, .newdata) predict(.x, newdata = .newdata, type = "variates"),
          numeric(nrow(x_test)),
          .newdata = x_test
        )
        # proj_bar <- plogis(rowMeans(qlogis(projs_i)))
        proj_bar <- rowMeans(projs_i)

        projs[[roi_val]] <- proj_bar

      }

      projs$trialnum <- d_trial$trialnum
      projs$hilo <- d_trial$hilo



      ## store predictions:
      nm <- paste0(session_val, "_", wave_val, "_", subj_val)
      l[[nm]] <- projs

      print(paste0("done ", subj_i, " ", wave_i, " ", session_i))



    }  ## session
  }  ## wave
}  ## subj



dat <- rbindlist(l, idcol = "wave_subj")
dat <- separate(dat, wave_subj, c("session", "wave", "subj"))
dat <- dat %>% melt(id.vars = c("session", "wave", "subj", "trialnum", "hilo"), value.name = "p")

saveRDS(
  dat,
  here(
    "out", "test_retest",
    paste0("stroop_projections_rda_lambda-", shrinkage_factor, "_cross-task", fname_rois, "-variates.RDS")
  )
)


