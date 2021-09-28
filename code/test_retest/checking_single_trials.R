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
source(here("..", "ub55", "code", "_funs.R"))  ## for read_betas_dmcc()
source(here("code", "_funs.R"))
source(here("code", "_read_behav_wave12.R"))



## input vars ----

atlas <- "schaefer"
subjs <- subjs_wave12
do_waves <- c(1, 2)
hi <- c(Axcpt = "BX", Cuedts = "InConInc", Stern = "LL5RN", Stroop = "biasInCon")
lo <- c(Axcpt = "BY", Cuedts = "ConInc", Stern = "LL5NN", Stroop = "biasCon")


## wrangle behavioral data:

behav_wave12$Axcpt <- rename(behav_wave12$Axcpt, rt = target.rt)
cols <- c("wave", "run", "session", "subj", "trialtype", "rt", "acc", "trial.num")
behav_wave12 <- lapply(behav_wave12, function(x) x[, ..cols])
behav_wave12 <- rbindlist(behav_wave12, idcol = "task")
behav_wave12 <- behav_wave12[order(wave, task, session, subj, run, trial.num), ]
behav_wave12[, trialnum := 1:.N, by = c("wave", "task", "subj", "session")]
behav_wave12[, c("run", "trial.num") := NULL]
behav_wave12[, wave := paste0("wave", wave)]

behav_wave12[session == "bas"]$session <- "baseline"
behav_wave12[session == "pro"]$session <- "proactive"
behav_wave12[session == "rea"]$session <- "reactive"



## execute ----

waves <- waves[do_waves]


# wave_i = 2
# task_i = 3
# session_i = 2
# task_val <- tasks[task_i]
# wave_val <- waves[wave_i]
# session_val <- sessions[session_i]


res <- enlist(combo_paste(waves, tasks, sessions))
for (wave_i in seq_along(waves)) {
  
  for (task_i in seq_along(tasks)) {
    
    for (session_i in seq_along(sessions)) {
     
      
      task_val <- tasks[task_i]
      wave_val <- waves[wave_i]
      session_val <- sessions[session_i]
      trialtype_val <- c(hi = hi[[task_val]], lo = lo[[task_val]])
      
      
      ## read trial-level residuals, condition-level betas
      
      resid <- read_results(
        wave_val, task_val, session_val, subjs,
        glmname = "null_2rpm", filename = "errts_trials_target_epoch.RDS",
        readRDS
      )
      names(resid) <- gsub(paste0(wave_val, "_", task_val, "_", session_val, "_"), "", names(resid))
      
      betas <- read_betas_dmcc(
        .subjs = subjs, 
        .task = task_val, 
        .glm = paste0(session_val, "_", name_glms_dmcc[task_val]), 
        .dir = file.path("/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS", wavedir_image[wave_val], "fMRIPrep_AFNI_ANALYSIS")
      )
      
      
      ## extract core32
      
      is_core32 <- schaefer10k %in% core32
      resid <- lapply(resid, function(x) x[, is_core32])
      betas <- betas[is_core32, , , ]
      
      
      ## aggregate over trials by trial-type
      
      l <- 
        behav_wave12[wave == wave_val & task == task_val & session == session_val] %>% 
        split(.$subj)
      
      A <- 
        lapply(
          l, 
          function(x) {
            A <- model.matrix(~ 0 + trialtype, x)
            A <- sweep(A, 2, colSums(A), "/")
            colnames(A) <- gsub("trialtype", "", colnames(A))
            A[, trialtype_val] %*% rbind(1, -1)  ## gives hilo contrast
            }
          )
      
      resid <- resid[names(A)]  ## make sure have same order
      
      resid_sum <- map2(A, resid, ~crossprod(.x, .y))
      #resid_sum <- map2(A, resid, ~as.data.table(crossprod(.x, .y), keep.rownames = "trialtype"))
      # resid_sum <- rbindlist(resid_sum, idcol = "subj")
      # resid_sum_d <- melt(resid_sum, id.vars = c("subj", "trialtype"), variable.name = "vertex")
      # resid_sum_d[, trialtype := gsub("trialtype", "", trialtype)]
      resid_sum <- abind(resid_sum, along = 1)
      
      
      # ## average trial-level residuals by trial-type
      # 
      # resid_d <- as.data.table(reshape2::melt(resid, varnames = c("trialnum", "vertex")))
      # resid_d[, subj := gsub(paste0(wave_val, "_", task_val, "_", session_val, "_"), "", L1)]  ## create subj col
      # 
      # l <- behav_wave12[wave == wave_val & task == task_val & session == session_val]
      # resid_d <- merge(resid_d, l, by = c("subj", "trialnum"))
      # 
      # resid_sum <- resid_d[, .(value = mean(value)), by = c("subj", "vertex", "trialtype")]
      
      
      
      ## extract target knots and average betas:
      
      b <- betas[, trialtype_val, target_trs[[task_val]], ]
      
      b <- aperm(b, c(1, 2, 4, 3))  ## put TR on 'outside'
      b <- rowMeans(b, dims = 3)
      b <- b[, trialtype_val["hi"], ] - b[, trialtype_val["lo"], ]
      # b <- as.data.table(reshape2::melt(b))
      # b <- rename(b, trialtype = reg)
      # b <- b[, vertex := paste0("V", vertex)]
      
      
      
      ## bind and correlate:
      resid_sum <- t(resid_sum)
      
      r <- colSums(scale2unit(meancenter(resid_sum)) * scale2unit(meancenter(b)))
      neg_rmse <- -sqrt(colMeans((resid_sum - b) * (resid_sum - b)))

      # d <- merge(resid_sum_d, b, by = c("vertex", "trialtype", "subj"))
      # r <- d[,
      #   .(
      #     r = cor(value.x, value.y),
      #     neg_mse = -sqrt(mean((value.x - value.y)^2))
      #     ),
      #   by = c("trialtype", "subj")
      #   ]
      
      nm <- paste0(wave_val, "_", task_val, "_", session_val)
      res[[nm]] <- data.table(subj = colnames(b), r = r, neg_rmse = neg_rmse)
      
      print(nm)

      
    }
    
  }
  
}


# res <- res[-grep("wave3", names(res))]


r <- rbindlist(res, idcol = "id")
r <- separate(r, id, c("wave", "task", "session"))



## plot ----


r %>%
  ggplot(aes(r, trialtype)) +
  geom_vline(xintercept = 0) +
  geom_boxplot(width = 0.2, fill = "grey40") +
  facet_wrap(vars(task), ncol = 1, scales = "free_y")


r %>%
  ggplot(aes(r, subj)) +
  geom_vline(xintercept = 0) +
  geom_boxplot(width = 0.2, fill = "grey40") +
  facet_wrap(vars(task), ncol = 1)
