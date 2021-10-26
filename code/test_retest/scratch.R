






# 
# 
# res <- as.matrix(projs[, -"trialnum"])
# 
# 
# 
# a <- plogis(
#   colMeans(qlogis(res[info_test$hilo == "hi", ])) - 
#     colMeans(qlogis(res[info_test$hilo == "lo", ]))
#   )
# 
# ar <- 
#   colMeans(res[info_test$hilo == "hi", ]) - 
#     colMeans(res[info_test$hilo == "lo", ])
# 
# 
# 
# head(sort(-a), 10)
# head(sort(ar, TRUE), 10)
# 
# 
# cla <- res > 0.5
# ac <- colMeans(cla == (info_test$hilo == "hi"))
# 
# head(sort(ac, TRUE), 10)
# head(sort(ac, FALSE), 10)
# 
# qlogis(as.matrix(projs[, -"trialnum"]))
# as.numeric(info_test$hilo) - 1.5
#   
# 
#   
#   
#   
#   
#   
#   
#   
#   
#   wave_sessions <- names(fmri)  ## iterate over
#   
#   m <- enlist(wave_sessions)
#   
#   for (wave_session_i in seq_along(wave_sessions)) {
#     
#     ## iterator values:
#     
#     # wave_session_val <- wave_sessions[wave_session_i]
#     # wave_session_split <- unlist(strsplit(wave_session_val, "_"))
#     # wave_val <- wave_session_split[1]
#     # session_val <- wave_session_split[3]
#     
#     ## extract:
#     
#     fmri_val <- fmri[[wave_session_i]]
#     fmri_val <- fmri_val[, trialnum := 1:.N]  ## make trial number column
#     behav_val <- behav[subj == subj_val & wave == wave_val & session == session_val]
#     
#       
#     d <- merge(fmri_val, behav_val, by = "trialnum")  ## bind to trial info
#     
#     ## remove missing/censored trials:
#     
#     vertex_cols <- paste0("V", 1:n_vert)
#     # print(noquote(paste0("removing ", sum(is.na(rowSums(d[, ..vertex_cols]))), " trials due to missing values")))
#     d <- d[!is.na(rowSums(d[, ..vertex_cols])), ]
#     
#     
#     ## remove non hilo trials
#     
#     d_hilo <- d[!is.na(d$hilo), ]
#     
#     
#     ## separate into two dts (corresponding rows)
#     
#     d_fmri <- d_hilo[, ..vertex_cols]
#     d_trial <- d_hilo[, -..vertex_cols]
#     
#     projs <- matrix(NA, nrow = nrow(d_trial), ncol = length(rois), dimnames = list(trial = NULL, roi = rois))
#     projs <- as.data.table(projs)
#     # options(warn=2)
#     # options(warn=1)
#     for (roi_i in seq_along(rois)) {
#       
#       roi_val <- rois[roi_i]
#       
#       
#       ## extract:
#       
#       which_parcels <- match(rois[[roi_i]], key_schaefer$parcel)  ## works with both network and parcel level
#       is_roi <- schaefer10k %in% which_parcels
#       d_roi <- d_fmri[, ..is_roi]
#       
#       
#       ## remove bad verts:
#       
#       is_good_vertex <- !is_equal(d_roi[, lapply(.SD, var)], 0)
#       d_roi <- d_roi[, ..is_good_vertex]
#       # d_roi <- as.data.table(t(scale(t(d_roi))))
#       
#       
#       ## get run inds:
#       
#       is_run1 <- d_trial$trialnum < (n_trialspr[paste0("Stroop_", session_val)] + 1)
#       
#       
#       ## fit lda:
#       
#       x_run1 <- d_roi[is_run1]
#       x_run2 <- d_roi[!is_run1]
#       
#       y_run1 <- d_trial$hilo[is_run1]
#       y_run2 <- d_trial$hilo[!is_run1]
# 
#       train_run1 <- rda(x = x_run1, grouping = y_run1, gamma = 0.5, lambda = 1, prior = 1)
#       train_run2 <- rda(x = x_run2, grouping = y_run2, gamma = 0.5, lambda = 1, prior = 1)
#       
#       proj_run1 <- predict(train_run2, newdata = x_run1)
#       proj_run2 <- predict(train_run1, newdata = x_run2)
#       
#       # mean(proj_run1$class == d_trial$hilo[is_run1])
#       # mean(proj_run2$class == d_trial$hilo[!is_run1])
#       
#       projs[[roi_val]] <- c(proj_run1$posterior[, "hi"], proj_run2$posterior[, "hi"])
# 
#     }
#     
#     projs$trialnum <- d_trial$trialnum
#     
#     m[[wave_session_val]] <- projs
#     
#     
#   }
#   
#   
#   l[[subj_val]] <- rbindlist(m, idcol = "wave_session")
#   
#   print(subj_val)
#   
# }
# # 
# # 
# 
# 
# ## save ----
# # RH_Limbic_OFC_3
# 
# allres <- rbindlist(l, idcol = "subj")
# # saveRDS(allres, here("out", "test_retest", "posterior_cross-run-cval_rlda_Stroop.RDS"))
# 
# 
# # 1: wave1 Stroop baseline      178647      114             LH_Vis_12     NaN
# 
# ## examine ----
# 
# allres <- separate(allres[, -"subj"], wave_session, c("wave", "task", "session", "subj"))
# allres_l <- melt(
#   allres, 
#   id.vars = c("wave", "task", "session", "subj", "trialnum"), 
#   variable.name = "roi", value = "p_incon"
#   )
# sum(is.nan(allres_l$p_incon))
# allres_l[is.nan(p_incon)]
# allres_l <- allres_l[!is.nan(p_incon)]  ## i assume issues with inversion
# d <- left_join(allres_l, behav, by = c("wave", "session", "subj", "trialnum"))
# 
# d$classification <- sign(d$p_incon - 0.5)
# d$true <- sign(as.numeric(d$hilo) - 1.5)
# d$correct <- d$classification*d$true > 0
# 
# 
# d_sum <- d[, .(acc = mean(correct)), by = c("wave", "session", "subj", "roi")]
# group_sum <- d_sum %>%
#   group_by(wave, session, roi) %>%
#   summarize(acc = mean(acc))
# 
# group_sum %>%
# 
#   filter(roi %in% rois[core32]) %>%
#   arrange(-acc) %>%
#   as.data.table %>%
#   head(100)
#   ggplot(aes(fill = roi %in% rois[core32])) +
#   geom_histogram(aes(x = acc))
# 
# 
# prob2logit <- function(p) ln(p / (1 - p))
# 
# d_sum <- d[, .(acc = mean(qlogis(p_incon))), by = c("wave", "session", "subj", "roi")]
# group_sum <- d_sum %>%
#   group_by(wave, session, roi) %>%
#   summarize(acc = mean(acc))
# 
# group_sum %>%
#   
#   # filter(session == "baseline", wave == "wave1") %>%
#   arrange(-acc) %>%
#   as.data.table %>%
#   head(100)



# 
# 
# d[,
#   .(
#     sign(p_incon - 0.5)
#     ),
#   by = c("")
#   ]


# atanh(pchisq(y_proj, df = ncol(Y)))
# crossprod(sign(as.matrix(projs)), sign(X[, 2]))



# 
# allres <- readRDS(here("out", "test_retest", "posterior_cross-run-cval_rlda_Stroop.RDS"))
# allres_l <- melt(allres, id.vars = c("subj", "wave_session", "trialnum"), variable.name = "roi", value = "p_incon")
# allres_l <- separate(allres_l[, -"subj"], wave_session, c("wave", "task", "session", "subj"))
# d <- merge(allres_l, behav, by = c("wave", "session", "subj", "trialnum"))








## scratch ----



## 
## 
## 
##     # X <- model.matrix(~hilo_all, d_trial)
# X <- model.matrix(~hilo, d_trial)
# colnames(X) <- c("b0", "b1")


##       ## fit model:
# 
# # Y <- as.matrix(d_roi)[!is.na(d_trial$hilo), ]
# Y <- as.matrix(d_roi)
# 
# fit <- .lm.fit(X, Y)
# B <- coef(fit)
# E <- resid(fit)
# b0 <- B[1, ]  ## mean of means
# b1 <- B[2, ]  ## hi vs lo difference vector
# 
# 
# ## estimate axis:
# 
# W2 <- cov.shrink(E)
# a <- scale2unit(W2 %*% b1)  ## hi-lo axis
# 
# 
# ## project:
# 
# Y_c <- sweep(Y, 2, b0)  ## demean Y
# y_proj <- Y_c %*% a
# 
# 
# projs[[roi_val]] <- y_proj


## reampling:
## 
# set.seed(0)
# 
# 
# reglda <- function(x, y, shrinkage_factor = 0.4, n_resamples = 1E4, n_cores = 1, quick = FALSE) {
#   # x <- d_roi
#   # y <- as.character(d_trial$hilo)
# 
#   ## downsample conditions:
# 
#   counts <- table(y)
#   min_count <- min(counts)
# 
#   ## inds for resampling:
# 
#   groups_list <- split(seq_along(y), y)
#   inds <- lapply(
#     seq_len(n_resamples),
#     function(.x) unlist(lapply(groups_list, sample, size = min_count), use.names = FALSE)
#   )
# 
#   ## fit models:
# 
#   fits <- mclapply(
#     inds,
#     function(.inds) {
#       rda(x = x[.inds], grouping = y[.inds], gamma = shrinkage_factor, lambda = 1)
#     },
#     mc.cores = n_cores
#   )
# 
#   ## generate predictions:
# 
#   projs_i <- vapply(
#     fits,
#     function(.x, newdata, quick) predict(.x, newdata = newdata, quick = quick)$posterior[, "hi"],
#     numeric(nrow(x)),
#     newdata = x,
#     quick = quick
#   )
# 
# 
#   projs_i
# 
# 
# }
# 
# 
# 
# 
# 
# 
# # projs <- reglda(d_roi, d_trial$hilo, n_cores = n_core - 4, quick = TRUE)
# # proj <- tanh(rowMeans(atanh(projs)))
## 
## 
## 
## 
##       ## fit model:
# 
# Y <- as.matrix(d_roi)[!is.na(d_trial$hilo), ]
# Y <- as.matrix(d_roi)
# 
# fit <- .lm.fit(X, Y)
# B <- coef(fit)
# E <- resid(fit)
# b0 <- B[1, ]  ## mean of means
# b1 <- B[2, ]  ## hi vs lo difference vector
# 
# 
# ## estimate axis:
# 
# W2 <- cov.shrink(E)
# a <- scale2unit(W2 %*% b1)  ## hi-lo axis
# 
# 
# ## project:
# 
# Y_c <- sweep(Y, 2, b0)  ## demean Y
# y_proj <- Y_c %*% a
# 
# 
# projs[[roi_val]] <- y_proj
# 

## functions
# 
# cvrda <- function(x, y, is_train, is_test, shrinkage_factor, lambda = 1, prior = 1, ...) {
#   
#   x_train <- x[is_train]
#   y_train <- y[is_train]
#   x_test <- x[is_test]
#   
#   if (shrinkage_factor == "estimate") {
#     
#     fit_lm <- .lm.fit(model.matrix(~ y_train), as.matrix(x_train))
#     E <- resid(fit_lm)
#     shrinkage_factor <- estimate.lambda(E)
#     
#   }
#   
#   
#   train <- rda(
#     x = x_train, 
#     grouping = y_train, 
#     gamma = shrinkage_factor,
#     lambda = lambda,
#     prior = prior,
#     ...
#   )
#   
#   test <- predict(train, newdata = x_test)
#   
#   test
# 
# }
# 


