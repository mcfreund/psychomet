#source(here::here("code", "_packages.R"))
#source(here("code", "read-behav.R"))
# library(here)
# library(dplyr)
# library(data.table)
# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(abind)
# library(doParallel)
# library(foreach)
#library(mikeutils)

source(here("code", "_constants.R"))
# source(here("code", "_constants.R"))
source(here::here("..", "ub55", "code", "_packages.R"))
source(here("..", "ub55", "code", "_vars.R"))
source(here("..", "ub55", "code", "_atlases.R"))
source(here("code", "_settings.R"))


subjs <- subjs_ub55[!subjs_ub55 %in% "432332"]

glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name_glm = c(
    "baseline_Cues_EVENTS_censored_shifted",
    "baseline_CongruencySwitch_EVENTS_censored_shifted",
    "baseline_ListLength_EVENTS_censored_shifted",
    "baseline_Congruency_EVENTS_censored_shifted"
  ),
  stringsAsFactors = FALSE
)
glminfo <- as.data.table(glminfo)
hi <- list(
  Axcpt = "BX",
  Cuedts = c("InConInc", "InConNoInc"),
  Stern = "LL5RN",
  Stroop = c("biasInCon", "PC50InCon")
)
lo <- list(
  Axcpt = "BY",
  Cuedts = c("ConInc", "ConNoInc"),
  Stern = "LL5NN",
  Stroop = c("biasCon", "PC50Con")
)



## 1. get contrast arrays
## 2. loop over tasks*ROIs
## 3. compute distance matrices and save.


## setup output directory

out_dir <- here("out", "id")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


for (glm_i in seq_len(nrow(glminfo))) {
	# glm_i = 1
	
	name_glm_i <- glminfo[glm_i]$name_glm
	name_task_i <- glminfo[glm_i]$task	  
	
	fig_dir <- file.path(out_dir, "dists_figs", paste0(name_task_i, "_", name_glm_i))
	if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
	
	## prepare betas
	
	betas_i <- readRDS(
		here("..", "ub55", "out", "glms", paste0("betas_", glminfo[glm_i]$task, "_", glminfo[glm_i]$name_glm,  ".RDS"))
	)
	betas_i <- betas_i[, , , !dimnames(betas_i)$subj %in% "432332", ]  ## remove subj with missing data
	
	## average across target TRs: 
	betas_i <- abind(
		apply(betas_i[, hi[[glm_i]], target_trs[[glm_i]], , ], c("vertex", "subj", "run"), mean),
		apply(betas_i[, lo[[glm_i]], target_trs[[glm_i]], , ], c("vertex", "subj", "run"), mean),
		along = 0
	)   ## condition, vertex, subj, run
	names(dimnames(betas_i)) <- c("condition", "vertex", "subj", "run")
	dimnames(betas_i)$condition <- c("hi", "lo")
	
	
	## estimate distance matrices
	
	cl <- makeCluster(n_cores/2)
	registerDoParallel(cl)
	time.start <- Sys.time()
	
	z <- foreach(
	  roi_i = seq_along(parcellation$key), .verbose = TRUE,
	  .packages = c("mikeutils", "here", "data.table", "ggplot2", "grid", "gridExtra", "dplyr", "abind")
	  ) %dopar% {
      # roi_i = 1
      
		is_roi <- schaefer10k == roi_i
		B <- betas_i[, is_roi, , ]
	  
		D <- 
			array(
				NA, 
				dim = c(length(subjs), length(subjs), 2), 
				dimnames = list(
					subj_run1 = subjs, subj_run2 = subjs, type = c("univariate", "multivariate")
				)
			)
			
		B_contrast <- B["hi", , , ] - B["lo", , , ]  ## get contrast

		B1 <- t(B_contrast[, , 1])  ## separate by runs
		B2 <- t(B_contrast[, , 2])

		v <- nrow(B_contrast)  ## number vertices
		B1_bar <- cbind(rowSums(B1))  ## across-vertex sum of b-coefficients. one per subj.
		B2_bar <- cbind(rowSums(B2))

		D[, , "multivariate"] <- pdist2(B1, B2) / v
		D[, , "univariate"] <- pdist2(B1_bar, B2_bar) / v
      
      
		## save 
		
		p <- arrangeGrob(
			mikeutils::matplot(D[, , "multivariate"]) + labs(title = "multivariate") + theme(legend.position = "left"),
			mikeutils::matplot(D[, , "univariate"]) + labs(title = "univariate") + theme(legend.position = "left"),
			top = parcellation$key[roi_i],
			nrow = 1
		)
      
		ggsave(
			file.path(fig_dir, paste0("euclidean_", parcellation$key[roi_i], ".pdf")), 
			p,
			width = 14, height = 8, units = "cm",
			device = "pdf"
			)
		
		D
		
	}
	
	stopCluster(cl)
	(time.run <- Sys.time() - time.start)
	
	
	## collate and save
	
	D <- 
		array(
			NA, 
			dim = c(length(subjs), length(subjs), length(parcellation$key), 2), 
			dimnames = list(
				subj_run1 = subjs, subj_run2 = subjs, roi = parcellation$key, type = c("univariate", "multivariate")
			)
	)
	
	for (roi_i in seq_along(parcellation$key)) D[, , roi_i, ] <- z[[roi_i]]

	saveRDS(D, file.path(out_dir, paste0("euclidean_",  name_task_i, "_", name_glm_i, ".RDS")))
	
	
}
	
  
res <- apply(D, c("roi", "type"), idi)

plot(res)
abline(0, 1)


m <- D[, , 78, 2]
m_rank <- m
m_rank[] <- rank(m)


m_row <- apply(m, 1, function(x) which.max(x))
m_col <- apply(m, 2, function(x) which.max(x))

mean(m_row == seq_along(m_row))
mean(m_col == seq_along(m_col))


arrangeGrob(
			mikeutils::matplot(D[, , 135, "multivariate"]) + labs(title = "multivariate") + theme(legend.position = "left"),
			mikeutils::matplot(D[, , "univariate"]) + labs(title = "univariate") + theme(legend.position = "left"),
			top = parcellation$key[roi_i],
			nrow = 1
		)
		
## network
## aggregate across tasks*ROIs
## scale by patterns
## implement fingerprinting version
## prewhiten

