## about ----
## 
## 
## mike freund, 2020-01-26

## setup ----

library(here)
library(magrittr)
library(gifti)
library(cifti)
library(abind)
library(data.table)
library(mikeutils)

## functions

whitening <- function(E, shrinkage = "ledoitwolf") {
  
  S <- cov(E)  ## sample cov
  # corrplot::corrplot(cov2cor(S), method = "color")
  H <- diag(nrow(S))  ## target to shrink towards (Ledoit-Wolf's 'F')
  
  if (shrink %in% c("lw", "ledoitwolf", "LW"))  {
    k <- tawny::shrinkage.intensity(E, H, S)
    lambda <- max(c(0, min(k / nrow(E), 1)))  ## shrinkage factor
  }
  
  S_hat <- lambda * H + (1 - lambda) * S  ## shrunken matrix
  
  solve(S_hat)  ## mahalanobis whitening matrix^2
  
}

## variables

parcellation <- read_atlas("schaefer400")

dir.ccp.hcp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
dir.subsubj <- "/data/nil-external/ccp/witzermanm/AFNI_ANALYSIS_SUBSUBJECT"
dir.results <- file.path(dir.subsubj, "RESULTS_RUNWISE")

subjs <- list.dirs(dir.results, recursive = FALSE, full.names = FALSE)
tasks <- c("Axcpt", "Cuedts", "Stern", "Stroop")
sessi <- c("baseline", "proactive", "reactive")
sessi.short <- c("Bas", "Pro", "Rea")
n.knots <- 8

## create "label table" (lt) for loop indices

lt <- rbind(
  data.table(
    task     = "Axcpt", 
    variable = c("AY", "BX", "BY", "AX", "Bng", "Ang", "error", "button1", "button2"),
    glm.name = c(rep("Cues_EVENTS_censored", 7), rep("Buttons_censored", 2))
  ),
  data.table(
    task     = "Cuedts", 
    variable = c(
      "InConInc", "InconNoInc", "ConInc", "ConNoInc", "error",
      "SwitchInc", "SwitchNoInc", "RepeatInc", "RepeatNoInc",
      "button1", "button2"
    ),
    glm.name = c(
      rep("CongruencyIncentive_EVENTS_censored", 5), 
      rep("SwitchIncentive_EVENTS_censored", 4),
      rep("Buttons_censored", 2)
    )
  ),
  data.table(
    task     = "Stern", 
    variable = c("LL5NP", "LL5NN", "LL5RN", "not5NP", "not5NN", "not5RN", "error", "button1", "button2"),
    glm.name = c(rep( "ListLength_EVENTS_censored", 7), rep("Buttons_censored", 2))
  ),
  data.table(
    task     = "Stroop", 
    variable = c("PC50Con", "PC50InCon", "biasCon", "biasInCon", "error", "button1", "button2"),
    glm.name = c(rep( "ListLength_EVENTS_censored", 5), rep("Buttons_censored", 2))
  )
)

contrast <- list(
  Axcpt  = c("Acue_Bcue", "HI_LO_conf", "Nogo_Go", "error_correct", "B1_B2"),
  Cuedts = c("Inc_NoInc", "InCon_Con", "Switch_Repeat", "error_correct", "B1_B2"),
  Stern  = c("RN_NN_all", "RN_NN_LL5", "RN_NN_not5", "not5_LL5", "error_correct", "B1_B2"),
  Stroop = c("InCon_Con_bias", "InCon_Con_PC50", "InCon_Con_PC50bias", "error_correct")
)

m <- list(
  Axcpt   = rbind(
    ##    ax  ay  ang  bx  by  bng  error
    cbind( 1/2,  1/2,    0, -1/2, -1/2,    0,    0),
    cbind(-1/2,  1/2,    0,  1/2, -1/2,    0,    0),
    cbind(-1/4, -1/4,  1/2, -1/4, -1/4,  1/2,    0),
    cbind(-1/4, -1/4,    0, -1/4, -1/4,    0,  1/2)
  ),
  Cuedts  = rbind(
    ##    ax  ay  ang  bx  by  bng  error
    cbind( 1/2,  1/2,    0, -1/2, -1/2,    0,    0),
    cbind(-1/2,  1/2,    0,  1/2, -1/2,    0,    0),
    cbind(-1/4, -1/4,  1/2, -1/4, -1/4,  1/2,    0),
    cbind(-1/4, -1/4,    0, -1/4, -1/4,    0,  1/2)
  ),
  Stern   = rbind(
  ),
  Stroop  = rbind(
  ),
)

# c("AX", "AY", "Ang", "BX", "BY", "Bng", "error")  ## in same order as in GLMs

for (ii in seq_along(m)) {
  dimnames(m[[ii]]) <- list(contrast = contrast[[task.i]], param = )
}


statistics <- combo.paste(c("cor", "euc"), c("raw", "unn", "mnn"))

## lists of arrays

r <- vector("list", length(tasks)) %>% setNames(tasks)  ## d within
u <- r
for (task.i in seq_along(tasks)) {
  
  r[[task.i]] <- array(
    NA,
    dim = c(
      subj  = length(subjs), 
      sess  = length(sessi),
      roi   = length(parcellation$key), 
      cont  = length(contrast[[task.i]]),
      stat  = length(statistics),
      knot  = n.knots,
      run   = 2
    ),
    dimnames = list(
      subj  = subjs, 
      sess  = sessi,
      roi   = parcellation$key, 
      cont  = contrast[[task.i]],
      stat  = statistics,
      knot  = paste0("knot", seq_len(n.knots)),
      run   = c("run1", "run2")
    )
  )
  
  u[[task.i]] <- array(
    NA,
    dim = c(
      subj  = length(subjs), 
      sess  = length(sessi),
      roi   = length(parcellation$key), 
      cont  = length(contrast[[task.i]]),
      knot  = n.knots,
      run   = 2
    ),
    dimnames = list(
      subj  = subjs, 
      sess  = sessi,
      roi   = parcellation$key, 
      cont  = contrast[[task.i]],
      knot  = paste0("knot", seq_len(n.knots)),
      run   = c("run1", "run2")
    )
  )
  
}


## loop ----


for (task.i in seq_along(tasks)) {
  # task.i = 1
  
  name.task.i <- tasks[task.i]
  
  for (sess.i in seq_along(sessi)) {
    # sess.i = 1
    
    name.sess.i <- sessi[sess.i]
    
    for (subj.i in seq_along(subjs)) {
      ## subj.i = 1
      
      name.subj.i <- subjs[subj.i]
      
      lt.i <- lt[task == name.task.i, c("variable", "glm.name")]
      
      for (name.glm.i in unique(lt.i$glm.name)) {
        # name.glm.i = unique(lt.i$glm.name)[1]
        
        vars.i <- lt.i[glm.name == name.glm.i]$variable
        
        dirs <- combopaste(
          file.path(dir.results, name.subj.i, name.task.i),
          paste0("/", name.sess.i, "/", name.sess.i, "_", name.glm.i)
        )
        
        betas <- vector("list", 4) %>% setNames(combo_paste(c("run1", "run2"), c("L", "R")))
        resid <- betas
        
        for (name.run.i in c("run1", "run2")) {
          # name.run.i = "run1"
          for (name.hemi.i in c("L", "R")) {
            # name.hemi.i = "L"
            
            f.betas <- file.path(paste0(dirs, "_", name.run.i), paste0("stats_", subjs[subj.i], "_", name.hemi.i, ".func.gii"))
            f.resid <- file.path(paste0(dirs, "_", name.run.i), paste0("wherr_", subjs[subj.i], "_", name.hemi.i, ".func.gii"))
            
            file.is.missing <- !file.exists(f.betas) | !file.exists(f.resid)
            if (file.is.missing) next  ## go to next task!
            
            name.run.hemi.i <- paste0(name.run.i, "_", name.hemi.i)
            
            betas[[name.run.hemi.i]] <- collate_surface_params(
              f.betas, 
              pattern = combopaste(vars.i, c("#[0-9]_Coef")) %>% paste0(collapse = "|")
            )
            
            resid[[name.run.hemi.i]] <- read_gifti2matrix(f.resid)
            
          }  ## end hemi loop
        }  ## end run loop
        
        betas <- abind(
          run1 = cbind(betas[["run1_L"]], betas[["run1_R"]]),
          run2 = cbind(betas[["run2_L"]], betas[["run2_L"]]),
          rev.along = 0
        )
        resid <- abind(
          run1 = cbind(resid[["run1_L"]], resid[["run1_R"]]),
          run2 = cbind(resid[["run2_L"]], resid[["run2_L"]]),
          rev.along = 0
        )
        
        
        ## estimation ----
        
        for (roi.i in seq_along(parcellation$key)) {
          # roi.i = 1
          name.roi.i <- parcellation$key[roi.i]
          betas.i <- betas[, parcellation$atlas == roi.i, ]
          
          for (knot.i in seq_len(n.knots)) {
            # knot.i = 1
            
            ## mask (for knot and region)
            
            rows.knot.i <- grep(paste0("#", knot.i - 1, "_Coef"), rownames(betas.i))  ## minus one b/c 0-based ind.
            betas.ii <- betas.i[rows.knot.i, , ]
            
            for (contrast.i in contrast[[task.i]]) {
              
              ## check order of rownames
              
              ## get contrasts
              
              
              contrasts.i <- abind(
                run1 = m %*% betas.ii[, , "run1"],
                run2 = m %*% betas.ii[, , "run2"],
                rev.along = 0
              )
              
              ## estimate
              
              means <- rowMeans(contrasts.i)
              relia <- apply(contrasts.i, 1, function(.) cor(.)[1, 2])
              
              
              
            }
            
            
            
            
          }
          
          
          
        }
        

      }  ## glm loop end
      
      
      
    }  ## subj loop end
    
    
  }  ## end session loop
  
  
  ## save task array
}




## save ----

saveRDS(r, here("out", "runwise", "reliability-contrasts_runwise_schaefer400.rds"))
saveRDS(u, here("out", "runwise", "mean-contrasts_runwise_schaefer400.rds"))


