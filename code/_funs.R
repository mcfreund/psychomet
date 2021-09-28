enlist <- function(x) setNames(vector("list", length(x)), x)

## NB: following function is obsolete... can use abind::abind() for this...
# list2array <- function(l, name_new_dim = NULL) {
# 	## input: a list of n-dimensional arrays, of identical dimension and dimnames.
# 	## binds arrays along new dimension.
# 	## output: a new, (n+1)-dimensional array.
# 	## 
# 	## works well for wrangling results of foreach::foreach() call, when n-D array is returned each iteration.
# 	
# 	
# 	## TODO: input validation --- check list, and all dims and dimnames equal
# 	
# 	## get info existing dims
# 	
# 	dim_existing <- dim(l[[1]]) ## vector with values equal to length of each dim
# 	dimnames_existing <- dimnames(l[[1]])  ## list of dim names
# 	dim_new <- length(l)  ## length of new dim
# 	dimnames_new <- setNames(list(names(l)), name_new_dim)
# 	
# 	## make new array:
# 	
# 	a <- array(
# 		NA,
# 		dim = c(dim_existing, dim_new),
# 		dimnames = c(dimnames_existing, dimnames_new)
# 	)
# 	
# 	## iteratively add elements of new dim to array:
# 	
# 	index_existing <- as.matrix(expand.grid(lapply(dim_existing, seq_len)))
# 	for (i in seq_len(dim_new)) a[cbind(index_existing, i)] <- l[[i]]
# 	
# 	a
# 	
# }



pdist2 <- function(A,B) {
  
  ## this function works on matrices A and B.
  ## A and B should be matrices with subjects as rows and vertices as columns.
  ## A and B should be from separate folds (i.e., "test" and "retest" or "run1" and "run2").
  
  ## this function computes the squared euclidean distances between each row of A and each row of B.
  ## the output is therefore a matrix of size nrow(A)*nrow(B), i.e., Sects * Sects.
  
  ## this is an efficient implementation of the pdist::pdist function.
  ## see links for more information:
  ## https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
  ## https://blog.smola.org/post/969195661/in-praise-of-the-second-binomial-formula
  
  
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  
  m = nrow(A)
  n = nrow(B)
  
  tmp = matrix(rep(an, n), nrow=m) 
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  
  tmp - 2 * tcrossprod(A,B)  ## squared euclidean distance
  
}




idi <- function(x) {
  
  ## this function works on square distance matrices, and computes the "intersubject discrimination index".
  ## it subtracts the mean diagonal from the mean of the off diagonals.
  
  offdiag <- x[row(x) != col(x)]
  diagonal <- diag(x)
  
  mean(offdiag) - mean(diagonal)
  
}




read_results <- function(waves, tasks, sessions, subjs, glmname, filename, .read) {
  
  
  l <- enlist(mikeutils::combo_paste(waves, tasks, sessions, subjs))
  
  for (wave_i in seq_along(waves)) {
    
    name_wave_i <- waves[wave_i]
    
    for (task_i in seq_along(tasks)) {
      
      name_task_i <- tasks[task_i]
      
      for (session_i in seq_along(sessions)) {
        
        name_session_i <- sessions[session_i]
        
        for (subj_i in seq_along(subjs)) {
          
          name_subj_i <- subjs[subj_i]
          
          filename_full <- 
            here::here(
              "out", "glms", name_subj_i, "RESULTS", name_task_i, 
              paste0(name_session_i, "_", glmname, "_", name_wave_i),
              filename
              )
          
          nm <- paste0(name_wave_i, "_", name_task_i, "_", name_session_i, "_", name_subj_i)
          l[[nm]] <- .read(filename_full)
          
          
        }
        
      }
      
    }
    
  }
  
  l
  
}



read_betas_dmcc <- function(
  .subjs,
  .task,
  .glm,
  .dir
) {
  # .subjs = "130518"
  # .task = "Axcpt"
  # .glm = "baseline_Cues_EVENTS_censored"
  # .dir = "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS"
  
  ## initialize array
  
  pick.a.file <- 
    file.path(.dir, .subjs[1], "1TRpK_SURFACE_RESULTS",  .task, paste0(.glm), paste0("STATS_", subjs[1], "_REML_L.func.gii"))
  labs <- afni("3dinfo", paste0("-label ", pick.a.file))
  labs <- unlist(strsplit(labs, "\\|"))
  is.reg <- !grepl("Full|block|Tstat|Fstat", labs)
  tab <- do.call(rbind, strsplit(gsub("_Coef", "", labs[is.reg]), "#"))
  trs <- as.numeric(unique(tab[, 2])) + 1
  regs <- unique(tab[, 1])
  
  n.vertex <- 10242
  n.tr <- length(trs)
  n.reg <- length(regs)
  n.subj <- length(subjs)
  
  betas <- array(
    NA,
    dim = c(n.vertex*2, n.reg, n.tr, n.subj),
    dimnames = list(vertex = NULL, reg = regs, tr = NULL, subj = subjs)
  )
  
  vertex.inds <- cbind(L = 1:n.vertex, R = (n.vertex + 1):(n.vertex * 2))
  
  for (subj.i in seq_along(subjs)) {
    # subj.i = 41; hemi.i = "L"
    
    for (hemi.i in c("L", "R")) {
      # hemi.i = "R"
      
      inds <- vertex.inds[, hemi.i]
      
      fname <- file.path(
        .dir, .subjs[subj.i], "1TRpK_SURFACE_RESULTS",  .task, paste0(.glm),  
        paste0("STATS_", subjs[subj.i], "_REML_", hemi.i, ".func.gii")
      )
      
      if (!file.exists(fname)) next
      
      B <- mikeutils::read_gifti2matrix(fname)[is.reg, ]
      
      is.ok.i <- isTRUE(all.equal(dim(B), c(n.reg * n.tr, n.vertex)))
      if (!is.ok.i) stop("mismatched beta array")
      
      
      for (reg.i in seq_len(n.reg)) {
        # reg.i = 1
        
        is.reg.i <- grepl(paste0("^", regs[reg.i], "#"), labs[is.reg])
        B.reg.i <- t(B[is.reg.i, ])
        
        is.ok.ii <- isTRUE(all.equal(dim(betas[inds, reg.i, , subj.i]), dim(B.reg.i)))
        if (!is.ok.ii) stop("mismatched regressor array")
        
        betas[inds, reg.i, , subj.i] <- B.reg.i
        
      }
      
    }
    
  }
  
  betas
  
}



# 
# 
# read_results_bluearc <- function(
#   waves, 
#   tasks , 
#   sessions, 
#   subjs, 
#   dir_child = "1TRpK_SURFACE_RESULTS", 
#   glmnames = 
#     c(
#       Axcpt = "Cues_EVENTS_censored",
#       Cuedts = "CongruencyIncentive_EVENTS_censored",
#       Stern = "ListLength_EVENTS_censored",
#       Stroop = "Congruency_EVENTS_censored"
#     ), 
#   filename, 
#   .read
#   ) {
#   
#   ## dir_child: subdirectory of subject's directory
#   ## glmname: if 
#   
#   
#   l <- enlist(mikeutils::combo_paste(waves, tasks, sessions, subjs))
#   
#   for (wave_i in seq_along(waves)) {
#     
#     name_wave_i <- waves[wave_i]
#     
#     for (task_i in seq_along(tasks)) {
#       
#       name_task_i <- tasks[task_i]
#       
#       for (session_i in seq_along(sessions)) {
#         
#         name_session_i <- sessions[session_i]
#         
#         for (subj_i in seq_along(subjs)) {
#           
#           name_subj_i <- subjs[subj_i]
#           
#           wavedir <- c(wave1 = "HCP_SUBJECTS_BACKUPS", wave2 = "DMCC_Phase3", wave3 = "DMCC_Phase4")[waves[wave_i]]
#           
#           list.files(
#             file.path(
#               "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS", wavedir,
#               "fMRIPrep_AFNI_ANALYSIS", name_subj_i, dir_child, paste0(name_session_i, glmnames[task_i])
#               )
#           )
#           
#           list.files(
#             file.path(
#               "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS", wavedir,
#               "fMRIPrep_AFNI_ANALYSIS"
#             )
#           )
#           
#           # filename_full <- 
#           #   here::here(
#           #     "out", "glms", name_subj_i, "RESULTS", name_task_i, 
#           #     paste0(name_session_i, "_", glmname, "_", name_wave_i),
#           #     filename
#           #   )
#           
#           nm <- paste0(name_wave_i, "_", name_task_i, "_", name_session_i, "_", name_subj_i)
#           l[[nm]] <- .read(filename_full)
#           
#           
#         }
#         
#       }
#       
#     }
#     
#   }
#   
#   l
#   
# }
# 
# 





get_network <- function(x) {
  gsub("^.H_(Vis|SomMot|Cont|Default|Limbic|SalVentAttn|DorsAttn)_.*", "\\1", x)
}


loads <- function(x, dims = c("PC1", "PC2")) {
  
  rot <- as.data.frame(x$rotation)[, dims]
  rot$variable <- rownames(rot)
  
  rot
  
}
