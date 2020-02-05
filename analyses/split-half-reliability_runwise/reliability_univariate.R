library(here)
library(magrittr)
library(gifti)
library(cifti)
library(abind)
library(mikeutils)
library(lme4)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(data.table)
library(tawny)  ## shrinkage.intensity()
library(expm)  ## sqrtm()
library(grid)
library(gridExtra)
library(cowplot)

theme_set(theme_classic(base_size = 8))

d <- readRDS(here("out", "runwise", "data", "estimates-win-run_hilo_baseline_schaefer400.rds"))

nodename <- Sys.info()["nodename"]
if (nodename == "ccplinux1") {
  dir.atlas <- "/data/nil-external/ccp/freund/atlases"
} else if (nodename == "CCP-FREUND") {
  dir.atlas <- "C:/local/atlases"
}

parcellation <- read_atlas("schaefer400")
hcp <- list(
  L  = readGIfTI(
    file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii")
  ),
  R = readGIfTI(
    file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii")
  )
)
over <- list(
  L = parcellation$atlas[1:(nrow(parcellation$atlas) / 2)], 
  R = parcellation$atlas[(nrow(parcellation$atlas) / 2):nrow(parcellation$atlas)]
)

## ----
dimnames(d)
u <- d[, , , "mean", ]  ## u for univariate

subj.has.na.u <- apply(u, "subj", function(.) any(is.na(c(.))))
subj.all.data <- names(subj.has.na.u)[!subj.has.na.u]
u <- u[subj.all.data, , , ]

## create dataframes

u.df <- reshape2::melt(u)
u.df <- cbind(u.df, reshape2::colsplit(u.df$roi, "_", c("hemi", "community", "parcel")))
u.df <- purrr::map_if(u.df, is.factor, as.character) %>% as.data.frame

u.df %<>%
  group_by(roi) %>%
  mutate(roi.num = match(roi, parcellation$key))

u.wide <- u.df %>% tidyr::spread(run, value)

## activation analysis

u.df.byroi <- split(u.df %>% group_by(task), u.df$roi)
# lm(value ~ subj, u.df.byroi[[1]])
mods.activation <- lapply(u.df.byroi, function(.) lmer(value ~ task - 1 + (1 | subj), .))
coef.activation <- lapply(mods.activation, function(.) as.data.frame(coef(summary(.))))
coef.activation %<>% lapply(function(.) tibble::rownames_to_column(., var = "task"))
coef.activation <- as.data.frame(do.call(rbind, coef.activation))
coef.activation$roi <- rep(names(mods.activation), each = 4)
coef.activation$task <- gsub("task", "", coef.activation$task)

coef.activation %<>% cbind(reshape2::colsplit(.$roi, "_", c("hemi", "community", "parcel")))
coef.activation %<>% rename(b = "Estimate", t = "t value", se = "Std. Error")

# coef.activation %>%
#   ggplot(aes(parcel, b)) +
#   facet_grid(vars(task), vars(community), scales = "free", space = "free") +
#   geom_hline(yintercept = 0) +
#   geom_point(aes(alpha = ifelse(abs(t) > 2, 1, 1/12))) +
#   geom_errorbar(
#     aes(
#       ymin = b - se * 1.96, ymax = b + se * 1.96,
#       alpha = ifelse(abs(t) > 2, 1, 1/12)
#     )
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, color = ifelse(abs(coef.activation[["t"]]) > 2, "black", "grey70")),
#     legend.position = "none"
#   )

conjunction <- coef.activation %>%
  group_by(roi) %>%
  summarize(sig = sum(t > 2))
conjunction %>% filter(sig == 2)  ## penumbra
core <- conjunction %>% filter(sig == 3)  ## core
# conjunction %>% filter(sig == 4)

control <- coef.activation %>%
  group_by(roi) %>%
  summarize(ns = sum(t > 0 & t < 1)) %>%
  filter(ns > 2)

cors <- u.wide %>%
  group_by(roi) %>%
  summarize(
    r.spearman = cor(run1, run2, method = "spearman"),
    r.pearson = cor(run1, run2)
  )

ul <- split(u.wide, interaction(u.wide$task, u.wide$roi))

icc.all <- lapply(ul, function(x) psych::ICC(x[c("run1", "run2")], lmer = FALSE)$results)

icc3 <- lapply(icc.all, function(x) x[x$type %in% c("ICC3"), ])
icc3 <- do.call(rbind, icc3)
icc3 <- select(icc3, p, lb = "lower bound", ub = "upper bound", icc = ICC)
icc3 <- cbind(icc3, reshape2::colsplit(rownames(icc3), "\\.", c("task", "roi")))

reliab <- full_join(coef.activation, cors)
reliab <- full_join(reliab, icc3)
reliab$control <- reliab$roi %in% control$roi
reliab$conjunc <- reliab$roi %in% core$roi

write.csv(reliab, here("out", "runwise", "data", "univariate_reliabilities.csv"), row.names = FALSE)

## read ---
reliab <- read.csv(here("out", "runwise", "data", "univariate_reliabilities.csv"))

## plot ----

reliab %>%
  ggplot(aes(t, r.pearson)) +
  geom_point() +
  facet_wrap(vars(task))
reliab %>%
  ggplot(aes(t, r.spearman)) +
  geom_point() +
  facet_wrap(vars(task))
reliab %>%
  ggplot(aes(t, icc)) +
  geom_point() +
  facet_wrap(vars(task))
reliab %>%
  ggplot(aes(icc, r.pearson)) +
  geom_point() +
  facet_wrap(vars(task))


reliab %>%
  ggplot(aes(x = icc)) +
  geom_histogram(aes(fill = t > 2))
reliab %>%
  ggplot(aes(x = r.pearson)) +
  geom_histogram(aes(fill = t > 2))

reliab %>%
  mutate(activated = t > 2.5) %>%
  ggplot(aes(x = activated, y = icc)) +
  geom_boxplot(aes(fill = activated), notch = TRUE) +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    # axis.text.y = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  labs(y = "split-half correlation")

reliab %>%
  ggplot(aes(y = icc)) +
  geom_boxplot(fill = "steelblue", notch = TRUE) +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  labs(y = "split-half correlation")

reliab %>%
  filter(t > 2) %>%
  ggplot(aes(x = community, y = icc)) +
  stat_summary(fun = "mean_cl_boot")
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  labs(y = "split-half correlation")


## workbench plots ----


## writing parcellated overlays for workbech:

# write.pscalar <- function(
#   path.wb = "C:/Program Files/workbench-windows64-v1.2.3/workbench/bin_windows64",
#   dir.to.write = here::here("out", "runwise", "wb"),
#   fname.overlay = "asdf",
#   values.overlay = as.numeric(as.factor(parcellation$key)),
#   name.atlas = "schaefer"
# ) {
#   ## adapted from Jo Etzel, 2019-02-26
#   ## see: http://mvpa.blogspot.com/2017/11/assigning-arbitrary-values-to-surface.html
#   
#   if (!dir.exists(dir.to.write)) stop(paste0("missing: ", dir.to.write))
#   
#   ## location of wb_command.exe, to call wb_command functions via system() within R:
#   dir.origwd <- getwd()
#   on.exit(setwd(dir.origwd))  ## return to original wd upon exit
#   setwd(path.wb)
#   
#   ## file names and dirs
#   dir.s1200      <- file.path(dir.atlas, "surf/HCP_S1200_GroupAvg_v1") ## HCP S1200 Group Average Data Release
#   fname.template <- file.path(dir.to.write, "template.pscalar.nii")  ## filename for the template we'll make
#   fname.text     <- file.path(dir.to.write, paste0(fname.overlay, ".txt"))
#   fname.cifti    <- file.path(dir.to.write, paste0(fname.overlay, ".pscalar.nii"))
#   
#   ## make a template pscalar.nii from the atlas
#   
#   name.atlas <- tolower(name.atlas)
#   
#   if (name.atlas == "schaefer400") {
#     
#     fname.dlabel <- file.path(
#       dir.atlas,
#       "Schaefer2018_Parcellations/HCP/fslr32k/cifti/Schaefer2018_400Parcels_7Networks_order.dlabel.nii"
#     )
#     
#   } else if (name.atlas == "glasser") {
#     
#     fname.dlabel <- file.path(
#       dir.s1200, 
#       "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii"
#       )
#     
#   } else if (name.atlas == "gordon") {
#     
#     fname.dlabel <- file.path(dir.atlas, "gordon/gordon_parcels/Parcels/Parcels_LR.dlabel.nii")
#     
#   }
#   
#   stdout.template <- system(
#     paste0(
#       "wb_command -cifti-parcellate ", dir.s1200, "/S1200.thickness_MSMAll.32k_fs_LR.dscalar.nii ", 
#       fname.dlabel, " COLUMN ",
#       fname.template
#     )
#   )
#   if (stdout.template != 0) stop("problem writing template")
#   if (!file.exists(fname.template)) stop(paste0("missing: ", fname.template))
# 
#   ## the text file needs to be arranged with the right hemisphere parcels in the first N/2 rows (in order),   
#   ## then the N/2 parcels for the left hemisphere.
#   write.table(values.overlay, fname.text, col.names = FALSE, row.names = FALSE)
#   
#   ## create a CIFTI from the text file for viewing in Workbench  
#   stdout.cifti <- system(
#     paste0("wb_command -cifti-convert -from-text ", fname.text, " ", fname.template, " ", fname.cifti)
#   )
#   if (stdout.cifti != 0) stop("problem writing cifti")
#   if (!file.exists(fname.cifti)) stop(paste("missing:", fname.cifti))
#   
#   ## remove text file
#   was.removed <- file.remove(file.path(dir.to.write, paste0(fname.overlay, ".txt")))
#   if (!was.removed) stop("trouble removing text file")
# 
# }
# 

cifti.parcellate <- function(
  path.wb = "C:/Program Files/workbench-windows64-v1.2.3/workbench/bin_windows64",
  dir.to.write = here::here("out", "runwise", "wb"),
  name.atlas = "schaefer400"
) {
  ## adapted from Jo Etzel, 2019-02-26
  ## see: http://mvpa.blogspot.com/2017/11/assigning-arbitrary-values-to-surface.html
  
  if (!dir.exists(dir.to.write)) stop(paste0("missing: ", dir.to.write))
  
  ## location of wb_command.exe, to call wb_command functions via system() within R:
  dir.origwd <- getwd()
  on.exit(setwd(dir.origwd))  ## return to original wd upon exit
  setwd(path.wb)
  
  ## file names and dirs
  name.atlas <- tolower(name.atlas)
  dir.s1200 <- file.path(dir.atlas, "surf/HCP_S1200_GroupAvg_v1") ## HCP S1200 Group Average Data Release
  
  if (name.atlas == "schaefer400") {
    
    fname.dlabel <- "Schaefer2018_400Parcels_7Networks_order"
    fname.dlabel.full <- file.path(
      dir.atlas, "Schaefer2018_Parcellations", "HCP", "fslr32k", "cifti", paste0(fname.dlabel, ".dlabel.nii")
    )
    
  } else if (name.atlas == "glasser") {
    
    fname.dlabel <- "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR"
    fname.dlabel.full <- file.path(dir.s1200, paste0(fname.dlabel, ".dlabel.nii"))
    
  } else if (name.atlas == "gordon") {
    
    fname.dlabel <- "Parcels_LR"
    fname.dlabel.full <- file.path(
      dir.atlas, "gordon", "gordon_parcels", "Parcels", paste0(fname.dlabel, ".dlabel.nii")
    )
    
  } else {stop("did not find name.atlas")}
  
  fname.template <- file.path(dir.to.write, paste0(fname.dlabel, ".pscalar.nii"))
  
  stdout.template <- system(
    paste0(
      "wb_command -cifti-parcellate ", dir.s1200, "/S1200.thickness_MSMAll.32k_fs_LR.dscalar.nii ", 
      fname.dlabel.full, " COLUMN ",
      fname.template
    )
  )
  
  if (stdout.template != 0) stop("problem writing template")
  if (!file.exists(fname.template)) stop(paste0("missing: ", fname.template))
  
  return(paste0("wrote pscalar for ", name.atlas))
  
}

# cifti.parcellate()
# cifti.parcellate(name.atlas = "gordon")
# cifti.parcellate(name.atlas = "glasser")

cifti.convert <- function(
  path.wb = "C:/Program Files/workbench-windows64-v1.2.3/workbench/bin_windows64",
  dir.to.write = here::here("out", "runwise", "wb"),
  fname.overlay = "asdf",
  values.overlay = as.numeric(as.factor(parcellation$key)),
  name.atlas = "schaefer400",
  dir.template
) {
  ## adapted from Jo Etzel, 2019-02-26
  ## see: http://mvpa.blogspot.com/2017/11/assigning-arbitrary-values-to-surface.html
  
  if (!dir.exists(dir.to.write)) stop(paste0("missing: ", dir.to.write))
  
  ## location of wb_command.exe, to call wb_command functions via system() within R:
  dir.origwd <- getwd()
  on.exit(setwd(dir.origwd))  ## return to original wd upon exit
  setwd(path.wb)
  
  ## file names and dirs
  dir.s1200 <- file.path(dir.atlas, "surf/HCP_S1200_GroupAvg_v1") ## HCP S1200 Group Average Data Release
  fname.text <- file.path(dir.to.write, paste0(fname.overlay, ".txt"))
  fname.cifti <- file.path(dir.to.write, paste0(fname.overlay, ".pscalar.nii"))
  
  if (name.atlas == "schaefer400") {
    
    fname.template <- file.path(dir.template, "Schaefer2018_400Parcels_7Networks_order.pscalar.nii")
    
  } else if (name.atlas == "glasser") {
    
    fname.template <- file.path(
      dir.template, 
      "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.pscalar.nii"
    )
    
  } else if (name.atlas == "gordon") {
    
    fname.template <- file.path(dir.template, "Parcels_LR.pscalar.nii")
    
  } else {stop("did not find name.atlas")}
  
  if (!file.exists(fname.template)) stop(paste0("missing: ", fname.template))
  
  ## the text file needs to be arranged with the right hemisphere parcels in the first N/2 rows (in order),   
  ## then the N/2 parcels for the left hemisphere.
  write.table(values.overlay, fname.text, col.names = FALSE, row.names = FALSE)
  
  ## create a CIFTI from the text file for viewing in Workbench  
  stdout.cifti <- system(
    paste0("wb_command -cifti-convert -from-text ", fname.text, " ", fname.template, " ", fname.cifti)
  )
  if (stdout.cifti != 0) stop("problem writing cifti")
  if (!file.exists(fname.cifti)) stop(paste("missing:", fname.cifti))
  
  ## remove text file
  was.removed <- file.remove(file.path(dir.to.write, paste0(fname.overlay, ".txt")))
  if (!was.removed) stop("trouble removing text file")
  
}


reliab.max <- reliab %>% 
  group_by(roi, hemi, control, conjunc) %>%
  summarize(icc = max(icc))
reliab.max$roi.num <- match(reliab.max$roi, parcellation$key)
reliab.max %<>% arrange(roi.num)


reliab.max.conjunc <- reliab.max
reliab.max.conjunc[!reliab.max.conjunc$conjunc, "icc"] <- 0
cifti.convert(
  fname.overlay = "runwise_reliability_icc_thresh_conjunc",
  values.overlay = reliab.max.conjunc$icc,
  dir.template = here::here("out", "runwise", "wb")
  )

reliab.max.noeffect <- reliab.max
reliab.max.noeffect[!reliab.max.noeffect$control, "icc"] <- 0
cifti.convert(
  fname.overlay = "runwise_reliability_icc_thresh_noeffect",
  values.overlay = reliab.max.noeffect$icc,
  dir.template = here::here("out", "runwise", "wb")
)

largest.t <- reliab %>%
  group_by(roi) %>%
  summarize(t = mean(t)) %>%
  arrange(-t) %>%
  .[1:25, ] %>%
  .$roi

reliab.max.largest <- reliab.max
reliab.max.largest[!reliab.max.largest$roi %in% largest.t, "icc"] <- 0
cifti.convert(
  fname.overlay = "runwise_reliability_icc_thresh_largest",
  values.overlay = reliab.max.largest$icc,
  dir.template = here::here("out", "runwise", "wb")
)



reliab.mean <- reliab %>% 
  group_by(roi, hemi, control, conjunc) %>%
  summarize(icc = tanh(mean(atanh(icc))))
reliab.mean$roi.num <- match(reliab.mean$roi, parcellation$key)
reliab.mean %<>% arrange(roi.num)

cifti.convert(
  fname.overlay = "runwise_reliability_icc_mean_unthresh",
  values.overlay = reliab.mean$icc,
  dir.template = here::here("out", "runwise", "wb")
)



##
library(magick)
brains <- image_read(
  here::here("out", "runwise", "wb", "split-half-icc-max_threshold-largest-smallest.png")
  )

theme_set(theme_classic(base_size = 12))

boxplot.icc <- reliab %>%
  ggplot(aes(y = icc)) +
  geom_boxplot(fill = "#1b9e77", notch = TRUE) +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  scale_y_continuous(breaks = c(round(min(reliab$icc), 2), round(max(reliab$icc), 2))) +
  labs(y = "split-half correlation (ICC)")

p <- grid.arrange(
  rasterGrob(brains),
  boxplot.icc,
  ncol = 2,
  widths = c(1, 0.15)
)

plot_grid(p)
ggsave(
  here::here("out", "runwise", "wb", "split-half-icc-max_threshold-largest-smallest.pdf"), 
  p,
  # device = "cairo",
  width = 6,
  height = 4
  )
