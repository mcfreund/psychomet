## about ----


## setup ----

library(here)
library(magrittr)
library(mikeutils)
library(cifti)
library(gifti)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(data.table)
library(grid)
library(gridExtra)
library(cowplot)
library(multcomp)
library(lme4)

theme_set(theme_classic(base_size = 8))


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

reliab <- read.csv(here("out", "runwise", "data", "univariate_reliabilities.csv"))

d <- readRDS(here("out", "runwise", "data", "estimates-win-run_hilo_baseline_schaefer400.rds"))
u <- d[, , , "mean", ]  ## u for univariate
subj.has.na.u <- apply(u, "subj", function(.) any(is.na(c(.))))
subj.all.data <- names(subj.has.na.u)[!subj.has.na.u]
u <- u[subj.all.data, , , ]

## create dataframes

u <- reshape2::melt(u)
u <- cbind(u, reshape2::colsplit(u$roi, "_", c("hemi", "community", "parcel")))
u <- purrr::map_if(u, is.factor, as.character) %>% as.data.frame
u %<>%
  group_by(roi) %>%
  mutate(roi.num = match(roi, parcellation$key))

u.wide <- u %>% tidyr::spread(run, value)


## unvariate analysis ----

stats <- u %>%
  filter(task != "Cuedts") %>%
  group_by(roi, subj) %>%
  summarize(value = mean(value)) %>%
  summarize(
    v = wilcox.test(value, alternative = "greater")$statistic,
    p = wilcox.test(value, alternative = "greater")$p.value
  ) %>%
  as.data.frame

stats %<>% cbind(., reshape2::colsplit(stats$roi, "_", c("hemi", "community", "parcel")))
stats %<>% mutate(p.adj = p.adjust(p, method = "holm"))

activ <- stats %>% filter(p.adj < 0.05)


## plot ICC

# forplot <- bind_rows(
#   reliab %>% 
#     filter(roi %in% activ$roi) %>%
#     group_by(roi) %>%
#     summarize(icc = tanh(mean(atanh(icc))), community = "activated"),
#   reliab %>%
#     group_by(community, roi) %>%
#     summarize(icc = tanh(mean(atanh(icc))))
# )

# maximums <- full_join(
#   forplot,
#   u %>% group_by(roi) %>% summarize(maximum = max(value))
# )

reliab.plus <- bind_rows(
  reliab %>% filter(roi %in% activ$roi) %>% mutate(community = "activated"), 
  reliab
  )

icc.means <- reliab.plus %>%
  group_by(roi, community) %>%
  summarize(icc = tanh(mean(atanh(icc))))

icc.means <- full_join(
  icc.means,
  reliab.plus %>% group_by(community) %>% summarize(maximum = max(icc))
)


# https://stackoverflow.com/questions/52034747/plot-only-one-side-half-of-the-violin-plot

yaxis.order <- icc.means %>% group_by(community) %>% summarize(icc = mean(icc)) %>% arrange(icc) %>% .$community
icc.means$community <- factor(icc.means$community, levels = yaxis.order)

icc.means %>%
  ggplot(aes(x = community, y = icc, color = community == "activated")) +
  ggbeeswarm::geom_quasirandom(width = 0.1, alpha = 0.05) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, size = 1) +
  geom_point(aes(y = maximum)) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(color = ifelse(sort(unique(icc.means$community)) == "activated", "firebrick", "black")),
  ) +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "black")) +
  labs(y = "relative consistency across scanning runs (ICC)") +
  coord_flip() +
  geom_segment(aes(y = min(icc), yend = maximum, x = -Inf, xend = -Inf)) +
  scale_y_continuous(breaks = c(round(min(icc.means$icc), 2), 0, 0.25, round(max(icc.means$maximum), 2)))
  


## conjunction (across tasks) ----

stats.bytask <- u %>%
  filter(task != "Cuedts") %>%
  group_by(roi, subj, task) %>%
  summarize(value = mean(value)) %>%
  group_by(roi, task) %>%
  summarize(
    v = wilcox.test(value, alternative = "greater")$statistic,
    p = wilcox.test(value, alternative = "greater")$p.value
  ) %>%
  as.data.frame

stats.bytask %<>% cbind(., reshape2::colsplit(stats.bytask$roi, "_", c("hemi", "community", "parcel")))
stats.bytask %<>% mutate(p.adj = p.adjust(p, method = "holm"))

conj <- stats.bytask %>%
  group_by(roi) %>%
  summarize(sig = sum(p < 0.05))
conj %<>% filter(sig == 3)


reliab.plus.conj <- bind_rows(
  reliab %>% filter(roi %in% conj$roi, task != "Cuedts") %>% mutate(community = "activated"), 
  reliab
)

icc.means.conj <- reliab.plus.conj %>%
  group_by(roi, community) %>%
  summarize(icc = tanh(mean(atanh(icc))))

icc.means.conj <- full_join(
  icc.means.conj,
  reliab.plus.conj %>% group_by(community) %>% summarize(maximum = max(icc))
)

yaxis.order.conj <- icc.means.conj %>% group_by(community) %>% summarize(icc = mean(icc)) %>% arrange(icc) %>% .$community
icc.means.conj$community <- factor(icc.means.conj$community, levels = yaxis.order.conj)

p.conj <- icc.means.conj %>%
  ggplot(aes(x = community, y = icc, color = community == "activated")) +
  geom_hline(yintercept = 0, color = "grey50", alpha = 0.5) +
  ggbeeswarm::geom_quasirandom(width = 0.1, alpha = 0.05) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, size = 1) +
  geom_point(aes(y = maximum)) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(
      color = ifelse(sort(unique(icc.means.conj$community)) == "activated", "firebrick", "black"),
      face = ifelse(sort(unique(icc.means.conj$community)) == "activated", "bold", "plain")
      ),
  ) +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "black")) +
  labs(y = "relative consistency across scanning runs (ICC)") +
  coord_flip() +
  geom_segment(aes(y = min(icc), yend = maximum, x = -Inf, xend = -Inf)) +
  scale_y_continuous(
    breaks = c(
      round(min(icc.means.conj$icc), 2), 0, 0.25, round(max(icc.means.conj$maximum), 2)
      )
    ) +
  scale_x_discrete(
    labels = c(
      DorsAttn = "dorsal attention", Cont = "fronto-parietal", SalVentAttn = "ventral attention",
      Default = "default", Limbic = "limbic", Vis = "visual", "SomMot" = "somato-motor"
    )
  )

ggsave(here("out", "runwise", "figs", "splithalf-icc_conjunction-sans-cts.pdf"), p.conj)



## TODO:
##  1. define activated (dmcc30?)
##  2. plot
##  3. workbench figure
##  
##   - perform contrast (anova) on each parcel
##   - remove cts?
##   - rank by effect in plot
##   - color activated
##   - wb, mean ICC, activated threshold with same
  
  
## workbench plots ----


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

