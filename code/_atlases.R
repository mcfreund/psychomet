## NB: assumes atlases located in path specified by dir_atlas variable


## get paths if haven't already:
if (!exists("dir_atlas")) source(here::here("code", "_constants.R"))


## read atlas keys ----

not_exist_keys <- 
  !file.exists(here::here("in", "atlas-key_schaefer400-07.csv")) | 
  !file.exists(here::here("in", "atlas-key_mmp.csv"))
if (not_exist_keys) stop("need to write atlas keys. run:\n  source(here::here('in', '_write_atlas_keys.R'))")

key_schaefer <- data.table::fread(here::here("in", "atlas-key_schaefer400-07.csv"))
key_mmp <- data.table::fread(here::here("in", "atlas-key_mmp.csv"))


## read voxel-wise giftis ----

if (Sys.info()["nodename"] == "ccplinux1") {
  
  schaefer10k <-
    c(
      gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_400Parcels_7Networks_order_10K_L.label.gii"))$data[[1]],
      gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_400Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 200
    )
  
  mmp10k <-
    c(
      gifti::read_gifti(
        file.path(dir_atlas, "HCP-MMP_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR_L.label.gii")
        )$data[[1]],
      gifti::read_gifti(
        file.path(dir_atlas, "HCP-MMP_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR_R.label.gii")
        )$data[[1]]
    )

}

  


## ROIs


## DMCC conjunction (baseline test) --- schaefer atlas:
# dmcc34 <- c(
#   22, 77, 78, 86, 87, 91, 93, 99, 101, 103, 105, 107, 110, 127, 130, 139, 140,
#   144, 148, 172, 175, 185, 189, 219, 301, 303, 306, 314, 340, 346, 347, 349, 350, 353
# )

schaefermd <- c(
  99, 127, 129, 130, 131, 132, 137, 140, 141, 142, 148, 163, 165, 182, 186, 300, 332, 333, 334, 335, 336, 337, 340, 345, 
  349, 350, 351, 352, 354, 361, 365, 387
)


## MD [from assem (2020); cereb cort.] --- mmp atlas:

md_core_names <- mikeutils::combo_paste(
  c("p9-46v", "a9-46v", "i6-8", "AVI", "8C", "IFJp", "IP2", "IP1", "PFm", "8BM", "SCEF"),
  c("L", "R")
)


md_core <- match(md_core_names, key_mmp$parcel)
names(md_core) <- md_core_names
md_core <- sort(md_core)
rm(md_core_names)

md_extended_names <- mikeutils::combo_paste(
  c(
    "a9-46v", "p10p", "a10p", "11l", "a47r", "p47r", "FOP5", "AVI", "p9-46v", "8C", "IFJp", "6r", 
    "s6-8", "i6-8", "SCEF", "8BM", "a32pr", "d32", "TE1m", "TE1p", "AIP", "IP2", "LIPd", "MIP", 
    "IP1", "PGs", "PFm", "POS2"
  ),
  c("L", "R")
)
md_extended <- match(md_extended_names, key_mmp$parcel)
names(md_extended) <- md_extended_names
md_extended <- sort(md_extended)
rm(md_extended_names)



