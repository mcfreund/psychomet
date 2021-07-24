
## get paths if haven't already:
if (!(exists("dir_schaefer") | !exists("dir_mmp"))) source(here("code", "_constants.R"))



# parcellation <- mikeutils::read_atlas("schaefer400", path.atlas = dir_schaefer)

## if need to write atlas info to ./in dir...

if (nodename == "ccplinux1") {
  
  
  
}

fname.atlas <- file.path(
  path.atlas,
  "Schaefer2018_Parcellations", "HCP", "fslr32k", "cifti", "Schaefer2018_400Parcels_7Networks_order_info.txt"
)
if (file.exists(fname.atlas)) {
  fin <- file(fname.atlas, 'rt')
  tmp <- readLines(fin);
  close(fin); unlink(fin);
  if (length(tmp) != 800) { stop("not expected Schaefer key."); }
  tmp <- tmp[seq(from=1, to=800, by=2)];   # every-other entry is a label
  atlas.key <- gsub("7Networks_", "", tmp);
}


## voxel-wise giftis

if (nodename %in% c("ccplinux1", "PUTER")) {
  
  schaefer10k <-
    c(
      gifti::read_gifti(file.path(dir_schaefer, "Schaefer2018_400Parcels_7Networks_order_10K_L.label.gii"))$data[[1]],
      gifti::read_gifti(file.path(dir_schaefer, "Schaefer2018_400Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 200
    )
  
  mmp10k <-
    c(
      gifti::read_gifti(file.path(dir_mmp, "HCP-MMP_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR_L.label.gii"))$data[[1]],
      gifti::read_gifti(file.path(dir_mmp, "HCP-MMP_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR_R.label.gii"))$data[[1]]
    )
  
  
}

## 