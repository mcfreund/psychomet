library(here)
library(dplyr)
library(data.table)

if (!exists("dir_atlas")) source(here("code", "_constants.R"))



## schaefer (L first, R second!) ---

fname_schaefer <- file.path(
  dir_atlas,
  "Schaefer2018_Parcellations", "HCP", "fslr32k", "cifti", "Schaefer2018_400Parcels_7Networks_order_info.txt"
)
if (file.exists(fname_schaefer)) {
  fin <- file(fname_schaefer, 'rt')
  tmp <- readLines(fin);
  close(fin); unlink(fin);
  if (length(tmp) != 800) { stop("not expected Schaefer key."); }
  tmp <- tmp[seq(from=1, to=800, by=2)];   # every-other entry is a label
  key_schaefer <- gsub("7Networks_", "", tmp);
}
key_schaefer <- data.table(
  parcel = key_schaefer,
  network = gsub("^.H_(Vis|SomMot|Cont|Default|Limbic|SalVentAttn|DorsAttn)_.*", "\\1", key_schaefer)
)


## mmp (R first, L second!) ----

key_mmp <- fread(
  file.path(dir_atlas, "HCP-MMP", "Glasser_et_al_2016_HCP_MMP1.0_RVVG", "MMP360ParcelsKey.csv")
)
key_mmp <- paste0(gsub("_ROI", "", key_mmp$Community), "_", key_mmp$Hem)[order(key_mmp$ParcelID)]

## network assignments:

coleanticevic <- RCurl::getURL(
  "https://raw.githubusercontent.com/ColeLab/ColeAnticevicNetPartition/master/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt"
)
coleanticevic <- data.table::fread(text = coleanticevic)
coleanticevic <- coleanticevic[!is.na(GLASSERLABELNAME), c("NETWORK", "GLASSERLABELNAME")]
coleanticevic$GLASSERLABELNAME <- gsub("(^.)_(.*)_ROI", "\\2_\\1", coleanticevic$GLASSERLABELNAME)
coleanticevic <- dplyr::rename(coleanticevic, parcel = GLASSERLABELNAME, network = NETWORK)

key_mmp <- data.table(full_join(data.frame(parcel = key_mmp), coleanticevic))  ## match to order from key_mmp


## write ----

fwrite(key_schaefer, here("in", "atlas-key_schaefer400-07.csv"))
fwrite(key_mmp, here("in", "atlas-key_mmp.csv"))
