
## get directories

if (Sys.info()["nodename"] == "ccplinux1") {
  
  dir.atlas <- "/data/nil-external/ccp/freund/atlases"
  dir.schaefer <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
  
} else if (Sys.info()["nodename"] == "CCP-FREUND") {
  ## mike freund's (i.e., ccp's) thinkpad
  ## reliant on box drive
  ## assumes box drive location at ./Users/mcf/Box
  
  dir.atlas <- "C:/local/atlases"
  dir.schaefer <- dir.atlas
  
} else if (Sys.info()["nodename"] == "PUTER") {
  
  dir.atlas <- "C:/Users/mcf/Documents/atlases"
  dir.schaefer <- file.path(dir.atlas, "ATLASES")
  
}


## read


library(cifti)
parcellation <- mikeutils::read_atlas("schaefer400", path.atlas = dir.schaefer)


if (Sys.info()["nodename"] %in% c("ccplinux1", "PUTER")) {
  
  schaefer10k <-
    c(
      gifti::read_gifti(file.path(dir.schaefer, "Schaefer2018_400Parcels_7Networks_order_10K_L.label.gii"))$data[[1]],
      gifti::read_gifti(file.path(dir.schaefer, "Schaefer2018_400Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 200
    )
  
}



## write:
##  _write_atlases.R
##  _atlases.R
## MMP-MD parcels
## Schaefer-network
## Schaefer-parcels
## ....
