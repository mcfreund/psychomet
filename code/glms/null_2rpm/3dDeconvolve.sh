#!/usr/bin/env bash

# assumes subject ID is defined as variable $subject
# https://stackoverflow.com/questions/17621798/linux-process-in-background-stopped-in-jobs
# for writing stdin/stdout to file, for background running


function deconvolve_2rpm {
	
    ## build xmat
    /usr/local/pkg/afni_18/3dDeconvolve \
    -local_times \
    -force_TR 1.2 \
    -x1D_stop \
    -input ${name_img_1}" "${name_img_2} \
    -polort A \
    -float \
    -censor ${dir_stimts}/movregs_FD_mask.txt \
    -num_stimts 0 \
    -ortvec ${dir_stimts}/motion_demean_${sessions[$session_i]}.1D movregs \
    -x1D ${dir_out}/X.xmat.1D \
    -xjpeg ${dir_out}/X.jpg \
    -nobucket

}

