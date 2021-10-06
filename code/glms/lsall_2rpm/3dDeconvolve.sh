#!/usr/bin/env bash

# assumes subject ID is defined as variable $subject
# https://stackoverflow.com/questions/17621798/linux-process-in-background-stopped-in-jobs
# for writing stdin/stdout to file, for background running


function deconvolve_2rpm_lsall_stroop {
	
    ## build xmat
    /usr/local/pkg/afni_18/3dDeconvolve \
    -local_times \
    -force_TR 1.2 \
    -x1D_stop \
    -input ${name_img_1}" "${name_img_2} \
    -polort A \
    -float \
    -censor ${dir_stimts}/movregs_FD_mask.txt \
    -stim_times_subtract 0.6 \
    -num_stimts 1 \
    -stim_times_IM 1 ${dir_stimts_local}/${subject}_Stroop_${sessions[$session_i]}_allTrials.txt 'BLOCK(1,1)' -stim_label 1 trial \
    -ortvec ${dir_stimts}/motion_demean_${sessions[$session_i]}.1D movregs \
    -x1D ${dir_out}/X.xmat.1D \
    -xjpeg ${dir_out}/X.jpg \
    -nobucket

}

