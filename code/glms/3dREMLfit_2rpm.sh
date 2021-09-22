#!/usr/bin/env bash

# assumes subject ID is defined as variable $subject
# https://stackoverflow.com/questions/17621798/linux-process-in-background-stopped-in-jobs
# for writing stdin/stdout to file, for background running


function remlfit_2rpm {
	
	local session=$1
	local hemi=$2
	
	/usr/local/pkg/afni_18/3dREMLfit \
	-matrix ${dir_out}/X.xmat.1D \
	-input ${name_img1}" "${name_img2} \
	-Rvar ${dir_out}/stats_var_${subject}_${hemi}_REML.func.gii \
	-Rbuck ${dir_out}/STATS_${subject}_${hemi}_REML.func.gii \
	-rwherr ${dir_out}/wherr_${subject}_${hemi}_REML.func.gii \
	-rerrts ${dir_out}/errts_${subject}_${hemi}_REML.func.gii \
	-GOFORIT \
	-fout \
	-tout \
	-nobout \
	-noFDR \
	-verb

}


