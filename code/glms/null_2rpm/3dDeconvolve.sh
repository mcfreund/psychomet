#!/usr/bin/env bash

wd=$(pwd)

for task_i in ${!tasks[@]}; do


	for session_i in ${!sessions[@]}; do

			## define paths and names
			sess=${sessions[$session_i]:0:3}  ## get short name
			sess=${sess^}  ## Namecase
			dir_stimts=${stimts}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}
			dir_out=${out}${subject}/RESULTS/${tasks[$task_i]}/${sessions[$session_i]}_${glm}
			name_img_1=${img}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}/lpi_scale_tfMRI_${tasks[$task_i]}${sess}1_AP_L.func.gii  ## either hemi OK
			name_img_2=${img}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}/lpi_scale_tfMRI_${tasks[$task_i]}${sess}2_PA_L.func.gii

			## make result dir
			mkdir -p ${dir_out}
			cd ${dir_out}
		
			## build xmat
			/usr/local/pkg/afni_18/3dDeconvolve \
			-local_times \
			-force_TR 1.2 \
			-x1D_stop \
			-input ${name_img_1}" "${name_img_2} \
			-polort A \
			-float \
			-censor ${dir_stimts}/censor_list.1D \
			-num_stimts 0 \
			-ortvec ${dir_stimts}/motion_demean_${sessions[$session_i]}.1D movregs \
			-x1D ${dir_out}/X.xmat.1D \
			-xjpeg ${dir_out}/X.jpg \
			-nobucket

	done



done

cd ${wd}
