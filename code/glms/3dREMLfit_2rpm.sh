#!/usr/bin/env bash

# assumes subject ID is defined as variable $subject
# https://stackoverflow.com/questions/17621798/linux-process-in-background-stopped-in-jobs
# for writing stdin/stdout to file, for background running


wd=$(pwd)

function remlfit {
	
	local session=$1
	local hemi=$2
	
	## define paths and names
	sess=${session:0:3}  ## get short name
	sess=${sess^}  ## Namecase
	dir_stimts=${stimts}${subject}/INPUT_DATA/${task}/${session}
	dir_out=${out}${subject}/RESULTS/${task}/${session}_${glm}${suffix}
	name_img1=${img}${subject}/INPUT_DATA/${task}/${session}/lpi_scale_tfMRI_${task}${sess}1_AP_${hemi}.func.gii
	name_img2=${img}${subject}/INPUT_DATA/${task}/${session}/lpi_scale_tfMRI_${task}${sess}2_PA_${hemi}.func.gii

	cd ${dir_out}

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


for task in ${tasks[@]}; do

	for session in ${sessions[@]}; do

		for hemi in ${hemis[@]}; do

			logpath=${out}${subject}/RESULTS/${task}/${session}_${glm}${suffix}

			#echo ${logpath}

			remlfit ${session} ${hemi} < /dev/null > ${logpath}/runtime.log 2>&1 &

		done

	done

done




cd ${wd}
