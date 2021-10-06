#!/usr/bin/env bash

#exec env --ignore-environment /bin/bash  ## new bash env

## <<<<<< BEGIN USER INPUT: define vars for master script >>>>>>

cd /data/nil-external/ccp/freund/psychomet

subjects_file=in/subjects_wave12_complete_2021-09-01.txt
do_single_subj=false  ## for dev/debugging

sessions=(baseline proactive reactive)
tasks=(Stroop)
waves=(1 2)

## <<<<<< END USER INPUT >>>>>>


## additional vars:

glm_prefix="lsall_2rpm"  ## name of subdirectory in scripts; will be added to subdirectory name in out
source code/glms/${glm_prefix}/3dDeconvolve.sh  ## for 3dDeconvolve function
hemis=(L R)

if [ $do_single_subj = false ]
then
    mapfile -t subjects < $subjects_file    
else  ## for dev
    unset subjects
    subjects=130518
    subject=130518
    wave=1
    task_i=0
    session_i=0
fi

## run:

wd=$(pwd)

for wave in ${waves[@]}; do

    echo "running wave "${wave}

    glm=$glm_prefix"_wave"$wave

    if [ ${wave} = 1 ]
    then
        wave_dir=HCP_SUBJECTS_BACKUPS
    else  ## wave 2 or greater
        wave_dir=DMCC_Phase$((wave+1))  ## add one to get phase
    fi

    # paths:

    stimts_local=/data/nil-external/ccp/freund/psychomet/out/glms/  ## WHERE EVENT STIMTIMES ARE LOCATED
    stimts=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/$wave_dir/fMRIPrep_AFNI_ANALYSIS/  ## JUST NEED MOVEMENT REGS
    out=/data/nil-external/ccp/freund/psychomet/out/glms/
    img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/$wave_dir/fMRIPrep_AFNI_ANALYSIS/
    scripts=/data/nil-external/ccp/freund/psychomet/code/glms/

    for subject in ${subjects[@]}; do

        echo ${subject}

        for task_i in ${!tasks[@]}; do

            for session_i in ${!sessions[@]}; do

                ## define paths and names
                sess=${sessions[$session_i]:0:3}  ## get short name
                sess=${sess^}  ## Namecase
                dir_stimts_local=${stimts_local}${subject}/INPUT_DATA
                dir_stimts=${stimts}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}
                dir_out=${out}${subject}/RESULTS/${tasks[$task_i]}/${sessions[$session_i]}_${glm}
                name_img_1=${img}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}/lpi_scale_tfMRI_${tasks[$task_i]}${sess}1_AP_L.func.gii  ## either hemi OK
                name_img_2=${img}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}/lpi_scale_tfMRI_${tasks[$task_i]}${sess}2_PA_L.func.gii

                ## make result dir and cd into it
                mkdir -p ${dir_out}
                cd ${dir_out}  ## IS THIS NECESSARY?

                ## build xmat

                deconvolve_2rpm_lsall_stroop < /dev/null > ${dir_out}/runtime_3dDeconvolve.log 2>&1 &  ## THIS MAY NOT WORK!!! NEEDS TESTING

                cd ${wd}  ## back to original dir

           	done

        done

        wait  ## run subjs serially, but each task and session in parallel

    done

done
