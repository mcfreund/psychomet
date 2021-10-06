#!/usr/bin/env bash
## defines vars, loops over subjs, sources ../3dREMLfit_2rpm.sh
## NB: assumes wd is /data/nil-external/ccp/freund/psychomet


conda deactivate  ## make sure python3 env is not attached

## <<<<<< BEGIN USER INPUT: define vars for master script >>>>>>

cd /data/nil-external/ccp/freund/psychomet

subjects_file=in/subjects_wave12_complete_2021-09-01.txt
do_single_subj=false  ## for dev/debugging
sessions=(baseline proactive reactive)
hemis=(L R)
tasks=(Stroop)
waves=(1 2)

## <<<<<< END USER INPUT >>>>>>


## additional vars:

source code/glms/3dREMLfit_2rpm.sh  ## for 3dREMLfit function
glm_prefix="lsall_2rpm"  ## name of subdirectory in scripts; will be added to subdirectory name in out
suffix=""


if [ $do_single_subj = false ]
then
    mapfile -t subjects < $subjects_file
else
    unset subjects
    subjects=130518
    subject=130518
    wave=1
    task_i=0
    session_i=0
    task=Stroop
    session=baseline
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

            for task in ${tasks[@]}; do

                for session in ${sessions[@]}; do

                    ## define paths and names, change dir
                        
                    sess=${session:0:3}  ## get short name
                    sess=${sess^}  ## Namecase
                    dir_stimts_local=${stimts_local}${subject}/INPUT_DATA
                    dir_stimts=${stimts}${subject}/INPUT_DATA/${task}/${session}
                    dir_out=${out}${subject}/RESULTS/${task}/${session}_${glm}${suffix}

                    cd ${dir_out}

                    for hemi in ${hemis[@]}; do

                        name_img1=${img}${subject}/INPUT_DATA/${task}/${session}/lpi_scale_tfMRI_${task}${sess}1_AP_${hemi}.func.gii
                        name_img2=${img}${subject}/INPUT_DATA/${task}/${session}/lpi_scale_tfMRI_${task}${sess}2_PA_${hemi}.func.gii

                        remlfit_2rpm ${session} ${hemi} < /dev/null > ${dir_out}/runtime_3dREMLfit.log 2>&1 &

                    done

                done

            done

        cd ${wd}

        wait  ## run subjs serially, lower vars in parallel

    done

done