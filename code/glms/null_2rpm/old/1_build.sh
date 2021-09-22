#!/usr/bin/env bash

#exec env --ignore-environment /bin/bash  ## new bash env

## <<<<<< BEGIN USER INPUT: define vars for master script >>>>>>

cd /data/nil-external/ccp/freund/psychomet

subjects_file=in/subjects_crosswave_complete_2021-09-01.txt
do_single_subj=false  ## for dev/debugging

sessions=(baseline proactive reactive)
tasks=(Axcpt Cuedts Stroop Stern)
#sessions=baseline
#tasks=Stroop
waves=(1 2)

## <<<<<< END USER INPUT >>>>>>


## additional vars:

glm_prefix="null_2rpm"  ## name of subdirectory in scripts; will be added to subdirectory name in out
hemis=(L R)

if [ $do_single_subj = false ]
then
    mapfile -t subjects < $subjects_file
else
    subjects=DMCC8033964
fi

## run:

for wave in ${waves[@]}; do

    echo "running wave "${wave}

    glm=$glm_prefix"_wave"$wave

    if [ ${wave} = 1 ]
    then
        wave_dir=HCP_SUBJECTS_BACKUPS
    elif [ ${wave} = 2 ]
    then
        wave_dir=DMCC_Phase$wave
    fi

    # paths:

    stimts=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/$wave_dir/fMRIPrep_AFNI_ANALYSIS/  ## JUST NEED MOVEMENT REGS
    out=/data/nil-external/ccp/freund/psychomet/out/glms/
    img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/$wave_dir/fMRIPrep_AFNI_ANALYSIS/
    scripts=/data/nil-external/ccp/freund/psychomet/code/glms/

    for subject in ${subjects[@]}; do

        echo ${subject}
        
        source code/glms/${glm_prefix}/3dDeconvolve.sh

    done

done
