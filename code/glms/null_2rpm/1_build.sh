#!/usr/bin/env bash


## <<<<<< BEGIN USER INPUT: define vars for master script >>>>>>

cd /data/nil-external/ccp/freund/psychomet

subjects_file=in/subjects_complete_2021-09-01.txt
do_single_subj=false  ## for dev/debugging
sessions=(baseline proactive reactive)
hemis=(L R)
tasks=(Axcpt Cuedts Stroop Stern)
wave=2

## <<<<<< END USER INPUT >>>>>>


## additional vars:

glm="null_2rpm"  ## name of subdirectory in scripts; will be added to subdirectory name in out


if [ ${wave} = 2 ]
then
    wave_dir=HCP_SUBJECTS_BACKUPS
else
    wave_dir=DMCC_Phase$wave
fi

if [ $do_single_subj = true ]
then
    subject='132017'
    session_i=0
    task_i=0
    run_i=0
else
    mapfile -t subjects < $subjects_file
fi


# paths:

stimts=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/$wave_dir/fMRIPrep_AFNI_ANALYSIS/  ## JUST NEED MOVEMENT REGS
out=/data/nil-external/ccp/freund/psychomet/out/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/$wave_dir/fMRIPrep_AFNI_ANALYSIS/
scripts=/data/nil-external/ccp/freund/psychomet/code/glms/


## run:

for subject in ${subjects[@]}; do

	echo ${subject}
	
	source code/glms/${glm}/3dDeconvolve.sh

done
