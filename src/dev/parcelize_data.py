import os
import sys
import shlex
import argparse
import cPickle
import numpy as np
import nibabel as nib
import multiprocessing as mp
from subprocess import Popen, PIPE

np.set_printoptions(precision=2, suppress=True)
np.set_printoptions(threshold=sys.maxsize)


### CREATE LIST OF SUBJECTS ###
parser = argparse.ArgumentParser(description='Process input files')
parser.add_argument('subjnums', type=argparse.FileType('r'))
args = parser.parse_args()

subjnums = []

for line in args.subjnums:
        subj = line.rstrip("\n\r")
        subjnums.append(subj)

### CREATE DATA MASKS ###
data_masks = np.empty([400, 64984])
mask_sizes = np.empty([400])

ci1_img = nib.cifti2.load("/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/Schaefer2018_Parcellations/HCP/fslr32k/cifti/Schaefer2018_400Parcels_7Networks_order.dscalar.nii")
atlas = ci1_img.dataobj.get_unscaled()[0]

for parcelation in range(1, 401):
        print("Parcelation " + str(parcelation))
        mask_count = 0
        for vertex in range(0, 64984):
                if atlas[vertex] == parcelation:
                        data_masks[parcelation - 1][vertex] = 1
                        mask_count += 1
                else:
                        data_masks[parcelation - 1][vertex] = 0

        mask_sizes[parcelation - 1] = mask_count

### FIND INDEX OF LABELS AND PUT INTO A LIST ###
labels = ["HI_LO_conf#0_Coef", "Acue_Bcue#0_Coef", "error_correct#0_Coef"]
label_indicies = []

for label in labels:
        cmd_3dinfo = "3dinfo -label2index " + label + " /data/nil-external/ccp/witzermanm/AFNI_ANALYSIS_SUBSUBJECT/RESULTS_RUNWISE/" + subjnums[0] + "/Axcpt/baseline/baseline_Cues_EVENTS_censored_run1/stats_" + subjnums[0] + "_L.func.gii"
        pipe_cmd = shlex.split(cmd_3dinfo)
        p = Popen(pipe_cmd, stdout = PIPE)
        darrayIndex = int(p.communicate()[0])
        label_indicies.append(darrayIndex)


### CREATE FULL DATASET ###
#takes a number as an input to represent a parcellation (0 - 399) and uses the corresponding data mask to create a dataset for all subjects
def create_dset(parcelation):
        parcelized_dset = np.empty([3, 118, int(mask_sizes[parcelation])])

        label_iterator = 0
        for label_index in label_indicies:
                for run_iterator in range(1, 3):
                        subj_iterator = 0
                        for subj in subjnums:
                                print("Label " + str(label_iterator + 1) + " | Run " + str(run_iterator) + " | Subj " + str((subj_iterator + 1) + (59 * (run_iterator - 1))))

                                ### GET SUBJET DATA FROM GIFTI FILE ###
                                gifti_path_L = "/data/nil-external/ccp/witzermanm/AFNI_ANALYSIS_SUBSUBJECT/RESULTS_RUNWISE/" + subj + "/Axcpt/baseline/baseline_Cues_EVENTS_censored_run" + str(run_iterator) + "/stats_" + subj + "_L.func.gii"
                                gifti_path_R = "/data/nil-external/ccp/witzermanm/AFNI_ANALYSIS_SUBSUBJECT/RESULTS_RUNWISE/" + subj + "/Axcpt/baseline/baseline_Cues_EVENTS_censored_run" + str(run_iterator) + "/stats_" + subj + "_R.func.gii"

                                if not os.path.exists(gifti_path_L):
                                        subj_iterator += 1
                                        continue

                                gifti_img_L = nib.gifti.giftiio.read(gifti_path_L)
                                gifti_img_R = nib.gifti.giftiio.read(gifti_path_R)


                                combined_hemispheres = np.concatenate((gifti_img_L.darrays[label_index].data, gifti_img_R.darrays[label_index].data))

                                ### RUN DATA THROUGH MASK ###
                                parcelation_iterator = 0
                                for vertex in range(0, 64984):
                                        ### IF VERTEX VALUE IN MASK IS 1, ADD THAT VERTEX'S DATA TO DATASET ###
                                        if data_masks[parcelation][vertex] == 1:
                                                parcelized_dset[label_iterator][subj_iterator + (59 * (run_iterator - 1))][parcelation_iterator] = combined_hemispheres[vertex]
                                                parcelation_iterator += 1

                                del gifti_img_L
                                del gifti_img_R

                                subj_iterator += 1

                        print(parcelized_dset[label_iterator][run_iterator * 58])

                label_iterator += 1

        print("Adding parcelation " + str(parcelation + 1) + " data to full_dset...")
        full_dset[parcelation] = parcelized_dset


#call create_dset function many times using multiprocessing
if __name__ == '__main__':
        full_dset = [None] * 400
        pool = mp.Pool(processes = 32)
        pool.map(create_dset, range(0, 5))

        print(full_dset[0])
        print("Pickling...")
        pkl_file = open('full_dset.pkl', 'wb')
        cPickle.dump(full_dset, pkl_file)
