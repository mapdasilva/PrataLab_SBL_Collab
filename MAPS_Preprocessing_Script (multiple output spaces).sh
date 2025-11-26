#!/bin/bash
#Runs fMRIPrep preprocessing on the subjcts listed in the subjects_list.txt file.
#Dataset must be bids compliant. 

#Paths to set
SUBJECTS_FOLDER=/mnt/f/MiguelWorkbench/BSCMR_BIDS_ISC
FS_LICENSE=/opt/freesurfer/license.txt
OUTPUT_PATH=/mnt/f/MiguelWorkbench/Output2
FS_SUBJECTS_DIR=/mnt/f/MiguelWorkbench/freesurfer

#Loop over subjects
for line in $(cat BSC_subj_list4.txt); 
do	
	line="${line%$'\r'}"
	docker run -v $SUBJECTS_FOLDER/:/SUBJECTS_BIDS \
	-v $FS_LICENSE:/license.txt \
	-v $FS_SUBJECTS_DIR:/freesurfer \
	-v $OUTPUT_PATH/$line:/FMRIPREP_OUTPUTS nipreps/fmriprep:latest /SUBJECTS_BIDS /FMRIPREP_OUTPUTS participant \
	--participant-label $line \
	--skip-bids-validation \
	--task-id ISC \
	--nthreads 16 \
	--mem 48000 \
	--output-spaces func MNI152NLin6Asym:res-2 MNI152NLin2009cAsym\
	--bold2anat-dof 12 \
	--no-submm-recon \
	--fallback-total-readout-time estimated \
	--fs-license-file /license.txt \
	--use-syn-sdc \
	--fs-subjects-dir /freesurfer \
	
done
