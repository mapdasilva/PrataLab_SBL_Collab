#!/usr/bin/env python
import os
database='/mnt/f/MiguelWorkbench/fMRIPrep_Output'
subjects=[i for i in os.listdir(database) if i.startswith('BSCMR')]
#MAIN SUBJECT LOOP
for subj in subjects: 
	print(subj)
	#main image file
	input_path=os.path.join(database,subj,'sub-'+ subj, 'func', 'sub-'+ subj + '_task-ISC_desc-preproc_bold.nii.gz')
	smoothed_path=input_path.replace('preproc','preproc6FWHM')
	mask_path=input_path.replace('preproc_bold','brain_mask')
	#smooth if not previously done
	if os.path.isfile(smoothed_path)==False:
		smooth_cmd=f"3dBlurToFWHM -quiet -input {input_path} -prefix {smoothed_path} -mask {mask_path} -FWHM 6"
		os.system(smooth_cmd)
        
  
