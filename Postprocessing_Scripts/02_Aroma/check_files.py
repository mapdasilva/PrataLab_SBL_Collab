import os


database='/mnt/f/MiguelWorkbench/fMRIPrep_Output' #, folder
subjects=[i for i in os.listdir(database) if i.startswith('BSCMR')]
#MAIN SUBJECT LOOP
for subj in subjects: 
    inFile= os.path.join(database, subj, f'sub-{subj}', 'func', f'sub-{subj}_task-ISC_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz')
    #inFile= os.path.join(database, subj, f'sub-{subj}', 'fmap', f'sub-{subj}_fmapid-auto00000_desc-preproc_fieldmap.json')

    if os.path.isfile(inFile)==False:
        print(f"Missing {inFile}")

