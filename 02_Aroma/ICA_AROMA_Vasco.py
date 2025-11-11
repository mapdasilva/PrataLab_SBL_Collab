#!/usr/bin/env python

# Import required modules
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
import os
import sys
import ICA_AROMA_functions as aromafunc
import pandas as pd
import nibabel as nib
import json
import numpy as np
from pqdm.processes import pqdm

def runAROMA(database_path, bids_path, subject, task, dummy_scans, FWHM, outDir, mc, confounds_path):
    # Change to script directory
    cwd = os.path.realpath(os.path.curdir)
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(scriptDir)
    # Define the FSL-bin directory
    fslDir = os.path.join(os.environ["FSLDIR"], 'bin', '')

    #-------------------------------------------- PARSER --------------------------------------------#
    #Set filepaths
    func_base=os.path.join(database_path, subject, f'sub-{subject}','func',f'sub-{subject}_task-{task}')
    inFile= f'{func_base}_desc-preproc{FWHM}FWHMminus{dummy_scans}_bold.nii.gz'
    inFile_MNI= f'{func_base}_space-MNI152NLin6Asym_res-2_preproc{FWHM}FWHMminus{dummy_scans}_bold.nii.gz'
    uncut= f'{func_base}_desc-preproc{FWHM}FWHM_bold.nii.gz'
    uncut_MNI=  f'{func_base}_space-MNI152NLin6Asym_res-2_desc-preproc{FWHM}FWHM_bold.nii.gz'
    mask=f'{func_base}_desc-brain_mask.nii.gz'
    mask_MNI=f'{func_base}_space-MNI152NLin6Asym_res-2_desc-brain_mask.nii.gz'
    # Pre flight check
    if os.path.isfile(inFile)==False:
        if os.path.isfile(uncut)==True:
            fsl_cmd=f'fslroi {uncut} {inFile} {dummy_scans} -1'
            os.system(fsl_cmd)
        else:
            print(f"Missing {uncut} ")
            return None  
    if os.path.isfile(inFile_MNI)==False:
        if os.path.isfile(uncut_MNI)==True:
            fsl_cmd=f'fslroi {uncut_MNI} {inFile_MNI} {dummy_scans} -1'
            os.system(fsl_cmd)
        else:
            print(f"Missing {uncut_MNI} ")
            return None  

    if os.path.isdir(outDir):
        print(f'Output directory {outDir} already exists.')
    else:
        os.makedirs(outDir)

    json_file=os.path.join(bids_path, f'sub-{subject}','func',f'sub-{subject}_task-{task}_bold.json')
    with open(json_file, 'r') as file:
        data = json.load(file)
        #to do get from func json
        bold_metadata={'RepetitionTime': data['RepetitionTime'],
                    'PhaseEncodingDirection': data['PhaseEncodingDirection'],
                    'TotalReadoutTime': data['EstimatedTotalReadoutTime']
            }
    file.close()
    
    #---------------------------------------- Run ICA-AROMA ----------------------------------------#
    ####apply fieldmap to inFile before aroma
    
    # STEP 1 - MELODIC
    melIC = os.path.join(outDir, 'melodic_IC_thr.nii.gz')
    if os.path.isfile(melIC)==False:
        print('Step 1) MELODIC')
        aromafunc.runICA(fslDir, inFile, outDir, mask, dim=0, TR=bold_metadata['RepetitionTime'])
  
    # STEP 2 - Classification
    print('Step 2) Automatic classification of the components')
    print('  - registering MELODIC to MNI') 
    melIC_MNI = os.path.join(outDir, 'melodic_IC_thr_MNI2mm.nii.gz')
    melmix = os.path.join(outDir, 'melodic.ica', 'melodic_mix')
    melFTmix = os.path.join(outDir, 'melodic.ica', 'melodic_FTmix')
    if os.path.isfile(melIC_MNI)==False:
        aromafunc.register2MNI(database_path, subject, melIC, bold_metadata, melIC_MNI, omp_nthreads=16, mem_gb=20)

        print('  - extracting the CSF & Edge fraction features')
        edgeFract, csfFract = aromafunc.feature_spatial(fslDir, outDir, scriptDir, melIC_MNI)

        print('  - extracting the Maximum RP correlation feature')
        maxRPcorr = aromafunc.feature_time_series(melmix, mc)

        print('  - extracting the High-frequency content feature')
        HFC = aromafunc.feature_frequency(melFTmix, TR=data['RepetitionTime'])

        print('  - classification')
        motionICs = aromafunc.classification(outDir, maxRPcorr, edgeFract, HFC, csfFract)
        
    # STEP 3 - Data Denoising
    print('Step 3) Data denoising')
    motionIDX=np.loadtxt(os.path.join(outDir,'classified_motion_ICs.txt'),delimiter=',', dtype=int)
    nonaggr=os.path.join(outDir,'denoised_func_data_nonaggr_MNI.nii.gz')
    if os.path.isfile(nonaggr)==False:
        allnoise_path=aromafunc.create_confounds(melmix, motionIDX, confounds_path)
        aromafunc.denoising_afni(inFile_MNI, mask_MNI, nonaggr, allnoise_path, highpass=0.006, lowpass=0.15, polort=2)
    
    # Revert to old directory
    os.chdir(cwd)
    print('\n----------------------------------- Finished -----------------------------------\n')


task='ISC'
dummy_scans=5
FWHM=6
database='/mnt/f/MiguelWorkbench/fMRIPrep_Output' #, folder
bids_path='/mnt/f/MiguelWorkbench/BSCMR_BIDS_ISC' 
subjects=[i for i in os.listdir(database) if i.startswith('BSCMR')]
subj_list=pd.read_csv('subject_list.csv')
failed_subjects=[]
#MAIN SUBJECT LOOP
for subj in subj_list['subjects']:
    try:
        #directories
        func_base=os.path.join(database,subj,f'sub-{subj}', 'func', f'sub-{subj}_task-{task}')
        outDir=os.path.join(database,subj,'aroma')

        #make mc file
        mc=func_base+ '_desc-mc_confounds.txt'
        if os.path.isfile(mc)==False:
            confounds_df=pd.read_csv(func_base +  '_desc-confounds_timeseries.tsv', sep='\t')
            mc_df=confounds_df[['trans_x','trans_y','trans_z','rot_x','rot_y','rot_z']]
            mc_df=mc_df.tail(n=int(len(mc_df)-dummy_scans)) #remove first n timepoints
            mc_df.to_csv(mc,sep=' ', header=False, index=False)

        #make csfwm file
        csf_wm=func_base+ '_desc-csfwm_confounds.txt'
        if os.path.isfile(csf_wm)==False:
            confounds_df=pd.read_csv(func_base +  '_desc-confounds_timeseries.tsv', sep='\t')
            csf_wm_df=confounds_df[['csf','white_matter']]
            csf_wm_df=csf_wm_df.tail(n=int(len(csf_wm_df)-dummy_scans)) #remove first n timepoints
            csf_wm_df.to_csv(csf_wm,sep=' ', header=False, index=False)
        
        runAROMA(database, bids_path, subj, task, dummy_scans, FWHM, outDir, mc, csf_wm)
    except:
        print(f"{subj}  failed AROMA")
        failed_subjects.append(subj)
        




    
        
 
    #ACTUALLY RUN AROMA
   




