import os
import glob
import warnings
import pandas as pd
import subprocess
from pqdm.processes import pqdm

#cmd function
def runproc(string):

    proc= subprocess.Popen(string,shell=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)     
    stdout, stderr= proc.communicate()  
    
    return proc.returncode, stdout, stderr


#%%Load fMRI data
#info_df=pd.read_csv(os.path.join(os.getcwd(),'PatientDetails.csv')) #info_df
database_path="/mnt/f/MiguelWorkbench/fMRIPrep_Output"
subj_list=os.listdir(database_path)

#%%Load atlas
atlas_path=os.path.join(os.getcwd(),'Yeo7resampled.nii')

#%% Set up experiment metadata
video_ids=['1001','1002','1006','1007',
           '2000','2002','2003','2009',
           '5000','5001','5002','5008']

def parcellate_subject(subj):
    if os.path.isdir(os.path.join(database_path,subj,'aroma')):
        os.makedirs("roistats", exist_ok=True)  # Creates roistats folder
        for video_id in video_ids:
            nii=os.path.join(database_path,subj,'**',f'denoised_func_data_nonaggr_MNI_video-{video_id}.nii.gz')
            output_file = os.path.join(os.getcwd(), "roistats", nii.split(os.sep)[-3] + "_v" + str(video_id) + ".1D")
            roi_cmd = "3dROIstats -mask " + atlas_path + " -1Dformat -nobriklab " + nii + " > " + output_file
            code, out, err = runproc(roi_cmd)
            if err!=0 and err!=None:
                print(f"Error in subject {subj} for video {video_id}")
    return None

result = pqdm(subj_list, parcellate_subject, n_jobs=16)



