import os
import nibabel as nib
import pandas as pd
import glob
from pqdm.processes import pqdm
import subprocess
def runproc(string):
    proc= subprocess.Popen(string,shell=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)     
    stdout, stderr= proc.communicate()  
    
    return proc.returncode, stdout, stderr

database_path="/mnt/f/MiguelWorkbench/fMRIPrep_Output"
niis=glob.glob(os.path.join(database_path, "**","denoised_func_data_nonaggr_MNI.nii.gz"), recursive=True)

#Loop subjects
window_len=19
dummy_scans=5
def split_bold(nii):
    subject=nii.split(os.sep)[-3].replace("BSCMR","")
    print(subject)
    behavioral_csv=os.path.join("/mnt/f/MiguelWorkbench/Behavioral",f"B{subject}ISC","LEARN_RESULTS.csv")
    times_df=pd.read_csv(behavioral_csv)
    for _, row in times_df.iterrows():
        video_id=row["video_id"][:-4]
        start_idx=int(row["TR_first"])-dummy_scans
        output=nii[:-7]+'_video-'+str(video_id)+'.nii.gz'
        fsl_cmd = f'fslroi {nii} {output} {start_idx} {window_len}'
        if os.path.isfile(output)==False:
            code, out, err = runproc(fsl_cmd)
            if err!=0 and err!=None:
                print(f"Error in subject {subject} for video {video_id}")

    return None

result = pqdm(niis, split_bold, n_jobs=20)

