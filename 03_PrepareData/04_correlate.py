import os
import glob
import warnings
import pandas as pd
import subprocess
import itertools
from pqdm.processes import pqdm

#cmd function
def runproc(string):
    proc= subprocess.Popen(string,shell=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)     
    stdout, stderr= proc.communicate()  
    
    return proc.returncode, stdout, stderr



def correlate_pair(pair):
    print(pair)
    for video_id in video_ids:
        file1 = os.path.join(os.getcwd(), "roistats", pair[0] + "_v" + str(video_id) + ".1D")
        file2 = os.path.join(os.getcwd(), "roistats", pair[1] + "_v" + str(video_id) + ".1D")
        if os.path.isfile(file1) and os.path.isfile(file2):
            output_file = os.path.join(os.getcwd(), "pairwise_correlations", pair[0] + "_" + pair[1] + "_v" + str(video_id) + ".1D")
            correlate_cmd = "3dTcorrelate -polort -1 " + "-prefix " + output_file + " " + file1 + "\\' " + file2 + "\\'" 
            code, out, err = runproc(correlate_cmd)
            if err!=0 and err!=None:
                print("Error " + pair + video_id)
    return None

#%%Load fMRI data
info_df=pd.read_csv('subject_list.csv') #info_df
database_path="/mnt/f/MiguelWorkbench/fMRIPrep_Output"
subj_list=info_df['subjects']


#%% Set up experiment metadata
video_ids=['1001','1002','1006','1007',
           '2000','2002','2003','2009',
           '5000','5001','5002','5008']

# Create output folder if it doesn't exist
os.makedirs("pairwise_correlations", exist_ok=True)

#Generate all subject pairs
subj_pairs=itertools.combinations(subj_list,2)
result = pqdm(tuple(subj_pairs), correlate_pair, n_jobs=16)

