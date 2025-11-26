'''
Compiles a series of potential exclusion criteria based on Framewise Displacement, DVARS, Rotation, and Translation.
'''
import nibabel as nib
import os
import glob
import pandas as pd
import numpy as np

database_path="/mnt/f/MiguelWorkbench/fMRIPrep_Output"
tsvs=glob.glob(os.path.join(database_path,'**','*'+'_task-ISC_desc-confounds_timeseries.tsv'), recursive=True)
main_list=[]
columns=['subject','mean_fd','mean_fd_noz',
         '80percentile_fd','80percentile_fd_noz',''
         'max_fd','max_fd_noz','percent_suprathreshold',
         'max_trans_x','max_trans_y','max_trans_z',
         'max_rot_x','max_rot_y','max_rot_z']
for tsv in tsvs:

    subject=tsv.split(os.sep)[-1][:12]
    print("Checking " + str(subject))

    #import confounds and make array
    confounds=pd.read_csv(tsv,sep='\t')
    confounds=confounds.tail(n=int(len(confounds)-5))
    confounds_par=confounds[['rot_x','rot_y','rot_z','trans_x','trans_y']]
    confounds_array = confounds_par.to_numpy()
    #calculate FD
    motion_diff = np.diff(confounds_array, axis=0, prepend=0)
    weighted_motion = np.concatenate([ 50 * np.abs(motion_diff[:, :3]), np.abs(motion_diff[:, 3:4]), np.abs(motion_diff[:, 5:6])], axis=1)
    FD = np.sum(weighted_motion, axis=1)


    #convert rot columns from radians to degrees
    confounds[['rot_x', 'rot_y', 'rot_z']] = confounds[['rot_x', 'rot_y', 'rot_z']] * (180 / np.pi)
    confounds['framewise_displacement_withouz'] = FD
    mean=confounds['framewise_displacement'].mean()
    mean_noz=confounds['framewise_displacement_withouz'].mean()
    percent_20=confounds['framewise_displacement'].quantile(0.8)
    percent20_noz=confounds['framewise_displacement_withouz'].quantile(0.8)
    max=confounds['framewise_displacement'].max()
    max_noz=confounds['framewise_displacement_withouz'].max()
    #trans/rot
    max_trans_x=confounds['trans_x'].abs().max()
    max_trans_y=confounds['trans_y'].abs().max()
    max_trans_z=confounds['trans_z'].abs().max()
    max_rot_x=confounds['rot_x'].abs().max()
    max_rot_y=confounds['rot_y'].abs().max()
    max_rot_z=confounds['rot_z'].abs().max()
    
    supra_thresh = confounds[(confounds["framewise_displacement_withouz"] > 0.25)] #threshold FD>0.5mm or DVARS>3% per Parkes et. al (2018)
    main_list.append([subject, mean, mean_noz, percent_20,percent20_noz, max, max_noz, len(supra_thresh)/len(confounds), max_trans_x, max_trans_y, max_trans_z, max_rot_x, max_rot_y, max_rot_z])


a=pd.DataFrame(main_list, columns=columns)
a.to_csv(os.path.join(os.getcwd(),'exclusion.csv'), index=False)
