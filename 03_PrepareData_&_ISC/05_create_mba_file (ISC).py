import os
import pandas as pd
import numpy as np
'''
Generate text file according to modified spec from:

https://afni.nimh.nih.gov/pub/dist/doc/program_help/MBA.html

 A data table in pure text format is needed as input for an MBA script. The
 data table should contain at least 4 columns that specify the information
 about subjects, region pairs and the response variable values with the
 following fixed header. The header labels are case-sensitive, and their order
 does not matter.

 ROI    Subj1   Subj2   Y       Video      
 Amyg   NGN001  NGN002     0.2643  v2003
 PFC    NGN001  NGN003     0.3762  v2004

Avoid using pure numbers to code the labels for categorical variables.
'''


pairwise_folder=os.path.join(os.getcwd(),'pairwise_correlations')
columns=['ROI','Subj1','Subj2','Y', 'Video']
main_df=pd.DataFrame(columns=columns)

#roi info
rois=pd.read_csv(os.path.join(os.getcwd(),'Yeo7.txt'), sep='\t')
i=0
for file in os.listdir(pairwise_folder):
    i+=1
    Subj1=file.split('_')[0]
    Subj2=file.split('_')[1]
    Video=file.split('_')[2][-8:-3]
    file_path=os.path.join(pairwise_folder, file)
    #load correlation values
    pair_arr=np.loadtxt(file_path)
    pair_df=pd.DataFrame(pair_arr, columns=['Y'])
    #Add info
    pair_df=pair_df.assign(Subj1=Subj1)
    pair_df=pair_df.assign(Subj2=Subj2)
    pair_df=pair_df.assign(Video=Video)
    pair_df['ROI']=rois['Name']
    #save

    main_df=pd.concat([main_df, pair_df], ignore_index=True)
    print(i)

main_df.to_csv(os.path.join(os.getcwd(), 'mba_text_file.txt'), sep='\t', index=False)