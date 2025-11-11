import os
import pandas as pd
import numpy as np

#%%Load  data
within_df=pd.read_csv(os.path.join(os.getcwd(), 'data_yeo7_modelwithin.txt'), sep='\t')
within_df=within_df[['ROI','Subj1','Subj2','Y','GroupStr','VideoType']]


#Subj1Group
#Make B/D
within_df['Yz']=np.arctanh(within_df['Y'])
within_df_new=within_df.groupby(['ROI','Subj1','Subj2', 'GroupStr','VideoType'], as_index=False).mean()
within_df_new['Y']=np.tanh(within_df_new['Yz'])
within_df_new.to_csv(os.path.join(os.getcwd(), 'data_yeo7_avg_modelwithin.txt'), sep='\t', index=False)


