import os
import pandas as pd


#%%Load  data
info_df=pd.read_csv(os.path.join(os.getcwd(),'subject_list.csv')) #info_df
subj_list=info_df['subjects'].to_list()
groups=info_df['Drug'].unique()
main_df=pd.read_csv(os.path.join(os.getcwd(), 'mba_text_file.txt'), sep='\t')
new_df=pd.DataFrame()

#Subj1Group
group_df=info_df.loc[info_df["Drug"]=="DrugA"]
subj_list=group_df['subjects'].to_list()
subj_list_DrugA=[subject for subject in subj_list]

group_df=info_df.loc[info_df["Drug"]=="DrugB"]
subj_list=group_df['subjects'].to_list()
subj_list_DrugB=[subject for subject in subj_list]

mapping_dict_DrugA={key: 'A' for key in subj_list_DrugA}
mapping_dict_DrugB={key: 'B' for key in subj_list_DrugB}
mapping_dict = {**mapping_dict_DrugA, **mapping_dict_DrugB}

main_df['Subj1Group']=main_df['Subj1'].map(mapping_dict)
main_df['Subj2Group']=main_df['Subj2'].map(mapping_dict)
main_df['GroupStr'] = main_df.apply(lambda row: row['Subj1Group'] + '_' + row['Subj2Group'], axis=1) 
mapping_dict_order={'A_B':'A_B', 'B_A':'A_B','A_A':'A_A','B_B':'B_B'}
main_df['GroupStr']=main_df['GroupStr'].map(mapping_dict_order)

def categorize(row):
    if row['Video'].startswith('v5'):
        return "neutral"
    elif row['Video'].startswith('v2'): 
        return "erotic"
    elif row['Video'].startswith('v1'): 
        return "horror"
        

#Categorize Socialness
main_df['VideoType'] = main_df.apply(lambda row: categorize(row),axis=1)

#Save

within_df=main_df.loc[main_df['GroupStr'].isin(['A_A','B_B'])]
#between_df=main_df.loc[main_df['GroupStr'].isin(['HC_HC','HC_SCZA','HC_SCZB'])]

within_df.to_csv(os.path.join(os.getcwd(), 'data_yeo7_modelwithin.txt'), sep='\t', index=False)
#between_df.to_csv(os.path.join(os.getcwd(), 'data_hammer_modelB.txt'), sep='\t', index=False)

#main_df.to_csv(os.path.join(os.getcwd(), 'data_aal_modelFull.txt'), sep='\t', index=False)


