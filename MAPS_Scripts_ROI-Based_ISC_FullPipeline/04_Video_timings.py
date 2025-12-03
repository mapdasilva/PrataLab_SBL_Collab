import os
import pandas as pd 
import glob

files=glob.glob(os.path.join(os.getcwd(),"Behavioral",'**',"RESULTS_VIDEO.txt"), recursive=True)

for file in files:

    subj_df=pd.read_csv(file, sep='\t')
    #Adjust timings
    subj_df['tITI2shown_ADJ']=subj_df.apply(lambda row: round(row['tITI2shown']-row['TimeFirst_ScannerPulse'])/1000, axis=1)
    subj_df['tVIDEO_START_ADJ']=subj_df.apply(lambda row: round(row['tVIDEO_START']-row['TimeFirst_ScannerPulse'])/1000, axis=1)
    subj_df['tITI1shown_ADJ']=subj_df.apply(lambda row: round(row['tITI1shown']-row['TimeFirst_ScannerPulse'])/1000, axis=1)

    for index, row in subj_df.iterrows():
        #Select TRs for each video 
        start_idx=round(row['tVIDEO_START_ADJ']/1.85)+4 
        #gets the first TR for which its middle point overlaps with the video, then adds 3 for hemodynamic \
        #compensation and 1 extra for indexing starting at 1 (first TR is TR number 1)     
        subj_df.at[index,'TR_first']=start_idx
        subj_df.at[index,'TR_last']=start_idx+18

        # warning for skipped video and dropped frames
        if row['DROPPED_FRAMES']!=0:
            print(f"Dropped {row['DROPPED_FRAMES']} in video {index} for subject {os.path.dirname(file).split(os.sep)[-1]}")
        if row['Skipped']!=False:
            print(f"Skipped video {index} for subject {os.path.dirname(file).split(os.sep)[-1]}")
    version=subj_df.at[0,"VERSION"].lower()
    subj_df["video_id"]=subj_df[f"video_id_v{version}"]
    subj_df["condition"]=subj_df[f"condition_v{version}"]
    subj_df=subj_df.drop(columns=["video_id_va","video_id_vb","condition_va","condition_vb"])

#Miguel Modifications
# ---- NEW PART: add subject ID to filename ----
    subj_folder = os.path.dirname(file)          # .../Behavioral/B012ISC
    folder_name = os.path.basename(subj_folder)  # "B012ISC"

    # assuming pattern B<subject>ISC, e.g. "B012ISC" -> "012"
    subject_id = folder_name.replace("B", "").replace("ISC", "")
    subject_id = subject_id.zfill(3)

    out_name = f"ISC_VideoTimings_BSCMR{subject_id}.csv"
    out_path = os.path.join(subj_folder, out_name)

    subj_df.to_csv(out_path, index=False)
    print(f"Saved {out_name} for subject {folder_name}")

#Previous Script from Vasco
    #subj_df.to_csv(os.path.join(os.path.dirname(file), "RESULTS_VIDEO_BSCMR{subj}.csv" ), index=False)
    #print(os.path.dirname(file).split(os.sep)[-1])
    #subj_df.to_csv(os.path.join(os.path.dirname(file), "LEARN_RESULTS.csv" ), index=False)
    #print(os.path.dirname(file).split(os.sep)[-1])
