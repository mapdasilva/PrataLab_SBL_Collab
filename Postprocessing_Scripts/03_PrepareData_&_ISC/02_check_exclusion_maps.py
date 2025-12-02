#Compiles a series of potential exclusion criteria based on Framewise Displacement, (DVARS??), Rotation, and Translation.

import os
import glob
import pandas as pd
import numpy as np

database_path = "/mnt/f/MiguelWorkbench/fMRIPrep_Output"

# Find all confounds files for the ISC task
tsvs = glob.glob(
    os.path.join(database_path, "**", "*_task-ISC_desc-confounds_timeseries.tsv"),
    recursive=True
)

rows = []
columns = [
    "subject",
    "mean_fd",
    "mean_fd_noz",
    "p80_fd",
    "p80_fd_noz",
    "max_fd",
    "max_fd_noz",
    "percent_suprathreshold",
    "max_trans_x",
    "max_trans_y",
    "max_trans_z",
    "max_rot_x",
    "max_rot_y",
    "max_rot_z",
]

FD_THRESHOLD = 0.5  # mm, consistent with your Parkes et al. comment - or 0.25 as was in Vasco's Script - More conservative 

for tsv in tsvs:

    # ----- subject ID -----
    fname = os.path.basename(tsv)
    # Keep your original approach if you know it's always 12 chars:
    subject = fname[:12]
    print(f"Checking {subject}")

    # ----- load confounds -----
    confounds = pd.read_csv(tsv, sep="\t")

    # Drop first 5 volumes (e.g. dummies) from the confounds as you intended
    confounds = confounds.iloc[5:].reset_index(drop=True)

    # ----- build motion parameter array (radians + mm) -----
    # include ALL 6 DOF
    rot_cols = ["rot_x", "rot_y", "rot_z"]
    trans_cols = ["trans_x", "trans_y", "trans_z"]

    confounds_par = confounds[rot_cols + trans_cols].to_numpy()

    # Compute framewise displacement WITHOUT fMRIPrep's own formulation ("noz")
    # First-difference in time of each parameter
    motion_diff = np.diff(confounds_par, axis=0, prepend=np.zeros((1, confounds_par.shape[1])))

    # Weight rotations by 50mm, translations left as mm
    weighted_motion = np.concatenate(
        [
            50 * np.abs(motion_diff[:, :3]),   # rotations (rad → mm with 50mm radius)
            np.abs(motion_diff[:, 3:])        # translations (mm)
        ],
        axis=1
    )

    FD_noz = np.sum(weighted_motion, axis=1)

    # Store your custom FD in the dataframe
    confounds["framewise_displacement_noz"] = FD_noz

    # ----- FD metrics (original fMRIPrep + your "noz") -----
    fd = confounds["framewise_displacement"]
    fd_noz = confounds["framewise_displacement_noz"]

    mean_fd = fd.mean()
    mean_fd_noz = fd_noz.mean()

    p80_fd = fd.quantile(0.80)
    p80_fd_noz = fd_noz.quantile(0.80)

    max_fd = fd.max()
    max_fd_noz = fd_noz.max()

    # ----- motion extrema -----
    max_trans_x = confounds["trans_x"].abs().max()
    max_trans_y = confounds["trans_y"].abs().max()
    max_trans_z = confounds["trans_z"].abs().max()

    # convert rotations to degrees for easier interpretation in the CSV
    confounds[rot_cols] = confounds[rot_cols] * (180 / np.pi)

    max_rot_x = confounds["rot_x"].abs().max()
    max_rot_y = confounds["rot_y"].abs().max()
    max_rot_z = confounds["rot_z"].abs().max()

    # ----- suprathreshold FD proportion (in %) -----
    supra_thresh = confounds[fd_noz > FD_THRESHOLD]
    percent_suprathreshold = 100 * len(supra_thresh) / len(confounds)

    # ----- collect row -----
    rows.append([
        subject,
        mean_fd,
        mean_fd_noz,
        p80_fd,
        p80_fd_noz,
        max_fd,
        max_fd_noz,
        percent_suprathreshold,
        max_trans_x,
        max_trans_y,
        max_trans_z,
        max_rot_x,
        max_rot_y,
        max_rot_z
    ])

# Save exclusion metrics
out_df = pd.DataFrame(rows, columns=columns)
out_df.to_csv(os.path.join(os.getcwd(), "exclusion.csv"), index=False)
print("Saved exclusion.csv")
