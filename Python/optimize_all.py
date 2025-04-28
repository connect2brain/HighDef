import os
from os import system as run_command
import pynibs
from datetime import datetime, timedelta
import time
import argparse

GLOBAL_START = time.time()

ROOT = "/home/bnplab-admin/TMS_localization/HighDef"

parser = argparse.ArgumentParser(description='Projects results for all specified subjects to fsavg')
parser.add_argument('-s', '--subject', help='id of subject, "all", or comma-separated list (e.g. -s "sub-001, sub-002, sub-005")', required=False, type=str, default="all")
args = parser.parse_args()


if args.subject == "all":
    subjects = sorted([s for s in os.listdir(ROOT) if s.startswith("sub") and not os.path.isfile(s)])
else:
    subjects = [s.strip() for s in args.subject.split(",")]

mesh_id = "mesh0"
n_cpu = 28

for subject_id in subjects:
    subject_folder = f"{ROOT}/{subject_id}"
    fn_subject = f"{ROOT}/{subject_id}/{subject_id}.hdf5"
    subject = pynibs.load_subject(fn_subject)
    for exp_id, exp_specification in subject.exp.items():
        hemisphere = exp_id.split('-')[1][0]
        roi_id = f'midlayer_{hemisphere.lower()}'
        for response in exp_specification["response_columns"]:
            print(f' ^  > {response}\t\t ({subject_id}/{mesh_id}/{roi_id}:{exp_id})')
            cmd = f'python TMS_localization/08_calc_opt_coil_pos.py -s {fn_subject} -m {mesh_id} -e {exp_id} -n {n_cpu} -l {response} -q "E_mag" -t {subject_folder}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{response}/sigmoid4/r2_roi_data.hdf5 -a "scalar" -c "/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd"'
            start = time.time()
            os.system(cmd)
            print(f'    < {response}\t\t ({subject_id}/{mesh_id}/{roi_id}:{exp_id}): Took {timedelta(seconds=time.time() - start)}\n')

