from glob import glob
import pynibs
import os
import time
from datetime import timedelta

n_cpu = 25
mesh_id = 'mesh0'
# Script 08 already skips a big part of the work load if the optimization has already been done (i.e.: data are present)

# 1  Get list of all subjects (name format: sub-...)
# 2  For each subject
#    3  Get the configured experiments
#    4  For each configured experiment:
#       5  For each configured Response Column:
#       6  Run script 8 from command-line

for fn_subject in glob('/home/bnplab-admin/TMS_localization/HighDef/sub-???/sub-???.hdf5'):
    subject = pynibs.load_subject(fn_subject)
    subject_folder, subject_id = os.path.split(fn_subject)
    subject_id = subject_id.split('.hdf5')[0]
    print(f'All optimization for subject: {subject_id}')
    for exp_id, exp_specification in subject.exp.items():
        print(f' > {subject_id}: exp_id={exp_id}')
        data_file = exp_specification["response_csv"]
        hemisphere = exp_id.split('-')[1][0] # map-R2L/map-R -> ['map', 'R...'] -> 'R...' -> 'R'
        roi_id = f'midlayer_{hemisphere.lower()}'

        for response in exp_specification["response_columns"]:
            print(f' ^  > {response}\t\t ({subject_id}/{mesh_id}/{roi_id}:{exp_id})')
            cmd = f'python TMS_localization/08_calc_opt_coil_pos.py -s {fn_subject} -m {mesh_id} -e {exp_id} -n {n_cpu} -l {response} -q "E_mag" -t {subject_folder}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{response}/sigmoid4/r2_roi_data.hdf5 -a "scalar" -c "/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd"'
            start = time.time()
            os.system(cmd)
            print(f'    < {response}\t\t ({subject_id}/{mesh_id}/{roi_id}:{exp_id}): Took {timedelta(seconds=time.time() - start)}\n')


