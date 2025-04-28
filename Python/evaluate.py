# This script is supposed to do evaluation (0506, 07, 08) for all available experiments of the given subject
import os
import time
from datetime import timedelta
import argparse
import pynibs

sep = '_'.join(40 * ['_'])
lightsep = ' '.join(40*['^'])


parser = argparse.ArgumentParser(description='Evaluates all experiments registered for given subject (also reruns create_subject)')
parser.add_argument('-f', '--fn_subject', help='path of subject', required=True, type=str)
parser.add_argument('-m', '--mesh_id', help='id of the mesh', default='mesh0', type=str)
parser.add_argument('-r', '--roi_id', help='id of the ROI', type=str, required=False)
parser.add_argument('-e', '--exp_id', help='id of the experiment', default='all', type=str)
parser.add_argument('-n', '--n_cpu', help='number of cpu cores to use', default=30, type=int)
args = parser.parse_args()

fn_subject = os.path.abspath(args.fn_subject)
subject_folder, subject_id = os.path.split(fn_subject)
subject_id = subject_id.split('.hdf5')[0]

mesh_id = args.mesh_id # global assumption
target_exp = args.exp_id

n_cpu = args.n_cpu


# (1) Run create_subject.py
os.system(f'python {subject_folder}/create_{subject_id}.py')

# (2) Load subject.hdf5 and loop over all experiments
subject = pynibs.load_subject(fn_subject)

for exp_id, exp_specification in subject.exp.items():
    if target_exp != "all" and exp_id != target_exp:
        continue
    print(f'\n{sep}\n\n Experiment: {exp_id}')
    #data_file = exp_specification["response_csv"]
    hemisphere = exp_id.split('-')[1][0] # map-R2L/map-R -> ['map', 'R...'] -> 'R...' -> 'R'
    if args.roi_id is None:
        roi_id = f'midlayer_{hemisphere.lower()}'
    else:
        roi_id = args.roi_id
    # (2/a) E-field simulation
    start = time.time()
    print(f'\nRunning E-field simulation for exp={exp_id}, mesh={mesh_id}, roi={roi_id}')
    os.system(f'python TMS_localization/0506_custom.py -f {fn_subject} -e {exp_id} -m {mesh_id} -r {roi_id} -n {n_cpu}')
    print(f'\nDone running E-field simulation: Took {timedelta(seconds=time.time() - start)}\n{lightsep}\n')

    # (2/b) Regression
    start = time.time()
    print(f'\nRunning Regression for exp={exp_id}, mesh={mesh_id}, roi={roi_id}')
    os.system(f'python TMS_localization/07_custom.py -f {fn_subject} -e {exp_id} -m {mesh_id} -r {roi_id} -n {n_cpu}')
    print(f'\nDone running Regression: Took {timedelta(seconds=time.time() - start)}\n\n\n')

    # (2/c) Find optimal coil location
    for response in exp_specification["response_columns"]:
        print(f' ^  > {response}\t\t ({subject_id}/{mesh_id}/{roi_id}:{exp_id})')
        cmd = f'python TMS_localization/08_calc_opt_coil_pos.py -s {fn_subject} -m {mesh_id} -e {exp_id} -n {n_cpu} -l {response} -q "E_mag" -t {subject_folder}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{response}/sigmoid4/r2_roi_data.hdf5 -a "scalar" -c "/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd"'
        start = time.time()
        os.system(cmd)
        print(f'    < {response}\t\t ({subject_id}/{mesh_id}/{roi_id}:{exp_id}): Took {timedelta(seconds=time.time() - start)}\n')


