# In this script:
# (1) Run create_subject again for each subject
# (2) Rename opt to "opt-pre-review-<timedate>"
#     and rename results/<exp_id>/r2 to results/<exp_id>/r2-pre-review-<timedate>
# (3) Run the regression step (07_custom)
# (4) Run the optimal coil position (08_custom)
# (5) Project back to fs_avg
# (6) Run cross_compare_all_spots.py

# (x7) Dont yet do FWHM rerun. Do this later!

import os
from os import system as run_command
import pynibs
from datetime import datetime, timedelta
import time
import logging

import pandas as pd
import h5py

GLOBAL_START = time.time()

ROOT = "/home/bnplab-admin/TMS_localization/HighDef"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler(f"{ROOT}/address-review-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}.log", encoding='utf-8')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)

subjects = sorted([s for s in os.listdir(ROOT) if s.startswith("sub") and not os.path.isfile(s)])
subjects.reverse()
logger.info(f'Subjects: {subjects}')
print(f'Subjects: {subjects}')

n_cpu = 30
mesh_id = "mesh0"


for subject_id in subjects:
    subject_start = time.time()
    
    # [0:i]   Rerun create_subject
    logger.info(f'Running create_{subject_id}.py')
    run_command(f'python {ROOT}/{subject_id}/create_{subject_id}.py')
    logger.info(f'Completed running create_{subject_id}.py')
    
    subject_folder = f"{ROOT}/{subject_id}"
    fn_subject = f"{subject_folder}/{subject_id}.hdf5"
    subject = pynibs.load_subject(fn_subject)

    # [0:ii]  Check that the subject is complete (i.e. all map-R, map-L, map-R2L, map-L2R are present in the subject-hdf5) -- otherwise skip    
    if not all([e in subject.exp.keys() for e in ["map-L", "map-R", "map-R2L", "map-L2R"]]):
        print(f">>] {subject_id} incomplete: Only found {subject.exp.keys()}, but need all 'map-L', 'map-R', 'map-L2R', 'map-R2L' --- thus skipping this subject")
        logger.warning(f'{subject_id} is incomplete: SKIPPING')
        continue

    
    # [2.a] Rename opt:
    if os.path.exists(f"{ROOT}/{subject_id}/opt"):
        os.rename(f"{ROOT}/{subject_id}/opt", f"{ROOT}/{subject_id}/opt-pre-review-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}")
    
    os.mkdir(f"{ROOT}/{subject_id}/opt")
    

    dont_redo = ["map-L", "map-R"] # only used for hotspot in subsequent session AND: the result has already been used. There is no purpose in rerunning these.
    for exp_id, exp_specification in subject.exp.items():
        if exp_id in dont_redo:
            continue

        hemisphere = exp_id.split('-')[1][0]
        print(f'Experiment: {exp_id}\nHemisphere: {hemisphere}')
        roi_id = f'midlayer_{hemisphere.lower()}'
        
        # [2.b] Rename results/<exp_id>/r2
        r2_path = f"{ROOT}/{subject_id}/results/exp_{exp_id}/r2"
        if os.path.exists(r2_path):
            os.rename(r2_path, f"{r2_path}-pre-review-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}")
        os.mkdir(r2_path)

        # CHECK if number of rows in new csv file matches existing n_trials in E-field, else recompute e-field
        recompute_e = False

        # Load csv data:
        if len(exp_id) > len("map-L2R"): # That is: map-L2R-plus-pi
            raw_data = pd.read_csv(f'{ROOT}/{subject_id}/{subject_id}_{exp_id}_raw.csv')
        else:
            raw_data = pd.read_csv(f'{ROOT}/{subject_id}/{subject_id}_{exp_id.split("-")[1]}_raw.csv')
        
        n_trials_csv = len(raw_data)

        # Load E-field data:
        fn_e_roi = f"{ROOT}/{subject_id}/results/exp_{exp_id}/electric_field/mesh_{mesh_id}/roi_{roi_id}/e.hdf5"
        with h5py.File(fn_e_roi, "r") as f:
            # Here: Neuron-model/Layer_id stuff was removed, bc. not needed for now.
            e_matrix = f["E_mag"][:]

        n_trials_E = e_matrix.shape[0]

        recompute_e = n_trials_csv != n_trials_E
        if recompute_e:
            logger.info(f"subject={subject_id}, exp={exp_id}: Csv file has {n_trials_csv} trials, but E-matrix has {n_trials_E}")
            logger.info(f"Thus recompute E!")
            e_path = f"{ROOT}/{subject_id}/results/exp_{exp_id}/electric_field"
            if os.path.exists(e_path):
                os.rename(e_path, f"{e_path}-pre-review-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}")
            os.mkdir(e_path)
        
            # (2/a) E-field simulation
            start = time.time()
            cmd = f'python /home/bnplab-admin/TMS_localization/0506_custom.py -f {fn_subject} -e {exp_id} -m {mesh_id} -r {roi_id} -n {n_cpu}'
            logger.info(f"E-field simulation step for {subject_id} exp={exp_id} roi={roi_id}:\t\t{cmd}")
            print(f'\nRunning E-field simulation for exp={exp_id}, mesh={mesh_id}, roi={roi_id}')
            run_command(cmd)
            print(f'\nDone running E-field simulation: Took {timedelta(seconds=time.time() - start)}\n')
            logger.info(f'Completed E-field simulation step for {subject_id}, took {timedelta(seconds=time.time() - start)}')

        # [3] Run regression
        cmd = f'python /home/bnplab-admin/TMS_localization/07_custom.py -f {fn_subject} -e {exp_id} -m {mesh_id} -r {roi_id} -n {n_cpu} --splits'
        logger.info(f"Regression step for {subject_id} exp={exp_id} roi={roi_id}:\t\t{cmd}")
        start = time.time()
        run_command(cmd)
        logger.info(f'Completed Regression step for {subject_id}, took {timedelta(seconds=time.time() - start)}')
        print(f"\n\n\n   Regression done for {subject_id} exp={exp_id} roi={roi_id}\n")



        # [4] Run optimization
        for response in exp_specification["response_columns"]:
            print(f' ^  > {response}\t\t ({subject_id}/{mesh_id}/{roi_id}:{exp_id})')
            cmd = f'python /home/bnplab-admin/TMS_localization/08_calc_opt_coil_pos.py -s {fn_subject} -m {mesh_id} -e {exp_id} -n {n_cpu} -l {response} -q "E_mag" -t {subject_folder}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{response}/sigmoid4/r2_roi_data.hdf5 -a "scalar" -c "/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd"'
            logger.info(f"Optimization step for {subject_id} exp={exp_id} roi={roi_id} response={response}:\t\t{cmd}")
            start = time.time()
            run_command(cmd)
            logger.info(f'Completed Optimization step for {subject_id} exp={exp_id} roi={roi_id} response={response}, took {timedelta(seconds=time.time() - start)}')

        print(f"\n\n\n   Optimization done for {subject_id} {exp_id}\n")


        
        # [5] Project back to fs_avg
        logger.info(f'Projecting back to avg. for {subject_id} exp={exp_id} roi={roi_id}')
        start = time.time()
        print(f'\nProjecting {subject_id} {exp_id} back to average brain')
        run_command(f"python /home/bnplab-admin/TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o CsE_FDI_in_uV -m mesh0")
        run_command(f"python /home/bnplab-admin/TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o SIHIscore_FDI -m mesh0")
        
        run_command(f"python /home/bnplab-admin/TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o CsE_FDI_in_uV -m mesh0 -p first_half")
        run_command(f"python /home/bnplab-admin/TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o SIHIscore_FDI -m mesh0 -p first_half")

        run_command(f"python /home/bnplab-admin/TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o CsE_FDI_in_uV -m mesh0 -p second_half")
        run_command(f"python /home/bnplab-admin/TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o SIHIscore_FDI -m mesh0 -p second_half")
        print(f'\nDone projecting {subject_id} {exp_id} back to average: Took {timedelta(seconds=time.time() - start)}\n\n')
        logger.info(f'Completed projecting back to avg. for {subject_id} exp={exp_id} roi={roi_id}, took {timedelta(seconds=time.time() - start)}')

    logger.info(f'Completed {subject_id} in {timedelta(seconds=time.time() - subject_start)}')
    print(f"\n\n\n\n < {subject_id} D O N E \n < {timedelta(seconds=time.time() - subject_start)} \n\n\n\n")


# [6] Do cross-comparison
print("Cross-comparison")
cmd = "python /home/bnplab-admin/TMS_localization/cross_compare_spots.py"
start = time.time()
logger.info(f'Cross-comparing spots:\t\t{cmd}')
run_command(cmd)
logger.info(f'Completed ross-comparing spots, took {timedelta(seconds=time.time() - start)}')



print(f"\n\nDone. Took {timedelta(seconds=time.time() - GLOBAL_START)}.")




