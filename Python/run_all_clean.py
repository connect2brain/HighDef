# TODO: In this script, run all the evaluation and analysis steps on all complete subjects
# (1) Create headmodel [optional: only if not yet done] and Regions of interest [obligatory!]
# (2) For all experiments:
#   (2.1) E field simulation
#   (2.2) Regression
#   (2.3) Optimal coil position
# (3) Project FDI hotspot and coldspot onto average brain, for map-R2L and map-L2R
# (4) Run FWHM analysis for FDI hotspot and coldspot (for map-R2L and map-L2R)

import os
from os import system as run_command
import pynibs
from datetime import datetime, timedelta
import time
import logging

GLOBAL_START = time.time()

ROOT = "/home/bnplab-admin/TMS_localization/HighDef"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler(f"{ROOT}/clean_run-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}.log", encoding='utf-8')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)

subjects = ["sub-006", "sub-004", "sub-003"] #sorted([s for s in os.listdir(ROOT) if s.startswith("sub") and not os.path.isfile(s)])
logger.info(f'Subjects: {subjects}')

for subject_id in subjects:
    # [0:i]   Rerun create_subject
    logger.info(f'Running create_{subject_id}.py')
    run_command(f'python {ROOT}/{subject_id}/create_{subject_id}.py')
    logger.info(f'Completed running create_{subject_id}.py')
    fn_subject = f"{ROOT}/{subject_id}/{subject_id}.hdf5"
    subject = pynibs.load_subject(fn_subject)

    # [0:ii]  Check that the subject is complete (i.e. all map-R, map-L, map-R2L, map-L2R are present in the subject-hdf5) -- otherwise skip    
    if not all([e in subject.exp.keys() for e in ["map-L", "map-R", "map-R2L", "map-L2R"]]):
        print(f">>] {subject_id} incomplete: Only found {subject.exp.keys()}, but need all 'map-L', 'map-R', 'map-L2R', 'map-R2L' --- thus skipping this subject")
        logger.warning(f'{subject_id} is incomplete: SKIPPING')
        continue

    # [0:iii] Reset ROIs!
    logger.info(f'Resetting ROI folder for {subject_id}')
    if os.path.exists(f"{ROOT}/{subject_id}/mesh/roi"):
        os.rename(f"{ROOT}/{subject_id}/mesh/roi", f"{ROOT}/{subject_id}/mesh/roi-outdated-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}")
    os.mkdir(f"{ROOT}/{subject_id}/mesh/roi")  
    logger.info(f'Completed resetting ROI folder for {subject_id}')
    
    # [1] Run headmodel script
    logger.info(f'Headmodel/ROI construction for {subject_id}')
    run_command(f"python Installation/tmsloc_proto-0.2023.8/scripts/02_make_msh_from_mri_simnibs4.py -s {fn_subject} -m mesh0 --charm_ini /mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/charm.ini --charm_cmd_params forceqform")
    logger.info(f'Completed headmodel/ROI construction for {subject_id}')

    # [2] Run evaluate.py with correct ROI
    
    logger.info(f'Full evaluation for {subject_id}')
    # Reset exp, results, opt
    if os.path.exists(f"{ROOT}/{subject_id}/exp"):
        os.rename(f"{ROOT}/{subject_id}/exp", f"{ROOT}/{subject_id}/exp-outdated-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}") 
    if os.path.exists(f"{ROOT}/{subject_id}/results"):
        os.rename(f"{ROOT}/{subject_id}/results", f"{ROOT}/{subject_id}/results-outdated-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}")
    if os.path.exists(f"{ROOT}/{subject_id}/opt"):
        os.rename(f"{ROOT}/{subject_id}/opt", f"{ROOT}/{subject_id}/opt-outdated-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}")
    os.mkdir(f"{ROOT}/{subject_id}/exp")
    os.mkdir(f"{ROOT}/{subject_id}/results")
    os.mkdir(f"{ROOT}/{subject_id}/opt")
    start = time.time()
    print(f'\nFull evaluation for {subject_id}')
    run_command(f'python /home/bnplab-admin/TMS_localization/evaluate.py -f {fn_subject}')
    print(f'\n\n    Done with full evaluation of {subject_id}:\n    Took {timedelta(seconds=time.time() - start)}\n\n\n\n\n')
    logger.info(f'Completed full evaluation for {subject_id}, took {timedelta(seconds=time.time() - start)}')

    for i_exp, exp_id in enumerate(["map-L2R", "map-R2L"]):
        roi_id = f'midlayer_{exp_id.split("-")[-1].split("2")[0].lower()}'
        logger.info(f'Post-evaluation for {subject_id} exp={exp_id} roi={roi_id}')
        
        # [3] Project to average headmodel (only for map-R2L and map-L2R)
        logger.info(f'Projecting back to avg. for {subject_id} exp={exp_id} roi={roi_id}')
        start = time.time()
        print(f'\nProjecting {subject_id} {exp_id} back to average brain')
        run_command(f"python TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o CsE_FDI_in_uV -m mesh0")
        run_command(f"python TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o SIHIscore_FDI -m mesh0")
        print(f'\nDone projecting {subject_id} {exp_id} back to average: Took {timedelta(seconds=time.time() - start)}\n\n')
        logger.info(f'Completed projecting back to avg. for {subject_id} exp={exp_id} roi={roi_id}, took {timedelta(seconds=time.time() - start)}')

        # [4] Run FWHM
        start = time.time()
        print(f'\nCollecting FWHM for {subject_id} {exp_id}')
        logger.info(f'Collecting SIHIscore_FDI FWHM for {subject_id} exp={exp_id} roi={roi_id}')
        run_command(f"python TMS_localization/simulate_FWHM_focality.py -s {subject_id} -e {exp_id} -o SIHIscore_FDI")
        logger.info(f'Completed collecting SIHIscore_FDI FWHM for {subject_id} exp={exp_id} roi={roi_id}')

        logger.info(f'Collecting CsE_FDI_in_uV FWHM for {subject_id} exp={exp_id} roi={roi_id}')
        run_command(f"python TMS_localization/simulate_FWHM_focality.py -s {subject_id} -e {exp_id} -o CsE_FDI_in_uV")
        logger.info(f'Completed collecting CsE_FDI_in_uV FWHM for {subject_id} exp={exp_id} roi={roi_id}')
        print(f'\nDone collecting FWHM for {subject_id} {exp_id}: Took {timedelta(seconds=time.time() - start)}\n\n')

    logger.info(f'Completed {subject_id}')

print(f"\n\nDone. Took {timedelta(seconds=time.time() - GLOBAL_START)}.")




